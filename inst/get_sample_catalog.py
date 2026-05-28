#!/usr/bin/env python3

import argparse
import re
from csv import DictWriter
from itertools import groupby
from pathlib import Path
from typing import Iterable

RAW_EXTS = [".fastq.gz", ".fasta", ".fq.gz", ".tar"]
RAW_REGEXPS = {"normal": ".*_B(_[0-9]+)?$", "tumor": ".*[_-]T?C([_-][0-9]+)?$"}

file = "/home/shannc/Downloads/files.txt"
root = Path("/volume1/cancer_ngs")


def get_group_key(path: Path):
    parents = list(path.parents)
    if parents and len(parents) >= 2:
        return parents[-2]
    return Path("/")


def split_multiplex(grouped_samples: str) -> list[Path]:
    """
    Split a multiplexed sample string of the form
    SAMPLE1,PREFIX{ID1,ID2,ID3},SAMPLE2 into
    SAMPLE1, PREFIXID1, PREFIXID2, PREFIXID3, SAMPLE
    """
    result = []
    cleaned = grouped_samples
    for brace_match in re.findall(r"(\w+\{.*?,.*?\}\w*)", grouped_samples):
        cleaned = cleaned.replace(brace_match, "")
        extracted = re.findall(r"(\w+)\{(.*,.*)\}(\w*)", brace_match)
        if extracted:
            prefix, sample_list, suffix = extracted[0]
            result.extend(
                [Path(f"{prefix}{s}{suffix}") for s in sample_list.split(",")]
            )
    splits = cleaned.split(",")
    for s in splits:
        if s:
            result.append(Path(s))
    return result


def group_iterator(parent: Path, parent_contents: Iterable[Path]):
    return groupby([p.relative_to(parent) for p in parent_contents], key=get_group_key)


def check_raw_type(file: Path, regexps: dict[str, str]) -> dict:
    for ext in RAW_EXTS:
        if file.name.endswith(ext):
            removed = file.name.removesuffix(ext)
            for k, v in regexps.items():
                if re.match(v, removed):
                    return {"type": k, "ext": v, "filename": removed}
            return {}
    return {}


def check_raw_dir(files: Iterable[Path]) -> tuple[bool, bool, int]:
    """
        Check the raw files `files` in a sample directory,
        returning a tuple of `has_pbmc` (normals), `has_tumor_data`,
            and the number raw tumor biopsies available

        Notes
        -----
    The regexps used below are based on the naming conventions outlined in the
     README.txt file at the base of "cancer_ngs"
    """
    has_tumor, has_normal = False, False
    seen_biopsy_count = set()
    n_biopsies = 0
    for file in files:
        raw_type = check_raw_type(file, RAW_REGEXPS)
        if raw_type.get("type") == "normal":
            has_normal = True
        elif raw_type.get("type") == "tumor":
            has_tumor = True
            biopsy = re.findall("([0-9]+)-T?C(_[0-9]+)?$", raw_type["filename"])
            if biopsy:
                seen_biopsy_count.add(biopsy[0])
    if has_tumor:
        n_biopsies = len(seen_biopsy_count) if seen_biopsy_count else 1
    return has_normal, has_tumor, n_biopsies


def catalog_helper(
    sample_dir: Path, contents: Iterable[Path], ttype: str, cohort: str, modality: str
) -> dict:
    """Catalog the data for the sample represented by
    `sample_dir`

    Returns
    -------
    A dictionary for use with csv DictWriter.


    """
    result = {
        "case_name": sample_dir.name,
        "cohort": cohort,
        "tumor_type": ttype,
        "modality": modality,
        "has_pbmc": False,
        "has_tumor": False,
        "has_processed": False,
        "has_raw": False,
        "n_biopsies": 0,
        "processed_files": [],
        "warnings": [],
    }
    for dir, dir_contents in group_iterator(sample_dir, contents):
        dir_contents = list(dir_contents)
        if dir.name == "raw":
            result["has_raw"] = True
            has_n, has_t, n_biopsies = check_raw_dir(dir_contents)
            result["has_pbmc"] = has_n
            result["has_tumor"] = has_t
            result["n_biopsies"] = n_biopsies
            if not has_n and not has_t:
                as_list: list = [str(d.relative_to(dir)) for d in dir_contents]
                result["warnings"].append(
                    f"non-empty raw directory without PBMC or tumor files. Files in raw: {','.join(as_list)}"
                )
        elif dir.name == "processed":
            result["has_processed"] = True
            for file in dir_contents:
                result["processed_files"].append(file.name)
    result["warnings"] = ";".join(result["warnings"])
    result["processed_files"] = ";".join(result["processed_files"])
    return result


def catalog_sample(
    sample_dir: Path, contents: Iterable[Path], ttype: str, cohort: str, modality: str
) -> dict | list[dict]:
    if re.match(r"\w+\{.*,.*\}\w*", sample_dir.name):
        contents = list(contents)
        return [
            catalog_helper(
                s,
                [Path(str(c).replace(str(sample_dir), str(s))) for c in contents],
                ttype,
                cohort,
                modality,
            )
            for s in split_multiplex(sample_dir.name)
        ]

    return catalog_helper(
        sample_dir, contents=contents, ttype=ttype, cohort=cohort, modality=modality
    )


def main(
    file_list: str,
    root: Path,
    modalities=("exome", "scrna_seq", "tcr_seq", "sc_atac_seq", "rna_seq"),
    blacklist=("summary",),
    ttypes=("pdac", "hcc", "cca", "brca"),
) -> list[dict]:
    files: list[Path] = [
        p
        for p in (
            Path(f).relative_to(root)
            for f in Path(file_list).read_text().strip().split("\n")
        )
        if p.suffix
    ]
    rows = []
    grouped = groupby(files, key=get_group_key)
    for modality, m_contents in grouped:
        if modality.name not in modalities:
            continue
        for ttype, t_contents in group_iterator(modality, m_contents):
            for cohort, c_contents in group_iterator(ttype, t_contents):
                for sample, s_contents in group_iterator(cohort, c_contents):
                    if sample.stem in blacklist or not sample.stem:
                        continue
                    rec = catalog_sample(
                        sample_dir=sample,
                        contents=s_contents,
                        modality=modality.name,
                        ttype=ttype.name,
                        cohort=cohort.name,
                    )
                    if isinstance(rec, dict):
                        rows.append(rec)
                    else:
                        rows.extend(rec)
    return rows


try:
    import pytest

    @pytest.mark.parametrize(
        "raw_type,fname",
        [
            ("normal", "A_B.fq.gz"),
            ("", "s13-B.fasta"),
            ("normal", "s24_B.fasta"),
            ("normal", "s24_B_12.fasta"),
            ("normal", "s24_B_2.fasta"),
            ("", "s24_CB_2.fasta"),
            ("tumor", "s24_C.fasta"),
            ("tumor", "s24_C_19.fasta"),
            ("tumor", "s24_TC_2.fasta"),
            ("tumor", "s24_C_9.fq.gz"),
            ("", "s24_C-9.fasta"),
            ("tumor", "s24_9-C_1.fq.gz"),
        ],
    )
    def test_raw_checks(raw_type, fname):
        assert check_raw_type(Path(fname), RAW_REGEXPS).get("type", "") == raw_type

    @pytest.mark.parametrize(
        "val,expected",
        [
            ("P{1,3,9}FOO", {"P1FOO", "P3FOO", "P9FOO"}),
            (
                "PDAC81,P{1,3,9}FOO,PDAC83",
                {"PDAC81", "PDAC83", "P1FOO", "P3FOO", "P9FOO"},
            ),
            ("PDAC81,PDAC83", {"PDAC81", "PDAC83"}),
            ("PDAC81,PDAC83", {"PDAC81", "PDAC83"}),
            ("P{1,3,9}FOO,S{9,2}BAR", {"P1FOO", "P3FOO", "P9FOO", "S9BAR", "S2BAR"}),
            (
                "baz,P{1,3,9}FOO,S{9,2}BAR,bat",
                {"baz", "bat", "P1FOO", "P3FOO", "P9FOO", "S9BAR", "S2BAR"},
            ),
        ],
    )
    def test_multiplex_split(val: str, expected: set[str]):
        split = split_multiplex(val)
        as_set = {str(s) for s in split}
        assert as_set == expected

except ImportError:
    pass


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        help="Input file listing absolute paths to all files in the repository, one file per line",
    )
    parser.add_argument("-r", "--root", help="Root of the repository")
    parser.add_argument("-c", "--config", default=None)
    parser.add_argument("-o", "--output", help="Output CSV file to write to")
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    args = parse_args()
    rows = main(args["input"], args["root"])
    with open(args["output"], "w") as f:
        writer = DictWriter(
            f,
            fieldnames=[
                "case_name",
                "cohort",
                "tumor_type",
                "modality",
                "has_pbmc",
                "has_tumor",
                "has_processed",
                "has_raw",
                "n_biopsies",
                "processed_files",
                "warnings",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
