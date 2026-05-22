#!/usr/bin/env python3

import re
from csv import DictWriter
from itertools import groupby
from pathlib import Path
from typing import Iterable

file = "/home/shannc/Downloads/files.txt"
root = Path("/volume1/cancer_ngs")

files = [
    p
    for p in (
        Path(f).relative_to(root) for f in Path(file).read_text().strip().split("\n")
    )
    if p.suffix
]


def get_group_key(path: Path):
    parents = list(path.parents)
    if parents and len(parents) >= 2:
        return parents[-2]
    return Path("/")


def group_iterator(parent: Path, parent_contents: Iterable[Path]):
    return groupby([p.relative_to(parent) for p in parent_contents], key=get_group_key)


RAW_EXTS = [".fastq.gz", ".fasta", ".fq.gz"]

NORMAL_RE: str = f".*_B_[12].({'|'.join(RAW_EXTS)})"
TUMOR_RE: str = f".*_C_[12].({'|'.join(RAW_EXTS)})"


def _check_raw_type(file: Path, regexp: str):
    for ext in RAW_EXTS:
        if file.name.endswith(ext):
            removed = file.name.removesuffix(ext)
            if re.match(regexp, removed):
                return True
            return False
    return False


def raw_is_normal(file: Path) -> bool:
    return _check_raw_type(file, ".*_B(_[0-9]+)?$")


def raw_is_tumor(file: Path) -> bool:
    return _check_raw_type(file, ".*[_-]T?C(_[0-9]+)?$")


def check_raw_dir(files: Iterable[Path]) -> tuple[bool, bool, int]:
    """
    Check the raw files `files` in a sample directory,
    returning a tuple of `has_pbmc` (normals), `has_tumor_data`,
        and the number raw tumor biopsies available
    """
    has_tumor, has_normal = False, False
    n_tumors = 0
    for file in files:
        if raw_is_normal(file):
            has_normal = True
    return has_normal, has_tumor, n_tumors


def catalog_sample(
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
        if dir.name == "raw":
            result["has_raw"] = True
            has_n, has_t, n_biopsies = check_raw_dir(dir_contents)
            result["has_pbmc"] = has_n
            result["has_tumor"] = has_t
            result["n_biopsies"] = n_biopsies
            if not has_n and not has_t:
                result["warnings"].append(
                    "nonempty raw directory without PBMC or tumor files"
                )
        elif dir.name == "processed":
            result["has_processed"] = True
            for file in dir_contents:
                result["processed_files"].append(file.name)
    result["warnings"] = ";".join(result["warnings"])
    result["processed_files"] = ";".join(result["processed_files"])
    return result


rows = []

grouped = groupby(files, key=get_group_key)
blacklist = ("summary", "/")
ttypes = ("pdac", "hcc", "cca", "brca")
modalities = ("exome", "scrna_seq", "tcr_seq", "sc_atac_seq", "rna_seq")
tmp = []
for modality, m_contents in grouped:
    if modality.name not in modalities:
        continue
    for ttype, t_contents in group_iterator(modality, m_contents):
        for cohort, c_contents in group_iterator(ttype, t_contents):
            for sample, s_contents in group_iterator(cohort, c_contents):
                if sample.stem in blacklist:
                    continue
                tmp.append((sample, list(s_contents)))
                rows.append(
                    catalog_sample(
                        sample_dir=sample,
                        contents=s_contents,
                        modality=modality.name,
                        ttype=ttype.name,
                        cohort=cohort.name,
                    )
                )

try:
    import pytest

    @pytest.mark.parametrize(
        "raw_type,fname,expectation",
        [
            ("n", "A_B.fq.gz", True),
            ("n", "s13-B.fasta", False),
            ("n", "s24_B.fasta", True),
            ("n", "s24_B_12.fasta", True),
            ("n", "s24_B_2.fasta", True),
            ("n", "s24_CB_2.fasta", False),
            ("t", "s24_C.fasta", True),
            ("t", "s24_C_19.fasta", True),
            ("t", "s24_TC_2.fasta", True),
            ("t", "s24_C_9.fq.gz", True),
            ("t", "s24_C-9.fasta", False),
            ("t", "s24_9-C_1.fq.gz", True),
        ],
    )
    def test_raw_checks(raw_type, fname, expectation):
        if raw_type == "n":
            assert raw_is_normal(Path(fname)) == expectation
        else:
            assert raw_is_tumor(Path(fname)) == expectation

except ImportError:
    pass
