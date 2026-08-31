"""
reader.py - iterate over SMILES and SDF files with minimal memory footprint.
"""

import csv
import io
import re
from pathlib import Path
from typing import Generator, Tuple, Optional

def detect_delimiter(line: str) -> Optional[str]:
    """
    Try to determine the delimiter (tab, comma, or whitespace) from a single
    line, based on which one splits it into exactly two fields.

    Returns None if the line is inconclusive - e.g. a bare SMILES with no
    name/ID column - rather than raising, since a single such line
    elsewhere in a large file shouldn't abort reading the whole file.
    """
    line = line.strip()
    if not line:
        return None
    for delim in ["\t", ",", " "]:
        parts = line.split(delim)
        if len(parts) == 2:
            return delim
    return None

def read_smiles(file_path: str) -> Generator[Tuple[str, str, int], None, None]:
    """
    Yield (smiles, name, line number) from a SMILE file.
    Handle comment (#) and blank lines.
    Automatically detects delimiter by scanning for the first line that
    unambiguously reveals it (exactly two fields); lines with only a bare
    SMILES and no name are tolerated anywhere in the file.
    """
    path = Path(file_path)
    delimiter = None
    found_any_data_line = False
    with path.open("r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            found_any_data_line = True
            delimiter = detect_delimiter(line)
            if delimiter is not None:
                break
            # Inconclusive line (e.g. a bare SMILES, no name column) -
            # keep scanning later lines for one that reveals the delimiter.

    if not found_any_data_line:
        return  # empty file

    if delimiter is None:
        # No line in the file yielded exactly two fields - the file is
        # likely entirely single-column SMILES (no name/ID column at all).
        # Whitespace is a safe default: such lines have nothing to split on.
        delimiter = " "

    with path.open("r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        line_number = 0
        for row in reader:
            line_number += 1
            if not row or row[0].strip().startswith("#"):
                continue
            if len(row) < 2:
                # only SMILE, name missing → skip later
                yield row[0].strip(), "", line_number
            else:
                yield row[0].strip(), row[1].strip(), line_number

def read_sdf_blocks(file_path: str) -> Generator[Tuple[str, int], None, None]:
    """
    Yield (mol_block_strings, entry_index) from a multiple-ligand SDF.
    The file is read line by line; entries are separated by '$$$'.
    """
    path = Path(file_path)
    with path.open("r", encoding="utf-8") as fh:
        block_lines = []
        entry_idx = 0
        for line in fh:
            block_lines.append(line)
            if line.strip() == "$$$$":
                block_str = "".join(block_lines)
                yield block_str, entry_idx
                entry_idx += 1
                block_lines = []

        # In case file ends without $$$ (Unlikely for valid SDF)
        if block_lines:
            block_str = "".join(block_lines)
            yield block_str, entry_idx
