"""
reader.py - iterate over SMILES and SDF files with minimal memory footprient.
"""

import csv
import io
import re
from pathlib import Path
from typing import Generator, Tuple

def detect_delimiter(first_line: str) -> str:
    """Return the delimiter that yields exactly two fields: tab, comma, or whitespace."""
    line = first_line.strip()
    if not line:
        return "\t"
    for delim in ["\t", ",", " "]:
        parts = line.split(delim)
        if len(parts) == 2:
            return delim

    # Falback to whitespace spliting (Collapse multiple space)
    parts = line.split()
    if len(parts) >= 2:
        return " "
    raise ValueError(f"Cannot determine delimiter from line: {line}")

def read_smiles(file_path: str) -> Generator[Tuple[str, str, int], None, None]:
    """
    Yield (smiles, name, line number) from a SMILE file.
    Handle comment (#) and blanket lines.
    Automatically detects delimiter
    """
    path = Path(file_path)
    with path.open("r", encoding="utf-8") as fh:
        # Skip blank/comment lines to find first data row
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            delimiter = detect_delimiter(line)
            break
        else:
            return # empty file

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