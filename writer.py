"""
writer.py - Write mol blocks to an SDF (append mode, thread/process safe).
"""
import io
from typing import Dict
from pathlib import Path
from rdkit import Chem

def mol_with_tags_to_sdf_block(mol: Chem.Mol, tags: Dict[str, str]) -> str:
    """Convert an RDKit Mol and a dict of tags into a single SDF block string.

    The SDF title line (line 1 of each block) is set to the Tautomer_Name
    if present (second SDF), otherwise to the Compound_Name (first SDF).
    """
    title = tags.get("Tautomer_Name") or tags.get("Compound_Name") or ""
    mol.SetProp("_Name", str(title))

    # set all properties
    for key, value in tags.items():
        mol.SetProp(key, str(value))
    sio = io.StringIO()
    w = Chem.SDWriter(sio)
    w.write(mol)
    w.flush()
    w.close()
    block = sio.getvalue()
    return block

def write_sdf_blocks(blocks, output_path: str, append: bool = False):
    """
    Write an iterable of SDF block strings to a file.
    Each block mught end with '$$$$\n'.
    """
    mode = "a" if append else "w"
    with open(output_path, mode, encoding="utf-8") as fh:
        for blk in blocks:
            # Ensure block ends with $$$$
            if not blk.strip().endswith("$$$$"):
                blk = blk.rstrip() + "\n$$$$\n"
            fh.write(blk)