"""
tautomer.py - Generate and rank tautomers using RDKit TautomerEnumerator.
"""
import logging
from typing import List, Tuple, Optional

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

logger = logging.getLogger("ligprep")

def enumerate_tautomers(mol: Chem.Mol) -> List[Tuple[Chem.Mol, bool]]:
    """
    Return list of (tautomer_mol, is_canonical) sorted with canonical first.
    """
    enumerator = rdMolStandardize.TautomerEnumerator()
    # Remove any existing hydrogen, enumerator adds its own.
    mol = Chem.RemoveHs(mol)
    result = enumerator.Enumerate(mol)
    if not result or not result.tautomers:
        logger.warning("No tautomers enumerated, returning input mol as canonical.")
        return[(mol, True)]

    tautomers = list(result.tautomers)
    canonical_mol = enumerator.Canonicalize(mol)
    canonical_smi = Chem.MolToSmiles(canonical_mol, canonical=True)
    # Create list with canonical flag
    out = []
    for tmol in tautomers:
        is_canon = Chem.MolToSmiles(tmol, canonical=True) == canonical_smi
        out.append((tmol, is_canon))

    out.sort(key=lambda x: not x[1]) # True (is Canonical) first
    return out