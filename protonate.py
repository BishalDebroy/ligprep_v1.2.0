"""
protonate.py - Generate major protonation states at a given pH using Dimorphite-DL.
"""

import logging
from typing import Optional, List

import dimorphite_dl
from rdkit import Chem
from rdkit.Chem import AllChem

logger = logging.getLogger("ligprep")

def protonate(
        mol: Chem.Mol,
        pH: float = 7.4,
        min_ph: Optional[float] = None,
        max_ph: Optional[float] = None,
        max_states: int = 1,
) -> Optional[Chem.Mol]:
    """
    Return the best protonated mol at the specified pH (or range).
    uses Dimorphite-DL; returns None if no valid state found.
    """
    if min_ph is None:
        min_ph = pH
    if max_ph is None:
        max_ph = pH

    # Dimorphite-DL operates on SMILES strings, not RDKit Mol objects, and
    # expects a pH range; set equal values for a single pH.
    try:
        smi = Chem.MolToSmiles(Chem.RemoveHs(mol))
        variants = dimorphite_dl.protonate_smiles(
            smi,
            ph_min=min_ph,
            ph_max=max_ph,
            max_variants=max_states,
        )
    except Exception as e:
        logger.error(f"Dimorphite-DL error: {e}")
        return None

    if not variants:
        # Fallback: add hydrogens at the input pH using RDKit only
        return Chem.AddHs(mol)

    # Take the first (dominant) state and rebuild a 3D, explicit-H mol,
    # since protonation via SMILES loses the original conformer.
    best_mol = Chem.MolFromSmiles(variants[0])
    if best_mol is None:
        logger.warning(f"Could not parse protonated SMILES: {variants[0]}")
        return None

    best_mol = Chem.AddHs(best_mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(best_mol, params)
    if status != 0:
        logger.warning(f"Embedding failed for protonated SMILES: {variants[0]}")
        return None

    return best_mol