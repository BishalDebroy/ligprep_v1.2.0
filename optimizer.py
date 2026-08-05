"""
optimizer.py - Optimize a single conformer with MMFF94/UFF and retrieve energy.
"""

import logging
from typing import Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

logger = logging.getLogger("ligprep")

def optimize_and_energy(
        mol: Chem.Mol,
        force_field: str = "MMFF94",
) -> Tuple[Optional[Chem.Mol], Optional[float], str]:
    """
    optimize mol (in place, return copy) with MMFF94 (fallback UFF).
    returns:
    - optimized mol, energy (kcal/mol), used_force_field
    Energy is None if optimisation failed.
    """
    mol = Chem.RWMol(mol)
    used_ff = ""
    energy = None

    try:
        if force_field.upper() == "MMFF94":
            if AllChem.MMFFHasAllMoleculeParams(mol):
                ff = AllChem.MMFFGetMoleculeForceField(mol)
                if ff:
                    ff.Minimize()
                    energy = ff.CalcEnergy()
                    used_ff = "MMFF94"
    except Exception as e:
        logger.debug(f"MMFF94 failed: {e}")

    if not used_ff:
        try:
            ff = AllChem.UFFGetMoleculeForceField(mol)
            ff.Minimize()
            energy = ff.CalcEnergy()
            used_ff = "UFF"
        except Exception as e:
            logger.warning(f"Uff optimization failed: {e}")
            return mol.GetMol(), None, ""

    return mol.GetMol(), energy, used_ff