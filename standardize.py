"""
Standardize.py - standardize SMILE, removes salts, normalize, generate 3D.
"""

import logging
from typing import Tuple, Optional

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

logger = logging.getLogger("ligprep")

def standardize_and_embed(
        smiles: str,
        force_field: str = "MMFF94",
) -> Tuple[Optional[Chem.Mol], str, str]:
    """
    Steps:
    - Remove salts, keep larger fragments
    - Normalize functional groups, aromatize, sanitize
    - Remove duplicate fragments (if any)
    - Generate 3D coords with ETKDGv3
    - Optimize with MMFF94 (fallback UFF)

    Returns:
    mol (or None), canonical_smiles, actual_forcefield_used
    """

    # -------------Step 1: Initial molecules from SMILES--------------------
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "", ""

    # ------------Step 2: Salt removal & keep largest fragment--------------
    chooser = rdMolStandardize.LargestFragmentChooser()
    mol = chooser.choose(mol)
    if mol is None:
        return None, "", ""

    # ------------Step 3: Normalize functional groups, aromatize, sanitize--
    normalizer = rdMolStandardize.Normalizer()
    mol = normalizer.normalize(mol)
    # Remove duplicate fragments
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        seen = set()
        unique_frags = []
        for f in frags:
            cansmi = Chem.MolToSmiles(f, canonical=True)
            if cansmi not in seen:
                seen.add(cansmi)
                unique_frags.append(f)
        if unique_frags:
            mol = unique_frags[0]
            for f in unique_frags[1:]:
                mol = Chem.CombineMols(mol, f)
    # Ensure aromaticity is perceived after all modifications
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Canonical SMILES for metadate
    canonical_smi = Chem.MolToSmiles(mol, canonical=True)

    # ------------Step 4: 3D ambedding--------------------------------------
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42 # reproducibility
    status = AllChem.EmbedMolecule(mol, params)
    if status !=0:
        logger.warning(f"embedding failes for {smiles}")
        return None, canonical_smi, ""


    # --------Step 5: Force field optimization------------------------------
    used_ff = ""
    try:
        if force_field.upper() == "MMFF94":
            if AllChem.MMFFHasAllMoleculeParams(mol):
                ff = AllChem.MMFFGetMoleculeForceField(mol)
                if ff:
                    ff.Minimize()
                    used_ff = "MMFF94"
    except Exception:
        pass

    if not used_ff:
        # Fallback to UFF
        try:
            ff = AllChem.UFFGetMoleculeForceField(mol)
            ff.Minimize()
            used_ff = "UFF"
        except Exception as e:
            logger.warning(f"UFF Optimization failed for {smiles}: {e}")
            return None, canonical_smi, ""

    mol = Chem.RemoveHs(mol) # hydrogens will be re-added during protonation
    return mol, canonical_smi, used_ff
