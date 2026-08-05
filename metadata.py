"""
metadata.py - Build SD tags for first and second SDF entries.
"""

from typing import Dict, Optional, Any
from importlib.metadata import version, PackageNotFoundError
from rdkit import Chem

try:
    _DIMORPHITE_VERSION = version("dimorphite-dl")
except PackageNotFoundError:
    _DIMORPHITE_VERSION = "unknown"

def first_sdf_tags(
        compound_name: str,
        original_smi: str,
        canonical_smi: str,
        force_field: str,
        ph: float,
        uuid: str,
        generation_date: str,
) -> Dict[str, str]:
    """Return dictionary of SD tad name → value for first SDF."""
    return{
        "Compound_Name": compound_name,
        "Original_SMILES": original_smi,
        "Canonical_SMILES": canonical_smi,
        "Generation_Date": generation_date,
        "Preparation_Software": "RDKit;Dimorphite-DL",
        "RDKit_version": Chem.rdBase.rdkitVersion,
        "Dimorphite-DL_version": _DIMORPHITE_VERSION,
        "Force_field": force_field,
        "pH": str(ph),
        "Processing_Step": "standardization+protonation",
        "UUID": uuid
    }

def second_sdf_tags(
        parent_tags: Dict[str, str],
        tautomer_name: str,
        tautomer_number: int,
        is_canonical: bool,
        energy: Optional[float],
        energy_unit: str,
        opt_force_field: str,
) -> Dict[str, str]:
    """Return tafs for a tautomer entry, retaining all first SDF fields."""
    new_tag = dict(parent_tags) # Copy all original fields
    new_tag.update({
        "Tautomer_Name": tautomer_name,
        "Tautomer_Number": str(tautomer_number),
        "Parent_compound": parent_tags["Compound_Name"],
        "Is_Canonical_Tautomer": "1" if is_canonical else "0",
        "Energy": f"{energy:.4f}" if energy is not None else "N/A",
        "Energy_Unit": energy_unit,
        "Optimization_ForceField": opt_force_field,
    })
    return new_tag