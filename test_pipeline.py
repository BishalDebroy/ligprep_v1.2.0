"""
Usage: python -m pytest test_pipeline.py
"""

import unittest
from rdkit import Chem
from rdkit.Chem import AllChem
from standardize import standardize_and_embed
from protonate import protonate
from tautomer import enumerate_tautomers
from optimizer import optimize_and_energy

class TestLigandPrep(unittest.TestCase):
    def test_standardize_success(self):
        mol, cansmi, ff = standardize_and_embed("CCO")
        self.assertIsNotNone(mol)
        self.assertEqual(cansmi, "CCO")
        self.assertIn(ff, ["MMFF94", "UFF"])

    def test_standardize_invalid(self):
        mol, _, _ = standardize_and_embed("X")
        self.assertIsNone(mol)

    def test_protonation(self):
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)
        prot = protonate(mol, pH=7.4)
        self.assertIsNotNone(prot)
        self.assertGreaterEqual(prot.GetNumAtoms(), mol.GetNumAtoms())

    def test_tautomer_enumeration(self):
        mol = Chem.MolFromSmiles("CC(=O)O")
        taut = enumerate_tautomers(mol)
        self.assertTrue(len(taut) >= 1)
        self.assertTrue(taut[0][1])

    def test_optimization_energy(self):
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        opt_mol, energy, ff = optimize_and_energy(mol)
        self.assertIsNotNone(opt_mol)
        self.assertIn(ff, ["MMFF94", "UFF"])

if __name__ == "__main__":
    unittest.main()