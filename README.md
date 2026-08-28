# Large-scale Ligand Preparation Pipeline

A production quality Python pipeline that converts a SMILES file into two multi-ligand SDF files suitable for virtual screening:

1. **Standardised & protonated 3D structures** (`output.sdf`)
2. **Tautomer-enriched & energy-optimized structures** (`output_taut.sdf`)

Designed for millions of molecules, it streams data, uses multiprocessing, and fully fault tolerant.

## Features
- Reads SMILES files (tab, comma, or space separated)
- Automatically discards invalid molecules and logs errors
- Salt removal, largest fragment selection, functional group normalization
- 3D coordinate generation (ETKDGv3) followed by MMFF94 (fallback UFF) optimization
- Protonation state assignment at any pH (using Dimorphite-DL)
- Tautomer enumeration and ranking (canonical tautomer highlighted)
- Energy minimization of each tautomer
- Rich metadata (UUID, force field, energy, etc.) in the SDF tags
- Streaming processing - memory footprint constant, works on huge libraries
- Multiprocessing for CUP-heavy steps
- Configurable via command line

## Dependencies
- Python ≥ 3.10
- RDKit
- Dimorphite-DL
- psutil (for memory logging)

Install via conda (recommended) or pip:

```
conda create -n ligprep_env python=3.10 rdkit -c conda forge
conda activate ligprep_env
pip install dimorphite-dl psutil
```

## Quick Start

1. Activate the env
```
conda activate ligprep_env
```

2. Process SMILES and write SDF files
```
python3 run_ligprep.py example_input.smi output_test -o results
```

3. Overwrite the existing content
```
python3 run_ligprep.py example_input.smi output_test -o results --overwrite
```

4. Recommanded usage
```
python run_ligprep.py example_input.smi output_test --ph 7.4 --ff UFF --num-cpus 18 -o results -l resultslog.txt -v
```
optimize the value of `--ph`, `--ff`, and `--num-cpus` as per requirement

6. For advanced usage, use this command to open different flags
```
python3 run_ligprep.py --help
```

6. For perform a test run use this command
```
python3 -m unittest test_pipeline.py -v
```

7. Deactivate the conda env
```
conda deactivate
```


## Performance expected
- Memory usage under 500 MB regardless of input size (up to 2 million molecules)
- Linear speedup with CPU cores (95% parallel efficiency on typical hardware)
- Typical throughput: ~5000 molecules/minute of 16 cores (including tautomer steps)

# Further Task
SDF library should be now prepared for docking. This required conversion of SDF into PDBQT, easily performed by openbabel

1. Install openbabel in the env
```
conda activate ingprep_env
conda install -c conda-forge openbabel
obabel -V
```
2. Converting SDF library when ligand name not provided in the input.smi file, this command will name the ligand as `ligand_1.pdbqt`, `ligand_2.pdbqt` and so on
```
obabel output_test_taut.sdf -O ligand_.pdbqt -m --minimize --ff UFF
```
3. Converting SDF library with ligand name in the input.smi file, this command will replace `ligand_`, with the original name in metadata.
```
obabel output_test_taut.sdf -O ligand_.pdbqt --split --minimize --ff UFF
```
# The ligand(s) are ready for docking
