"""
main.py - Ligand Preparation Pipeline.
Basic usage: python run_ligprep.py example_input.smi output_name --ph 7.4 -o results
"""
import argparse
import io
import logging
import multiprocessing
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional, List, Tuple, Dict, Any

from rdkit import Chem

from reader import read_smiles, read_sdf_blocks
from standardize import standardize_and_embed
from protonate import protonate
from tautomer import enumerate_tautomers
from optimizer import optimize_and_energy
from metadata import first_sdf_tags, second_sdf_tags
from writer import mol_with_tags_to_sdf_block, write_sdf_blocks
from utils import setup_logging, unique_id, current_date_tag, memory_usage_mb, log_memory, start_queue_listener, worker_logging_init

__version__ = "1.2.0"

logger = logging.getLogger("ligprep")

# ---------------------------------------------------------------------
# Worker functions for parallel stages
# ---------------------------------------------------------------------
def process_smiles_worker(args: Tuple[str, str, int, float, Optional[float], Optional[float], int, str]) -> Optional[Tuple[str, None]]:
    """
    Process a single SMILES entry.
    Returns (mol_block_sdf, None) or None on failure.
    """
    smi, name, line_no, pH, min_pH, max_pH, max_states, force_field = args
    try:
        mol, cansmi, used_ff = standardize_and_embed(smi, force_field)
        if mol is None:
            raise ValueError ("Standardization or embedding failed")

        # Protonation
        prot_mol = protonate(mol, pH, min_pH, max_pH, max_states)
        if prot_mol is None:
            raise ValueError("Protonation failed")

        # Generate metadata
        uid = unique_id()
        date_tag = current_date_tag()
        tags = first_sdf_tags(
            compound_name=name if name else f"mol_{line_no}",
            original_smi=smi,
            canonical_smi=cansmi,
            force_field=used_ff,
            ph=pH,
            uuid=uid,
            generation_date=date_tag
        )

        block = mol_with_tags_to_sdf_block(prot_mol, tags)
        return block, None
    except Exception as e:
        logger.error(f"Failed molecule at line {line_no}: {smi} - {e}")
        return None

def process_tautomer_worker(args: Tuple[str, int, str]) -> List[str]:
    """
    Read one SDF block, enumerate tautomers, optimize, return list of SDF blocks.
    """
    block_str, entry_idx, force_field = args
    try:
        supplier = Chem.ForwardSDMolSupplier(io.BytesIO(block_str.encode("utf-8")), removeHs=False)
        mol = next(supplier, None)
    except Exception as e:
        logger.error(f"Could not parse SDF entry {entry_idx}: {e}")
        return []
    if mol is None:
        logger.error(f"Could not parse SDF entry {entry_idx}")
        return []

    # Recover parent tags from the SDF properties
    parent_tags = mol.GetPropsAsDict()
    parent_compound = parent_tags.get("Compound_Name", f"mol_{entry_idx}")

    # Generate tautomers
    taut_list = enumerate_tautomers(mol)
    blocks_out = []
    for taut_num, (taut_mol, is_canon) in enumerate(taut_list, start=1):
        # Optimize tautomer
        opt_mol, energy, opt_ff = optimize_and_energy(taut_mol, force_field)
        if opt_mol is None:
            logger.warning(f"Optimization failed for tautomer {taut_num} of entery {entry_idx}")
            continue
        taut_name = f"{parent_compound}_taut{taut_num:04d}"
        tags = second_sdf_tags(
            parent_tags=parent_tags,
            tautomer_name=taut_name,
            tautomer_number=taut_num,
            is_canonical=is_canon,
            energy=energy,
            energy_unit="kcal/mol",
            opt_force_field=opt_ff,
        )
        blk = mol_with_tags_to_sdf_block(opt_mol, tags)
        blocks_out.append(blk)
    return blocks_out

# ---------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Ligand preparation pipeline")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("input_smi", help="input SMILES file")
    parser.add_argument("output_basename", help="Basename for output SDF files")
    parser.add_argument("--ph", type=float, default=7.4, help="pH for protonation (default 7.4)")
    parser.add_argument("--min_ph", type=float, help="Minimum pH for protonation range")
    parser.add_argument("--max_ph", type=float, help="Maximum pH for protonation range")
    parser.add_argument("-ff", default="MMFF94", choices=["MMFF94", "UFF"], help="Prefered force field (default MMFF94)")
    parser.add_argument("--num-cpus", type=int, default=multiprocessing.cpu_count(), help="Number of CPU cores to use")
    parser.add_argument("-o", default=".", help="Directory for output files")
    parser.add_argument("--overwrite", action="store_true", help="Overwriting existing output files")
    parser.add_argument("-l", dest="log_file", default=None, help="Log file path (default: <output_dir>/<output_basename>_<DDMONYY>.log)")
    parser.add_argument("-v", action="store_true", help="Verbose logging")
    parser.add_argument("--max_protonation_states", type=int, default=1, help="Number of protonation states to keep (default 1)")
    parser.add_argument("--batch_size", type=int, default=1000, help="Batch size for parallel processing")
    args = parser.parse_args()

    # Output directory must exist before we can default the log path into it
    out_dir = Path(args.o)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    if args.log_file:
        log_path = args.log_file
    else:
        log_path = str(out_dir / f"{args.output_basename}_{current_date_tag()}.log")
    logger = setup_logging(log_path, args.v)
    logger.info(f"Ligand Preparation Pipeline v{__version__}")
    logger.info("Ligand Preparation Started")
    start_time = datetime.now()
    log_memory(logger, "Start")

    # Multiprocessing-safe logging: workers route log records through this
    # queue instead of writing directly to the file/console handles they'd
    # otherwise inherit via fork (which causes concurrent-write/fd errors).
    log_manager = multiprocessing.Manager()
    log_queue = log_manager.Queue(-1)
    log_listener = start_queue_listener(logger, log_queue)

    first_sdf_path = out_dir / f"{args.output_basename}.sdf"
    second_sdf_path = out_dir / f"{args.output_basename}_taut.sdf"
    invalid_path = out_dir / "invalid_smiles.csv"

    # Overwrite check
    if not args.overwrite:
        for p in [first_sdf_path, second_sdf_path]:
            if p.exists():
                logger.error(f"Output file {p} already exists. Use --overwrite to replace.")
                sys.exit(1)

    # ------- Stage 1: SMILE → standardize + protonated SDF ---------
    logger.info("Stage 1: SMILE to 3D generation")

    # prepare pool
    pool = multiprocessing.Pool(
        processes=args.num_cpus,
        initializer=worker_logging_init,
        initargs=(log_queue, args.v),
    )
    # Collect invalid lines
    invalid_entries = []

    # We'll firsy use imap_unordered to keep memory low
    # first, create a generator that yield task tuples
    def task_generator():
        for smi, name, lineno in read_smiles(args.input_smi):
            yield (smi, name, lineno, args.ph, args.min_ph, args.max_ph, args.max_protonation_states, args.ff)

    # Process in batches to control memory
    batch = []
    total_processed = 0
    total_skipped = 0

    def write_sdf_block(blk):
        nonlocal total_processed
        if blk:
            write_sdf_blocks([blk], str(first_sdf_path), append=(total_processed > 0))
            total_processed += 1

    # We'll collect results from pool and write them sequentially.
    for result in pool.imap_unordered(process_smiles_worker, task_generator(), chunksize=args.batch_size):
        if result is None:
            total_skipped += 1
            continue
        block, _= result
        write_sdf_block(block)

    pool.close()
    pool.join()

    # Write invalid SMILES file
    with open(invalid_path, "w") as f:
        f.write("line,SMILE,NAME,Error\n")
        # In this simple implementation we don't collect per-line errors in the workers;
        # we log them and count skipped. For completeness you could extent.
        f.write(f"Total skipped: {total_skipped}\n")
    logger.info(f"Stage 1 finished. processed {total_processed} molecules, skipped {total_skipped} molecules.")
    log_memory(logger, "after stage 1")


    # --------------- Stage 2: Tautomer generation ------------------
    logger.info("Stage 2: Tautomer enumeration and optimizaton")

    if not first_sdf_path.exists() or first_sdf_path.stat().st_size == 0:
        logger.warning("No multi-ligand 3D SDF found, skipping tautomer generation.")
    else:
        pool2 = multiprocessing.Pool(
            processes=args.num_cpus,
            initializer=worker_logging_init,
            initargs=(log_queue, args.v),
        )
        # Task generator from first SDF
        def sdf_task_gen():
            for blk, idx in read_sdf_blocks(str(first_sdf_path)):
                yield (blk, idx, args.ff)

        total_tautomers = 0
        with open(second_sdf_path, "w") as fh_out:
            # Use imap_unordered and write blocks as they come.
            # we need to write blocks sequencially; we can collect and write in order?
            # Since order doesn't matter, we write immediately.
            for taut_blocks in pool2.imap_unordered(process_tautomer_worker, sdf_task_gen(), chunksize=args.batch_size):
                if not taut_blocks:
                    continue
                for blk in taut_blocks:
                    fh_out.write(blk)
                    total_tautomers += 1

        pool2.close()
        pool2.join()
        logger.info(f"Stage 2 finished. Total tautomers written: {total_tautomers}")

    end_time = datetime.now()
    elapsed = end_time - start_time
    logger.info(f"Ligand Preparation completed. elapsed time: {elapsed}")
    log_memory(logger, "final")
    logger.info(f"Finished at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    log_listener.stop()


if __name__ == "__main__":
    main()