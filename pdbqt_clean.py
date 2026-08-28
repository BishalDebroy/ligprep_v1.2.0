import os
import glob

def remove_model_endmdl_lines(directory=".", inplace=False):
    """
    For every .pdbqt file in the given directory:
      - Delete the first line if it starts with 'MODEL' (case‑sensitive)
      - Delete the last line if it is exactly 'ENDMDL' (ignoring surrounding whitespace)
    
    Parameters
    ----------
    directory : str
        Directory to scan (default: current working directory).
    inplace : bool
        If True, overwrite original files.
        If False (default), save cleaned files as '<original>_cleaned.pdbqt'.
    """
    pattern = os.path.join(directory, "*.pdbqt")
    pdbqt_files = glob.glob(pattern)

    if not pdbqt_files:
        print("No .pdbqt files found in the current directory.")
        return

    for filepath in pdbqt_files:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        if not lines:
            print(f"Skipping empty file: {filepath}")
            continue

        # Determine which lines to keep
        start_idx = 0
        if lines[0].lstrip().startswith("MODEL"):
            start_idx = 1  # drop first line

        end_idx = len(lines)
        if lines[-1].strip() == "ENDMDL":
            end_idx -= 1   # drop last line

        cleaned_lines = lines[start_idx:end_idx]
        modified = (start_idx != 0) or (end_idx != len(lines))

        # Set output path
        if inplace:
            outpath = filepath
        else:
            base, ext = os.path.splitext(filepath)
            outpath = f"{base}_cleaned{ext}"

        # Write the result
        with open(outpath, 'w') as f:
            f.writelines(cleaned_lines)

        if modified:
            actions = []
            if start_idx != 0:
                actions.append("MODEL line removed")
            if end_idx != len(lines):
                actions.append("ENDMDL line removed")
            print(f"{', '.join(actions)}: {filepath} -> {outpath}")
        else:
            print(f"No changes needed: {filepath} (unchanged, saved to {outpath})")


if __name__ == "__main__":
    # Set inplace=True to overwrite originals (use with caution)
    remove_model_endmdl_lines(directory=".", inplace=False)