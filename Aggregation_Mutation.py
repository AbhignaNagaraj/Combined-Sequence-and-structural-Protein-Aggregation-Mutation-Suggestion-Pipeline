import os
import pandas as pd
import matplotlib.pyplot as plt
import freesasa
from Bio import SeqIO

# ==========================================================
# Hydrophobicity definitions
# ==========================================================
KYTE_DOOLITTLE = {
    'A': 1.8,  'C': 2.5,  'D': -3.5, 'E': -3.5,
    'F': 2.8,  'G': -0.4, 'H': -3.2, 'I': 4.5,
    'K': -3.9, 'L': 3.8,  'M': 1.9,  'N': -3.5,
    'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8,
    'T': -0.7, 'V': 4.2,  'W': -0.9, 'Y': -1.3
}

STRONGLY_HYDROPHOBIC = set("AILMFWVY")
SURFACE_THRESHOLD = 0.25   # Relative SASA cutoff


# ==========================================================
# 1. Sequence-based aggregation
# ==========================================================
def sequence_aggregation(sequence, min_run=4):
    flags = [0] * len(sequence)
    start = None

    for i, aa in enumerate(sequence):
        if aa in STRONGLY_HYDROPHOBIC:
            start = i if start is None else start
        else:
            if start is not None and i - start >= min_run:
                for j in range(start, i):
                    flags[j] = 1
            start = None

    if start is not None and len(sequence) - start >= min_run:
        for j in range(start, len(sequence)):
            flags[j] = 1

    return flags


# ==========================================================
# 2. Structure-based aggregation (AlphaFold-safe)
# ==========================================================
def structural_aggregation(pdb_file, seq_len):
    if not os.path.exists(pdb_file):
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    structure = freesasa.Structure(pdb_file)
    result = freesasa.calc(structure)

    # Group atoms by residue in order of appearance
    residues = {}
    residue_order = []

    for i in range(structure.nAtoms()):
        key = (
            structure.chainLabel(i),
            structure.residueNumber(i)
        )
        if key not in residues:
            residues[key] = []
            residue_order.append(key)
        residues[key].append(i)

    struct_flags = [0] * seq_len
    ABSOLUTE_SASA_THRESHOLD = 30.0  # Å²

    for idx, key in enumerate(residue_order):
        if idx >= seq_len:
            break

        atom_indices = residues[key]
        sasa = sum(result.atomArea(i) for i in atom_indices)

        if sasa >= ABSOLUTE_SASA_THRESHOLD:
            struct_flags[idx] = 1

    return struct_flags

# ==========================================================
# 3. Mutation suggestion engine
# ==========================================================
def suggest_mutation(residue):
    mutation_map = {
        "L": "D", "I": "D", "V": "D", "M": "D",
        "F": "N", "W": "N", "Y": "N",
        "A": "S"
    }
    return mutation_map.get(residue)


# ==========================================================
# 4. Combined pipeline
# ==========================================================
def run_combined_pipeline(fasta_file, pdb_file, outdir="results"):
    os.makedirs(outdir, exist_ok=True)

    record = next(SeqIO.parse(fasta_file, "fasta"))
    sequence = str(record.seq)
    protein_id = record.id.replace(" ", "_")

    seq_flags = sequence_aggregation(sequence)
    struct_flags = structural_aggregation(pdb_file, len(sequence))

    rows = []
    for i, aa in enumerate(sequence):
        high_conf = seq_flags[i] == 1 and struct_flags[i] == 1
        mutation = suggest_mutation(aa) if high_conf else ""

        rows.append({
            "Position": i + 1,
            "Residue": aa,
            "Sequence_Aggregation": seq_flags[i],
            "Structural_Aggregation": struct_flags[i],
            "High_Confidence_Aggregation": high_conf,
            "Suggested_Mutation": f"{aa}{i+1}{mutation}" if mutation else ""
        })

    df = pd.DataFrame(rows)

    # ======================================================
    # Save unified CSV
    # ======================================================
    csv_path = os.path.join(
        outdir, f"{protein_id}_aggregation_mutation_analysis.csv"
    )
    df.to_csv(csv_path, index=False)

    # ======================================================
    # Save unified PNG plot
    # ======================================================
    plt.figure(figsize=(10, 3))
    plt.bar(
        df["Position"],
        df["High_Confidence_Aggregation"].astype(int),
        color="red"
    )
    plt.xlabel("Residue Position")
    plt.ylabel("Aggregation Risk")
    plt.title(f"Sequence + Structural Aggregation: {protein_id}")
    plt.tight_layout()

    png_path = os.path.join(
        outdir, f"{protein_id}_aggregation_mutation_analysis.png"
    )
    plt.savefig(png_path, dpi=300)
    plt.close()

    print("Generated files:")
    print(csv_path)
    print(png_path)


# ==========================================================
# 5. Run
# ==========================================================
if __name__ == "__main__":
    run_combined_pipeline(
        fasta_file="input_protein.fasta",
        pdb_file="input_protein.pdb"
    )
