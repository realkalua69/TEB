"""
Extract tripeptide windows from STRIDE secondary structure output.
"""
import sys
import re
import csv
import os

AA_MAP = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def parse_stride_file(path):
    entries = []
    with open(path) as fh:
        for ln in fh:
            if not ln.startswith("ASG"):
                continue
            parts = ln.split()
            if len(parts) < 10:
                continue
            rnum_match = re.match(r"\d+", parts[3])
            if rnum_match is None:
                continue
            try:
                entries.append({
                    "res": parts[1],
                    "chain": parts[2],
                    "resnum": int(rnum_match.group()),
                    "ordinal": int(parts[4]),
                    "ss_code": parts[5],
                    "ss_name": parts[6],
                    "phi": float(parts[7]),
                    "psi": float(parts[8]),
                    "sasa": float(parts[9]),
                })
            except (ValueError, IndexError):
                pass
    return entries


def extract_windows(entries, target_res, pdb_tag):
    output = []
    for idx in range(1, len(entries) - 1):
        center = entries[idx]
        if center["res"] != target_res:
            continue

        prev = entries[idx - 1]
        after = entries[idx + 1]

        triplet_seq = "".join(AA_MAP.get(e["res"], "?") for e in [prev, center, after])
        ss_pattern = prev["ss_code"] + center["ss_code"] + after["ss_code"]

        location = f"{center['chain']}:{prev['resnum']},{center['chain']}:{center['resnum']},{center['chain']}:{after['resnum']}"

        for e in [prev, center, after]:
            output.append([
                e["res"], e["chain"], e["resnum"], e["ordinal"],
                e["ss_code"], e["ss_name"], e["phi"], e["psi"], e["sasa"],
                AA_MAP.get(e["res"], "?"),
                triplet_seq, ss_pattern, pdb_tag, location,
            ])
    return output


if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    target = sys.argv[3]

    tag = os.path.basename(infile).replace(".stride", "").replace(".ss.out", "").split(".")[0]
    records = parse_stride_file(infile)
    windows = extract_windows(records, target, tag)

    with open(outfile, "w", newline="") as f:
        wr = csv.writer(f, delimiter="\t")
        for row in windows:
            wr.writerow(row)
