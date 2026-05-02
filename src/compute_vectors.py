"""
Compute signed inter-residue angles using CA->sidechain centroid vectors.
Filters for all-helix tripeptide windows (HHH pattern).
Groups results by size category of the preceding residue.
"""
import os
import sys
import gzip
import csv
import numpy as np
from Bio.PDB import PDBParser
from tqdm import tqdm

SIZE_CATEGORY = {
    "G": "Tiny",   "A": "Tiny",
    "V": "Small",  "P": "Small",  "S": "Small",  "T": "Small",  "C": "Small",
    "N": "Mid",    "L": "Mid",    "I": "Mid",    "D": "Mid",
    "K": "Large",  "E": "Large",  "M": "Large",  "Q": "Large",  "H": "Large",
    "R": "Bulky",  "F": "Bulky",  "Y": "Bulky",  "W": "Bulky",
}

STRUCT_DIR = os.environ.get("PDB_PATH", "data/pdb_files")


def vector_angle(u, v, ref_axis):
    """Signed angle between u and v around ref_axis (degrees)."""
    u = np.asarray(u, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    ref_axis = np.asarray(ref_axis, dtype=np.float64)
    u = u / np.linalg.norm(u)
    v = v / np.linalg.norm(v)
    ref_axis = ref_axis / np.linalg.norm(ref_axis)
    cp = np.cross(u, v)
    dp = np.dot(u, v)
    return np.degrees(np.arctan2(np.dot(cp, ref_axis), dp))


def open_structure(pdb_id):
    fpath = os.path.join(STRUCT_DIR, f"{pdb_id}.pdb.gz")
    if not os.path.exists(fpath):
        return None
    p = PDBParser(QUIET=True)
    with gzip.open(fpath, "rt") as handle:
        return p.get_structure(pdb_id, handle)


def find_residue(structure, chain_id, resnum):
    for mdl in structure:
        if chain_id not in mdl:
            continue
        for r in mdl[chain_id]:
            if r.id[1] == resnum:
                return r
    return None


def alpha_carbon(residue):
    if residue and "CA" in residue:
        return residue["CA"].get_coord()
    return None


def sidechain_center(residue):
    bb_names = {"N", "CA", "C", "O"}
    atoms = [a.get_coord() for a in residue if a.get_name() not in bb_names]
    if not atoms:
        return None
    return np.mean(atoms, axis=0)


def run(window_dir, amino_acid, output_file):
    tsv_files = [f for f in os.listdir(window_dir) if f.endswith(".tsv")]
    data = []

    for fname in tqdm(tsv_files, desc="computing angles"):
        fp = os.path.join(window_dir, fname)
        if os.path.getsize(fp) == 0:
            continue
        try:
            with open(fp) as fh:
                rows = list(csv.reader(fh, delimiter="\t"))
        except Exception:
            continue
        if not rows:
            continue

        pdb_id = rows[0][12]
        struct = open_structure(pdb_id)
        if struct is None:
            continue

        pos = 0
        while pos + 2 < len(rows):
            r_prev, r_mid, r_next = rows[pos], rows[pos+1], rows[pos+2]
            pos += 3

            if r_mid[0] != amino_acid:
                continue
            if r_mid[11] != "HHH":
                continue

            chain = r_mid[1]
            try:
                n_prev = int(r_prev[2])
                n_mid = int(r_mid[2])
                n_next = int(r_next[2])
            except (ValueError, IndexError):
                continue

            res_p = find_residue(struct, chain, n_prev)
            res_m = find_residue(struct, chain, n_mid)
            res_n = find_residue(struct, chain, n_next)
            if not all([res_p, res_m, res_n]):
                continue

            ca_p = alpha_carbon(res_p)
            ca_m = alpha_carbon(res_m)
            ca_n = alpha_carbon(res_n)
            sc_p = sidechain_center(res_p)
            sc_m = sidechain_center(res_m)

            if any(x is None for x in [ca_p, ca_m, ca_n, sc_p, sc_m]):
                continue

            vec_prev = sc_p - ca_p
            vec_mid = sc_m - ca_m
            backbone_dir = ca_m - ca_p

            if np.linalg.norm(backbone_dir) < 1e-9:
                continue

            try:
                ang = vector_angle(vec_prev, vec_mid, backbone_dir)
            except Exception:
                continue

            neighbor_aa = r_prev[9]
            cat = SIZE_CATEGORY.get(neighbor_aa)
            if cat is None:
                continue

            data.append([pdb_id, neighbor_aa, cat, f"{ang:.4f}"])

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["pdb_id", "neighbor", "size_cat", "angle_deg"])
        for rec in data:
            w.writerow(rec)

    # write contributing PDB IDs
    id_file = os.path.join(os.path.dirname(output_file), "structures_used.txt")
    ids = sorted(set(r[0] for r in data))
    with open(id_file, "w") as f:
        f.writelines(i + "\n" for i in ids)

    print(f"wrote {len(data)} measurements -> {output_file}")
    print(f"{len(ids)} structures contributed data")


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], sys.argv[3])
