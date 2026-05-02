"""Compute vector angles from extracted windows."""
import glob

configfile: "settings.yml"

RESIDUE = config["residue"]
WINDOWS = [f.split("/")[-1].replace(".tsv", "") for f in glob.glob("tmp/windows/*.tsv")]

rule all:
    input:
        "output/angles.tsv"

rule compute:
    input:
        expand("tmp/windows/{w}.tsv", w=WINDOWS)
    output:
        "output/angles.tsv"
    params:
        res = RESIDUE
    shell:
        "python3 src/compute_vectors.py tmp/windows {params.res} {output}"
