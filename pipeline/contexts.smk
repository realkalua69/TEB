"""Extract tripeptide windows from STRIDE output."""
import glob

configfile: "settings.yml"

RESIDUE = config["residue"]
STRIDE_FILES = [f.split("/")[-1].replace(".stride", "") for f in glob.glob("tmp/stride/*.stride")]

rule all:
    input:
        expand("tmp/windows/{res}_{pdb}.tsv", pdb=STRIDE_FILES, res=[RESIDUE])

rule extract:
    input:
        "tmp/stride/{pdb}.stride"
    output:
        "tmp/windows/{res}_{pdb}.tsv"
    params:
        res = "{res}"
    shell:
        "python3 src/stride_parser.py {input} {output} {params.res}"
