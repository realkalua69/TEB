"""Run STRIDE on a batch of PDB structures."""

configfile: "settings.yml"

STRUCT_PATH = config.get("pdb_path", "data/pdb_files")
BATCH = config.get("batch_file", "batch.txt")

with open(BATCH) as f:
    TARGETS = [l.strip() for l in f if l.strip()]

rule all:
    input:
        expand("tmp/stride/{name}.stride", name=TARGETS)

rule unzip:
    input:
        STRUCT_PATH + "/{name}.pdb.gz"
    output:
        temp("tmp/pdb/{name}.pdb")
    shell:
        "gzip -dc {input} > {output}"

rule stride:
    input:
        "tmp/pdb/{name}.pdb"
    output:
        "tmp/stride/{name}.stride"
    shell:
        "stride {input} > {output}"
