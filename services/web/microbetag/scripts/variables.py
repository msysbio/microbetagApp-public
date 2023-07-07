"""
Global variables
"""
import os
import yaml
from datetime import datetime

BASE = "/".join(os.path.abspath(__file__).split("/")[:-2])
IO_PATH = "/tmp"

# temporar output dir
now = datetime.now()
dt_string = now.strftime("%d-%m-%Y--%H.%M.%S")
OUT_DIR = os.path.join(IO_PATH, dt_string)
os.mkdir(OUT_DIR)

# Abundance table
ABUNDANCE_TABLE = os.path.join(IO_PATH, "abundance_table.tsv")
ABUNDANCE_TABLE_DELIM = "\t"
TAX_COL = "taxonomy"
TAX_DELIM = ";"
OTU_COL = "seqId"
COM_CHAR = "#"

# DATABASE
HOST = "db"
USER_NAME = "msysbio"
PASSWORD = "pass"
DB_NAME = "microbetagDB"
DB_PORT = 3306

# Paths to reference database and mapping files
# ----------------------------------------------
REF_DBS = os.path.join(BASE, "ref-dbs")
MAPPINGS = os.path.join(BASE, "mappings")

# GTDB
GTDB_METADATA = os.path.join(
    MAPPINGS, "GTDB_QUALITY_REPRESENTATIVE_GENOMES_v202")

# NCBI Ids paths
SPECIES_NCBI_IDS = os.path.join(MAPPINGS, "species2ncbiId.tsv")
GENERA_NCBI_IDS = os.path.join(MAPPINGS, "genera2ncbiId.tsv")
FAMILIES_NCBI_IDS = os.path.join(MAPPINGS, "families2ncbiId.tsv")

# KEGG paths - Pathway complmementarity module
KMODULES_DEFINITIONS = os.path.join(
    REF_DBS, "kegg_mappings/module_definitions.tsv")
KMODULES_DEFINITIONS_PARSED = os.path.join(
    REF_DBS, "kegg_mappings/module_definition_map.json")
ALL_GENOMES_MODULES = os.path.join(REF_DBS, "all_genomes_modules")

# Tools
# -----
TOOLS = os.path.join(BASE, "tools")
FLASHWEAVE_SCRIPT = os.path.join(TOOLS, "flashweave/flashweave.jl")
FAPROTAX_SCRIPT = os.path.join(TOOLS, "faprotax/collapse_table.py")
FAPROTAX_DB = os.path.join(TOOLS, "faprotax/FAPROTAX.txt")

# FlashWeave
FLASHWEAVE_OUTPUT_DIR = os.path.join(OUT_DIR, "flashweave")
FLASHWEAVE_TMP_INPUT = os.path.join(
    FLASHWEAVE_OUTPUT_DIR,
    "otu_table_flashweave_format.tsv")
FLASHWEAVE_EDGELIST = os.path.join(
    FLASHWEAVE_OUTPUT_DIR,
    "network_output.edgelist")

# FAPROTAX
FAPROTAX_OUTPUT_DIR = os.path.join(OUT_DIR, "faprotax")
FAPROTAX_FUNCT_TABLE = os.path.join(
    FAPROTAX_OUTPUT_DIR, "functional_otu_table.tsv")
FAPROTAX_SUB_TABLES = os.path.join(FAPROTAX_OUTPUT_DIR, "sub_tables")

# PhenDB
PHEN_OUTPUT_DIR = os.path.join(OUT_DIR, "phen_predictions")


# Read configuration file
check = True
try:
    with open("/tmp/config.yml", "r") as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
except BaseException:
    check = False

if check:

    if cfg["taxonomy"] == "GTDB":
        GTDB_ACCESSION_NCBI_TAX_IDS = os.path.join(
            MAPPINGS, "species2ncbiId2accession.tsv")
    elif cfg["taxonomy"] == "dada2":
        GTDB_ACCESSION_NCBI_TAX_IDS = os.path.join(
            MAPPINGS, "dada2ncbi2accession.tsv")
    elif cfg["taxonomy"] == "qiime2":
        GTDB_ACCESSION_NCBI_TAX_IDS = os.path.join(
            MAPPINGS, "qiime2species2ncbi2accession.tsv")
    else:
        GTDB_ACCESSION_NCBI_TAX_IDS = os.path.join(
            MAPPINGS, "gtdbSpecies2ncbiId2accession.tsv")

    if cfg["column_names_are_in"]:
        COM_HEAD = '"' + "last_comment_line" + '"'
    else:
        COM_HEAD = ""

    EDGE_LIST = os.path.join(
        IO_PATH, cfg["edge_list"]) if cfg["edge_list"] is not None else False

    USERS_TAXONOMY = cfg["taxonomy"]
