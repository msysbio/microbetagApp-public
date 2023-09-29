"""
Global variables
"""
import os
import sys

BASE = "/".join(os.path.abspath(__file__).split("/")[:-2])
IO_PATH = "/tmp"

# Abundance table
ABUNDANCE_TABLE_DELIM = "\t"
TAX_COL = "taxonomy"
TAX_DELIM = ";"
SEQ_COL = "seqId"
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
