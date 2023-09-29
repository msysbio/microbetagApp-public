#!/usr/bin/env python3

from microbetag.scripts.utils import *
from microbetag.scripts.db_functions import *
from microbetag.scripts.logger import *
from microbetag.scripts.variables import *
import os
import sys
import json

version = "v0.1.0"
license = ("GPLv3",)
packages = ["microbetag"]
description = """
            microbetag is a microbial network annotator that exploits several software and databases to
            annotate microbial, co-occurence networks.
            Functional traits like whether a species is methanotroph, fermentative etc are assigned to
            each node of the network.

            For the associations that include 2 taxa of the species/strain level, microbetag also performs a
            pathway complementarity step; species are assigned to their corresponding GTDB representative
            genomes and using those genomes, complementarity for KEGG modules are investigated.

            Finally, microbetag supports a series of metabolic models covering the GTDB representative genomes.
            Metabolic analysis methods such as FBA and flux sampling are performed to highlight potential metabolic
            interactions among the taxa.
"""
author = "Haris Zafeiropoulos"
author_email = "haris.zafeiropoulos@kuleuven.be"
name = "microbetag"


def main(out_dir, cfg, abundance_table_file=None, edge_list_file=None):
    """
    Setting logging; the logger variable comes from the logger.py script
    """
    if cfg["taxonomy"] == "GTDB":
        gtdb_accession_ncbi_tax_id = os.path.join(
            MAPPINGS, "species2ncbiId2accession.tsv")
    elif cfg["taxonomy"] == "dada2":
        gtdb_accession_ncbi_tax_id = os.path.join(
            MAPPINGS, "dada2ncbi2accession.tsv")
    elif cfg["taxonomy"] == "qiime2":
        gtdb_accession_ncbi_tax_id = os.path.join(
            MAPPINGS, "qiime2species2ncbi2accession.tsv")
    else:
        gtdb_accession_ncbi_tax_id = os.path.join(
            MAPPINGS, "gtdbSpecies2ncbiId2accession.tsv")

    # Using FileHandler writing log to file
    logfile = os.path.join(out_dir, "log.txt")
    open(logfile, "w")
    fh = logging.FileHandler(logfile)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)

    # Using StreamHandler writing to console
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)

    # Using error
    eh = logging.StreamHandler()
    eh.setLevel(logging.ERROR)
    eh.setFormatter(formatter)

    # Add the two Handlers
    logger.addHandler(ch)
    logger.addHandler(fh)
    logger.addHandler(eh)

    """
    Welcome message and arguments values
    """
    logging.info("Hello microbe-fun! microbetag is about to start!")
    logging.info("Your command was: {}".format(" ".join(sys.argv)))

    # STEP 0: Load vars

    flashweave_output_dir = os.path.join(out_dir, "flashweave")
    flashweave_tmp_input = os.path.join(
        flashweave_output_dir,
        "abundance_table_flashweave_format.tsv")
    flashweave_edgelist = os.path.join(
        flashweave_output_dir,
        "network_output.edgelist")

    faprotax_output_dir = os.path.join(out_dir, "faprotax")
    faprotax_funct_table = os.path.join(
        faprotax_output_dir, "functional_otu_table.tsv")
    faprotax_sub_tables = os.path.join(faprotax_output_dir, "sub_tables")

    phen_output_dir = os.path.join(out_dir, "phen_predictions")
    pathway_complementarity = os.path.join(out_dir, "path_compl")

    # STEP 1: abundance table preprocess
    if abundance_table_file:

        abundance_table = is_tab_separated(abundance_table_file, TAX_COL)
        logging.info(
            "STEP 1: Assign NCBI Tax Id and GTDB reference genomes".center(80, "*"))
        if not edge_list_file:

            logging.info(
                "The user has not provided an edge list. microbetag will build one using FlashWeaeve.")
            if not os.path.exists(flashweave_output_dir):
                os.mkdir(flashweave_output_dir)

            ext = ensure_flashweave_format(abundance_table, TAX_COL, SEQ_COL, flashweave_output_dir)

        else:
            ext = abundance_table.copy()
            ext["microbetag_id"] = abundance_table[SEQ_COL]

        # Map taxonomies to ontology ids
        logging.info(
            "Get the NCBI Taxonomy id for those OTUs that have been assigned either at the species, the genus or the family level.")
        try:
            seqID_taxid_level_repr_genome, repr_genomes_present = map_seq_to_ncbi_tax_level_and_id(ext, TAX_COL, SEQ_COL, gtdb_accession_ncbi_tax_id)

        except BaseException:
            logging.error("""microbetag was not able to map your table's taxonomies to NCBI and GTDB ids.
                           Check on how you describe your table.
                           Also check whether you have set the edge_list_file parameter in the config file; if there is not an edge list, leave it blank.""")
            sys.exit(0)

        seqID_taxid_level_repr_genome.to_csv(
            os.path.join(
                out_dir,
                "taxa_with_ndbiId_and_repr_gtdbIds.csv"),
            "\t")

    # STEP 2: Get co-occurrence network
    logging.info("STEP 2: Get co-occurrence network".center(80, "*"))
    if not edge_list_file:

        flashweave_params = [
            "julia", FLASHWEAVE_SCRIPT, flashweave_output_dir, flashweave_tmp_input
        ]
        flashweave_command = " ".join(flashweave_params)
        logging.info("Run FlashWeave")
        os.system(flashweave_command)

        # Taxa pairs as NCBI Tax ids
        logging.info(
            "Map your edge list to NCBI Tax ids and keep only associations that both correspond to a such."
        )
        edge_list = edge_list_of_ncbi_ids(flashweave_edgelist, seqID_taxid_level_repr_genome)

        # Example:
        # [{'node_a': 'microbetag_17', 'ncbi_tax_id_node_a': 77133, 'gtdb_gen_repr_node_a': 'GCA_903925685.1', 'ncbi_tax_level_node_a': 'mspecies',
        #   'node_b': 'microbetag_21', 'ncbi_tax_id_node_b': 136703, 'gtdb_gen_repr_node_b': nan, 'ncbi_tax_level_node_b': 'species'},
        # {'node_a': 'microbetag_17', 'ncbi_tax_id_node_a': 77133, 'gtdb_gen_repr_node_a': 'GCA_903925685.1', 'ncbi_tax_level_node_a': 'mspecies',
        #  'node_b': 'microbetag_74', 'ncbi_tax_id_node_b': 77133, 'gtdb_gen_repr_node_b': 'GCA_903925685.1', 'ncbi_tax_level_node_b': 'mspecies'}
        edgelist_as_a_list_of_dicts = edge_list.applymap(lambda x: str(x) if pd.notna(x) else 'null').to_dict(orient="records")

        edge_list.to_csv("edgelist.csv", sep="\t")

        microb_id_taxonomy = ext[["microbetag_id", "taxonomy"]]

        base_network = build_base_graph(edgelist_as_a_list_of_dicts, microb_id_taxonomy)

        with open(os.path.join(out_dir, "basenet.json"), "w") as f:
            json.dump(base_network, f)

        logging.info("Base network has been built and saved.")

    # STEP 3: PhenDB
    logging.info("STEP 3: PhenDB ".center(80, "*"))
    if cfg["phenDB"]:

        # Get phen traits for each GTDB genome present on your table
        feats = get_column_names()
        feats.insert(1, "NCBI_ID")
        feats.insert(1, "Species")
        traits = []
        for gtdb_id in set(repr_genomes_present):
            check = True
            try:
                r = get_phenDB_traits(gtdb_id)
            except BaseException:
                check = False
            if not check:
                if gtdb_id[:3] == "GCA":
                    try:
                        r = get_phenDB_traits(gtdb_id.replace("GCA", "GCF"))
                        check = True
                    except BaseException:
                        pass
                else:
                    try:
                        r = get_phenDB_traits(gtdb_id.replace("GCF", "GCA"))
                    except BaseException:
                        pass
            if check:
                tr = seqID_taxid_level_repr_genome.loc[seqID_taxid_level_repr_genome["gtdb_gen_repr"] == gtdb_id, [
                    "Species", "ncbi_tax_id"]]
                sp = tr.iloc[0, 0]
                ncbiId = str(int(tr.iloc[0, 1]))
                r.insert(1, ncbiId)
                r.insert(1, sp)
                traits.append(r)
            else:
                print(
                    "Genome :",
                    gtdb_id,
                    "is not present in the phenDB version of microbetag.")

        # Save phen traits as a .tsv file
        if not os.path.exists(phen_output_dir):
            os.mkdir(phen_output_dir)
        outfile = os.path.join(phen_output_dir, "phen_traits.tsv")
        export_phen_traits_to_file(
            column_names=feats,
            rows=traits,
            filename=outfile)

    # STEP 4: FAPROTAX
    logging.info(
        "STEP 4: FAPROTAX database oriented analaysis".center(80, "*"))
    if cfg["faprotax"] and abundance_table_file:

        if not os.path.exists(faprotax_output_dir):
            os.mkdir(faprotax_output_dir)

        faprotax_params = [
            "python3", FAPROTAX_SCRIPT,
            "-i", abundance_table_file,
            "-o", faprotax_funct_table,
            "-g", FAPROTAX_DB,
            "-c", '"' + COM_CHAR + '"',
            "-d", '"' + TAX_COL + '"',
            "-v",
            "--force",
            "-s", faprotax_sub_tables,
        ]

        faprotax_command = " ".join(faprotax_params)

        # Run FAPROTAX
        # In the sub tables files, in column 1 we have the taxonomy and in column 2 the OTU ID.
        if os.system(faprotax_command) != 0:
            logging.error("""Something went wrong when running the FAPROTAX analysis!
                            Check how you describe your input table.
                            Also, please make sure you have set the column_names_are_in parameter properly.""")
            sys.exit(0)

    # STEP 5: PATHWAY COMPLEMENTARITY
    logging.info("STEP 5: Pathway complementarity".center(80, "*"))
    if cfg["pathway_complement"]:

        if edge_list_file:
            edge_list = edge_list_file.copy()

        ncbi_id_pairs_with_complements = get_complements_of_list_of_pair_of_ncbiIds(edgelist_as_a_list_of_dicts)
        os.mkdir(pathway_complementarity)

        with open(os.path.join(pathway_complementarity, "complements.json"), "w") as f:
            ncbi_id_pairs_with_complements = {tuple_to_str(key): value for key, value in ncbi_id_pairs_with_complements.items()}
            json.dump(ncbi_id_pairs_with_complements, f)

    logging.info("Pathway complementarity has been completed successfully.".center(80, "*"))

    # STEP : BUILD OUTPUT ANNOTATED NETWORK
    annotated_network = annotate_network(base_network, cfg, seqID_taxid_level_repr_genome, out_dir)

    with open(os.path.join(out_dir, "annotated_network.json"), "w") as f:
        json.dump(annotated_network, f)

    return annotated_network

#     STEP 7: Net Cooperate
#     logging.info("STEP 7: NetCooperate between species/strains paired nodes".center(80, "*"))
#     if NETCOOPERATE:

    """
    # STEP 5: BugBase

    logging.info("STEP 5: BugBase database oriented analaysis".center(80, "*"))
    if BUGBASE and abundance_table_file:

        # Make a copy of the otu table without the taxonomy column
        f = open(abundance_table_file, "r")
        g = open(BUGBASE_TMP, "w")
        for line in f:
            g.write("\t".join(line.split("\t")[:-1]) + "\n")

        bugbase_params = [
            "Rscript", BUGBASE_SCRIPT,
            "-i", BUGBASE_TMP,
            "-o", BUGBASE_OUTPUT,
            "-a",
        ]

        if METADATA_FILE:
            bugbase_params = bugbase_params + ["-m", METADATA_FILE]

        bugbase_command = " ".join(bugbase_params)
        print(bugbase_command)
        print(">>>", BUGBASE_PATH)

        # Run BugBase
        logging.info(["Command to run: ", bugbase_command])
        if os.system(bugbase_command) != 0:
            logging.error(
                "\nSomething went wrong when running the BugBase analysis!")
            sys.exit(0)


        # [TODO]: - Parse the bugbase/otu_contributions/contributing_otus.txt to assign features in the OTUs

        os.remove(BUGBASE_TMP)

        # for tr_file in os.listdir(BUGBASE_OUTPUT)

    sys.exit(0)
    """
