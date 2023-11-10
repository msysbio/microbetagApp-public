#!/usr/bin/env python3

from microbetag.scripts.utils import *
from microbetag.scripts.db_functions import *
from microbetag.scripts.logger import *
from microbetag.scripts.variables import *
import os
import sys
import json
import networkx as nx

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


def main(out_dir, cfg, abundance_table_file=None, edge_list_file=None, metadata_file=None):
    """
    Setting logging; the logger variable comes from the logger.py script
    """
    # Using FileHandler writing log to file
    logfile = os.path.join(out_dir, "log.txt")
    open(logfile, "w")
    fh = logging.FileHandler(logfile); fh.setLevel(logging.DEBUG); fh.setFormatter(formatter)

    # Using StreamHandler writing to console
    ch = logging.StreamHandler(); ch.setLevel(logging.INFO); ch.setFormatter(formatter)
    # Using error
    eh = logging.StreamHandler(); eh.setLevel(logging.ERROR); eh.setFormatter(formatter)
    # Add the two Handlers
    logger.addHandler(ch); logger.addHandler(fh); logger.addHandler(eh)

    # Load vars
    flashweave_output_dir = os.path.join(out_dir, "flashweave")
    flashweave_tmp_input = os.path.join(flashweave_output_dir, "abundance_table_flashweave_format.tsv")
    flashweave_edgelist = os.path.join(flashweave_output_dir, "network_output.edgelist")

    faprotax_output_dir = os.path.join(out_dir, "faprotax")
    faprotax_funct_table = os.path.join(faprotax_output_dir, "functional_otu_table.tsv")
    faprotax_sub_tables = os.path.join(faprotax_output_dir, "sub_tables")

    phen_output_dir = os.path.join(out_dir, "phen_predictions")
    pathway_complementarity_dir = os.path.join(out_dir, "path_compl")
    seed_scores_dir = os.path.join(out_dir, "seed_scores")

    # Fix arguments
    default_args = {
        "get_children": True,
        "sensitive": "true",
        "heterogeneous": "false",
        "phenDB": True,
        "faprotax": True,
        "pathway_complement": True,
        "seed_scores": False,
        "manta": False
    }
    for i in list(default_args.keys()):
        if i not in list(cfg.keys()):
            if i == "delimiter":
                return "You have not provided a taxonomy delimiter. microbetag will exit. Please set one and try again."
            cfg[i] = default_args[i]

        elif i == "sensitive" or i == "heterogeneous":
            if cfg[i]:
                cfg[i] = "true"
            else:
                cfg[i] = "false"

    if metadata_file is None:
        metadata_file = "false"

    # Abundance table preprocess
    if abundance_table_file:

        logging.info("STEP: Assign NCBI Tax Id and GTDB reference genomes".center(80, "*"))

        # Abundance table as a pd dataframe
        abundance_table = is_tab_separated(abundance_table_file, TAX_COL)

        if not edge_list_file:

            logging.info(
                "The user has not provided an edge list. microbetag will build one using FlashWeaeve.")
            if not os.path.exists(flashweave_output_dir):
                os.mkdir(flashweave_output_dir)

            # ext remains a pd dataframe
            ext = ensure_flashweave_format(abundance_table, TAX_COL, SEQ_COL, flashweave_output_dir)

        else:
            ext = abundance_table.copy()
            ext["microbetag_id"] = abundance_table[SEQ_COL]

        # Map taxonomies to ontology ids
        logging.info(
            "Get the NCBI Taxonomy id for those OTUs that have been assigned either at the species, the genus or the family level.")

        if cfg["get_children"]:
            seqID_taxid_level_repr_genome, repr_genomes_present, children_df = map_seq_to_ncbi_tax_level_and_id(
                ext,
                TAX_COL,
                SEQ_COL,
                cfg["taxonomy"],
                cfg["delimiter"],
                cfg["get_children"]
            )

        else:
            seqID_taxid_level_repr_genome, repr_genomes_present = map_seq_to_ncbi_tax_level_and_id(
                ext,
                TAX_COL,
                SEQ_COL,
                cfg["taxonomy"],
                cfg["delimiter"],
                cfg["get_children"]
            )
            children_df = pd.DataFrame()

        children_genomes = set(children_df.explode("gtdb_gen_repr")["gtdb_gen_repr"].to_list())

        # [ATTENTION!] seqID not exploded!
        seqID_taxid_level_repr_genome["taxonomy"] = ext["taxonomy"]
        exploded_seqID_taxid_level_repr_genome = seqID_taxid_level_repr_genome.explode("gtdb_gen_repr")

        seqID_taxid_level_repr_genome.to_csv(
            os.path.join(
                out_dir,
                "taxa_with_ndbiId_and_repr_gtdbIds.csv"),
            "\t")

    # Get co-occurrence network
    if not edge_list_file:

        logging.info("STEP: Build co-occurrence network".center(80, "*"))

        flashweave_params = [
            "julia", FLASHWEAVE_SCRIPT, flashweave_output_dir, flashweave_tmp_input, str(cfg["sensitive"]), str(cfg["heterogeneous"]), metadata_file
        ]
        flashweave_command = " ".join(flashweave_params)
        logging.info("Run FlashWeave")
        os.system(flashweave_command)

        # Taxa pairs as NCBI Tax ids
        logging.info("Map your edge list to NCBI Tax ids and keep only associations that both correspond to a such.")
        edge_list = edge_list_of_ncbi_ids(flashweave_edgelist)

        """Example:
        [{'node_a': 'microbetag_17', 'ncbi_tax_id_node_a': 77133, 'gtdb_gen_repr_node_a': 'GCA_903925685.1', 'ncbi_tax_level_node_a': 'mspecies',
         'node_b': 'microbetag_21', 'ncbi_tax_id_node_b': 136703, 'gtdb_gen_repr_node_b': nan, 'ncbi_tax_level_node_b': 'species'},.. {..}] """
        edgelist_as_a_list_of_dicts = edge_list.applymap(lambda x: str(x) if pd.notna(x) else 'null').to_dict(orient="records")

        # Save the edgelist returned from flashweave to a file
        edge_list.to_csv("edgelist.csv", sep="\t")

        # Make a set with all unique sequence ids that are nodes in the network
        all_nodes_ids = set(pd.concat([edge_list["node_a"], edge_list["node_b"]]).tolist())

    # Phen annotations
    if cfg["phenDB"]:
        """
        NOTE: Starting from Python 3.7, dictionaries in Python maintain the order of insertion.
        So, if you create a dictionary d using elements from a list as keys and then use list(d.keys()),
        the list will have the same order as the keys were inserted into the dictionary.
        """
        logging.info("STEP: PhenDB ".center(80, "*"))
        # Get phen traits for each GTDB genome present on your table
        feats = get_column_names("phenDB")
        feats.insert(0, "gtdb_id"); feats.insert(1, "NCBI_ID"); feats.remove("gtdbId")
        traits = []
        for gtdb_id in set(repr_genomes_present):
            phen_traits_dict = get_phendb_traits(gtdb_id)
            if phen_traits_dict != 0:
                if gtdb_id not in children_genomes:
                    ncbiId = exploded_seqID_taxid_level_repr_genome.loc[
                        exploded_seqID_taxid_level_repr_genome["gtdb_gen_repr"] == gtdb_id, "ncbi_tax_id"
                    ].to_list()[0]
                else:
                    # [NOTE] Remember that the actual NCBI Tax Id of the genome in this case, is the child, not the one written in the file.
                    ncbiId = children_df.explode("gtdb_gen_repr")[
                        children_df.explode("gtdb_gen_repr")["gtdb_gen_repr"] == gtdb_id]["parent_ncbi_tax_id"].to_list()[0]
                ncbiId = str(ncbiId)
                q = list(phen_traits_dict.values())
                q.insert(1, ncbiId)
                q = [str(x) for x in q]
                traits.append(q)
            else:
                print("Genome:", gtdb_id, "is not present in the phenDB version of microbetag.")
        # Save phen traits as a .tsv file
        if not os.path.exists(phen_output_dir):
            os.mkdir(phen_output_dir)
        outfile = os.path.join(phen_output_dir, "phen_traits.tsv")
        export_phen_traits_to_file(
            column_names=feats,
            rows=traits,
            filename=outfile)

    # FAPROTAX
    if cfg["faprotax"] and abundance_table_file:

        logging.info("STEP: FAPROTAX database oriented analaysis".center(80, "*"))

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

    # Get species/strain to species/strain associations
    if cfg["pathway_complement"] or cfg["seed_scores"]:

        logging.info("""For the pathway complementarity and the seed scores modules, we focus only on
            to species/strain to species/strain level associations of the network. Let's find those pairs!""")

        """
        species_level_associations: set of sets of NCBI ids
        genomes_of_species_nodes: a dictionary with the genomes assigned (value; type: list) to each ncbi id (key)
        """
        species_level_associations, seqIds_of_species_level_associations, genomes_of_species_nodes, genomes_of_children = export_species_level_associations(
            edgelist_as_a_list_of_dicts,
            seqID_taxid_level_repr_genome,
            children_df
        )

    # Pathway complementarity
    if cfg["pathway_complement"]:

        logging.info("STEP: Pathway complementarity".center(80, "*"))

        if edge_list_file:
            edge_list = edge_list_file.copy()

        ncbi_id_pairs_with_complements = get_complements_of_list_of_pair_of_ncbiIds(species_level_associations, genomes_of_species_nodes)

        if not os.path.exists(pathway_complementarity_dir):
            os.mkdir(pathway_complementarity_dir)

        with open(os.path.join(pathway_complementarity_dir, "complements.json"), "w") as f:
            ncbi_id_pairs_with_complements = {tuple_to_str(key): value for key, value in ncbi_id_pairs_with_complements.items()}
            ncbi_id_pairs_with_complements_str = convert_to_json_serializable(ncbi_id_pairs_with_complements)  # convert_tuples_to_strings
            json.dump(ncbi_id_pairs_with_complements_str, f)

        logging.info("Pathway complementarity has been completed successfully.".center(80, "*"))

    # Seed - based scores
    if cfg["seed_scores"]:

        logging.info("""STEP: Seed scores based on draft genome-scale reconstructions.""".center(80, "*"))

        ncbi_id_pairs_with_seed_scores = get_seed_scores_for_list_of_pair_of_ncbiIds(species_level_associations, genomes_of_species_nodes)

        if not os.path.exists(seed_scores_dir):
            os.mkdir(seed_scores_dir)

        ncbi_id_pairs_with_seed_scores = remove_query(ncbi_id_pairs_with_seed_scores)
        ncbi_id_pairs_with_seed_scores = {tuple_to_str(key): value for key, value in ncbi_id_pairs_with_seed_scores.items()}
        ncbi_id_pairs_with_seed_scores_str = convert_to_json_serializable(ncbi_id_pairs_with_seed_scores)  # convert_tuples_to_strings

        with open(os.path.join(seed_scores_dir, "seed_scores.json"), "w") as f:
            json.dump(ncbi_id_pairs_with_seed_scores_str, f)

        logging.info("Seed scores have been assigned successfully.". center(80, "*"))

    # Run manta
    if cfg["manta"]:

        logging.info("""STEP: network clustering using manta and the abundance table""".center(80, "*"))

        # Build base network; no annotations added
        base_network = build_base_graph(edgelist_as_a_list_of_dicts, seqID_taxid_level_repr_genome, cfg)
        base_network = convert_to_json_serializable(base_network)
        base_network_file = os.path.join(out_dir, "basenet.cyjs")
        with open(base_network_file, "w") as f:
            json.dump(base_network, f, indent=4)

        logging.info("Base network has been built and saved.")

        # Build the manta command
        manta_output_file = "/".join([out_dir, 'manta_annotated'])
        manta_params = [
            "manta",
            "-i", base_network_file,
            "-f", "cyjs",
            "-o", manta_output_file,
            "--layout"
        ]
        manta_command = " ".join(manta_params)

        # Run manta
        if os.system(manta_command) != 0:
            logging.error("""manta failed. Check network format""")
            sys.exit(0)

    # Build output; an annotated graph
    logging.info("""STEP: Constructing the annotated network""".center(80, "*"))
    annotated_network = build_cx_annotated_graph(
        edge_list,
        edgelist_as_a_list_of_dicts,
        seqID_taxid_level_repr_genome,
        cfg,
        out_dir,
        genomes_of_children
    )
    # Save annotated network to file
    net_file = "/".join([out_dir, 'annotated.cx'])
    with open(net_file, "w") as f:
        annotated_network2file = convert_to_json_serializable(annotated_network)
        json.dump(annotated_network2file, f)

    return annotated_network
