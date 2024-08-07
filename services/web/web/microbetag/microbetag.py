#!/usr/bin/env python3

"""
Aim:
    This is the main of the microbetag workflow.
    It gets as input user's config file as a dictionary (cfg) and their corresponding data (abundance table, metadata file, network)
    and returns an annotated network file in .cx format

Author:
    Haris Zafeiropoulos

Licensed under GNU LGPL.3, see LICENCE file
"""
from microbetag.scripts.utils import *
from microbetag.scripts.db_functions import *
from microbetag.scripts.logger import *
from microbetag.scripts.variables import *
from microbetag.scripts.buid_cx_annotated_graph import *
import os
import subprocess
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


def main(out_dir, cfg, abundance_table_file=None, edge_list_file=None, metadata_file=None, logger=None):
    """
    microbetagApp main routine
    Gets an abundance table in .tsv format and optionally a co-occurrence netwok also in .tsv format.
    Returns an annotated network in .cx format, based on the annotations asked through the confige var (cfg)
    """
    # Using FileHandler writing log to file
    if logger is None:
        # Set logger
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter(
            '%(asctime)s [%(levelname)-8s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')
        logfile = os.path.join(out_dir, "log.txt")
        open(logfile, "w")
        fh = logging.FileHandler(logfile); fh.setLevel(logging.DEBUG); fh.setFormatter(formatter)
        # Using StreamHandler writing to console, error and the two Handlers
        ch = logging.StreamHandler(); ch.setLevel(logging.INFO); ch.setFormatter(formatter)
        eh = logging.StreamHandler(); eh.setLevel(logging.ERROR); eh.setFormatter(formatter)
        logger.addHandler(ch); logger.addHandler(fh); logger.addHandler(eh)

    # Load vars
    flashweave_output_dir = os.path.join(out_dir, "flashweave")
    flashweave_tmp_input = os.path.join(flashweave_output_dir, "abd_tab_flashweave_format.tsv")
    flashweave_edgelist = os.path.join(flashweave_output_dir, "network_output.edgelist")

    faprotax_output_dir = os.path.join(out_dir, "faprotax")
    faprotax_funct_table = os.path.join(faprotax_output_dir, "functional_otu_table.tsv")
    faprotax_sub_tables = os.path.join(faprotax_output_dir, "sub_tables")

    phen_output_dir = os.path.join(out_dir, "phen_predictions")
    pathway_complementarity_dir = os.path.join(out_dir, "path_compl")
    seed_scores_dir = os.path.join(out_dir, "seed_scores")

    # Fix arguments
    default_args = {
        "get_children": False,
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
                logger.error(
                    "You have not provided a taxonomy delimiter. microbetag will exit. Please set one and try again."
                )
            cfg[i] = default_args[i]

        elif i == "sensitive" or i == "heterogeneous":
            if cfg[i]:
                cfg[i] = "true"
            else:
                cfg[i] = "false"

    if metadata_file is None:
        metadata_file = "false"
        metadata = "false"
        env_edges = []
    else:
        metadata = "true"

    # Abundance table preprocess
    if abundance_table_file is None:
        e = """microbetag requires an abundance table."""
        logger.error(e)
        raise FileNotFoundError(e)

    """
    Parse abundance table to match taxa to microbetag genomes
    """
    logger.info("STEP: Assign NCBI Tax Id and GTDB reference genomes".center(80, "*"))

    # Abundance table as a pd dataframe
    abundance_table = is_tab_separated(abundance_table_file, TAX_COL)

    if not edge_list_file:
        logger.info("The user has not provided an edge list. microbetag will build one using FlashWeaeve.")
        if not os.path.exists(flashweave_output_dir):
            os.mkdir(flashweave_output_dir)
        # ext remains a pd dataframe
        ext = ensure_flashweave_format(abundance_table, TAX_COL, SEQ_COL, flashweave_output_dir)

    else:
        ext = abundance_table.copy()
        ext["microbetag_id"] = abundance_table[SEQ_COL]

    # Map taxonomies to ontology ids
    logger.info(
        """Get the NCBI Taxonomy id for those OTUs that have been assigned either at the species, the genus or the family level."""
    )

    if cfg["get_children"]:
        logger.info(ext)
        seqID_taxid_level_repr_genome, repr_genomes_present, children_df = map_seq_to_ncbi_tax_level_and_id(
            ext,
            TAX_COL,
            SEQ_COL,
            cfg["taxonomy"],
            cfg["delimiter"],
            cfg["get_children"]
        )
        children_genomes = set(children_df.explode("gtdb_gen_repr")["gtdb_gen_repr"].to_list())

    else:
        seqID_taxid_level_repr_genome, repr_genomes_present = map_seq_to_ncbi_tax_level_and_id(
            ext,
            TAX_COL,
            SEQ_COL,
            cfg["taxonomy"],
            cfg["delimiter"],
            cfg["get_children"]
        )

    # Save annotated taxonomies to file
    seqID_taxid_level_repr_genome["taxonomy"] = ext["taxonomy"]
    seqID_taxid_level_repr_genome.to_csv(
        os.path.join(
            out_dir,
            "taxa_with_ndbiId_and_repr_gtdbIds.csv"
        ), "\t")

    # IMPORTANT! Explode df so 1 row for 1 genome
    exploded_seqID_taxid_level_repr_genome = seqID_taxid_level_repr_genome.explode("gtdb_gen_repr")

    # Get co-occurrence network
    if not edge_list_file:

        logger.info("STEP: Build co-occurrence network".center(80, "*"))

        flashweave_params = [
            "julia", FLASHWEAVE_SCRIPT, flashweave_output_dir, flashweave_tmp_input, str(cfg["sensitive"]), str(cfg["heterogeneous"]), metadata, metadata_file
        ]
        flashweave_command = " ".join(flashweave_params)

        logger.info("Run FlashWeave")
        if os.system(flashweave_command) != 0:
            e = """ \
                FlashWeave failed.
                Please check on the content of your abundance table and the format of your metadata file if one provided.
                We suggest you use FlashWeave or any other software to get the network and then perform microbetag providing the network returned.
                You may find an implementation of FlashWeave in a Docker image in microbetag's preprocessing image:
                https://hub.docker.com/r/hariszaf/microbetag_prep
            """
            logger.error(e)
            raise ValueError(e)

        # Taxa pairs as NCBI Tax ids
        logger.info("""Map your edge list to NCBI Tax ids and keep only associations that both correspond to a such.""")
        try:
            if metadata == "true":
                edge_list, env_edges = edge_list_of_ncbi_ids(flashweave_edgelist, metadata_file=metadata_file)
            else:
                edge_list = edge_list_of_ncbi_ids(flashweave_edgelist)
        except ValueError:
            e = """ \
                No edges in the cooccurrence network produced. Please check your abundance table and the FlashWeave settings used.
                Most likely a co-occurrence network cannot be built using your abundance table and this set of parameters.
                We suggest you use FlashWeave or any other software to get the network and then perform microbetag providing the network returned.
                You may find an implementation of FlashWeave in a Docker image in microbetag's preprocessing image:
                https://hub.docker.com/r/hariszaf/microbetag_prep
            """
            logger.error(e)
            raise ValueError(e)

        # Save the edgelist returned from flashweave to a file
        edge_list.to_csv("edgelist.tsv", sep="\t")

    else:
        logger.info("STEP: Load user's co-occurrence network".center(80, "*"))
        try:
            edge_list = edge_list_of_ncbi_ids(edge_list_file)
        except ValueError:
            e = """The network file provided is either empty or of bad format and could not be loaded."""
            logger.error(e)
            raise ValueError(e)

        all_seqids_on_edgelist = list(edge_list["node_a"].values) + list(edge_list["node_b"].values)
        all_seqids_on_edgelist = set(all_seqids_on_edgelist)
        all_seqids_in_abundance_table = set(seqID_taxid_level_repr_genome["microbetag_id"].values)
        if not all_seqids_on_edgelist.issubset(all_seqids_in_abundance_table):
            e = """Not all sequence identifiers in the co-occurrence network provided are also menbers of abundance table!"""
            logger.warn(e)
            raise ValueError(e)

    """
    Example edgelist_as_a_list_of_dicts:
    [ {'node_a': 'microbetag_17', 'ncbi_tax_id_node_a': 77133, 'gtdb_gen_repr_node_a': 'GCA_903925685.1', 'ncbi_tax_level_node_a': 'mspecies',
      'node_b': 'microbetag_21', 'ncbi_tax_id_node_b': 136703, 'gtdb_gen_repr_node_b': nan, 'ncbi_tax_level_node_b': 'species'} ]
    """
    edgelist_as_a_list_of_dicts = edge_list.map(lambda x: str(x) if pd.notna(x) else 'null').to_dict(orient="records")

    if not isinstance(env_edges, list):
        # [{'node_a': 'ENV4_A', 'node_b': 'ENV4_C', 'score': -0.774596631526947}, .. ]
        env_edges = env_edges.to_dict(orient="records")

    # # Make a set with all unique sequence ids that are nodes in the network
    # all_nodes_ids = set(pd.concat([edge_list["node_a"], edge_list["node_b"]]).tolist())

    # Phen annotations
    if cfg["phenDB"]:
        """
        NOTE: Starting from Python 3.7, dictionaries in Python maintain the order of insertion.
        So, if you create a dictionary d using elements from a list as keys and then use list(d.keys()),
        the list will have the same order as the keys were inserted into the dictionary.
        """
        logger.info("STEP: PhenDB ".center(80, "*"))

        # Get phen traits for each GTDB genome present on your table
        feats = get_column_names("phenDB")
        feats.insert(0, "gtdb_id"); feats.insert(1, "NCBI_ID"); feats.remove("gtdbId")
        traits = []

        for gtdb_id in set(repr_genomes_present):

            phen_traits_dict = get_phendb_traits(gtdb_id)

            if phen_traits_dict != 0:

                if cfg["get_children"]:

                    if gtdb_id in children_genomes:
                        # print("this is a child genome:", gtdb_id)
                        # [NOTE] Remember that the actual NCBI Tax Id of the genome in this case, is the child, not the one written in the file.
                        ncbiId = children_df.explode("gtdb_gen_repr")[
                            children_df.explode("gtdb_gen_repr")["gtdb_gen_repr"] == gtdb_id
                        ]["parent_ncbi_tax_id"].to_list()[0]

                    else:
                        # print("not a child genome:", gtdb_id)
                        ncbiId = exploded_seqID_taxid_level_repr_genome.loc[
                            exploded_seqID_taxid_level_repr_genome["gtdb_gen_repr"] == gtdb_id, "ncbi_tax_id"
                        ].to_list()[0]
                else:
                    ncbiId = exploded_seqID_taxid_level_repr_genome.loc[
                        exploded_seqID_taxid_level_repr_genome["gtdb_gen_repr"] == gtdb_id, "ncbi_tax_id"
                    ].to_list()[0]

                ncbiId = str(ncbiId)
                q = list(phen_traits_dict.values())
                q.insert(1, ncbiId)
                q = [str(x) for x in q]
                traits.append(q)

            else:
                logger.info(" ".join(["Genome:", gtdb_id, "is not present in the phenDB version of microbetag."]))

        # Save phen traits as a .tsv file
        if not os.path.exists(phen_output_dir):
            os.mkdir(phen_output_dir)
        outfile = os.path.join(phen_output_dir, "phen_traits.tsv")
        export_phen_traits_to_file(
            column_names=feats,
            rows=traits,
            filename=outfile
        )

    # FAPROTAX
    if cfg["faprotax"] and abundance_table_file:

        logger.info("STEP: FAPROTAX database oriented analaysis".center(80, "*"))

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
        process = subprocess.Popen(faprotax_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            e = """ \
                Something went wrong when running the FAPROTAX analysis!
                Check how you describe your input table.
                Also, please make sure you have set the column_names_are_in parameter properly.
            """
            logger.error(stderr.decode())
            logger.error()
            raise ValueError(e)

    # Get species/strain to species/strain associations
    if cfg["pathway_complement"] or cfg["seed_scores"]:

        logger.info(
            """STEP: Extract edges where both nodes represent taxa annotated at the species/strain level.""".center(80, "*")
        )
        """
        species_level_associations: set of sets of NCBI ids
        genomes_of_species_nodes: a dictionary with the genomes assigned (value; type: list) to each ncbi id (key)
        """
        logger.info(edgelist_as_a_list_of_dicts)
        if cfg["get_children"]:
            species_level_associations, genomes_of_species_nodes, child_genomes_in_net = export_species_level_associations(
                edgelist_as_a_list_of_dicts,
                seqID_taxid_level_repr_genome,
                children_df
            )
        else:
            species_level_associations, genomes_of_species_nodes = export_species_level_associations(
                edgelist_as_a_list_of_dicts,
                seqID_taxid_level_repr_genome,
            )

    # Pathway complementarity
    if cfg["pathway_complement"]:

        logger.info("STEP: Pathway complementarity".center(80, "*"))
        try:
            ncbi_id_pairs_with_complements = get_complements_of_list_of_pair_of_ncbiIds(species_level_associations, genomes_of_species_nodes)
        except ValueError:
            e = """ \
                microbetag failed to get pathway complementarities using your abundance table.
                This is probably due to potential bugs on our side.
                Please consider contacting us at haris.zafeiropoulos@kuleuven.be sharing your abundance table.
            """

        if not os.path.exists(pathway_complementarity_dir):
            os.mkdir(pathway_complementarity_dir)

        with open(os.path.join(pathway_complementarity_dir, "complements.json"), "w") as f:
            ncbi_id_pairs_with_complements = {tuple_to_str(key): value for key, value in ncbi_id_pairs_with_complements.items()}
            ncbi_id_pairs_with_complements_str = convert_to_json_serializable(ncbi_id_pairs_with_complements)
            json.dump(convert_tuples_to_strings(ncbi_id_pairs_with_complements_str), f)

        logger.info("Pathway complementarity has been completed successfully.".center(80, "*"))

    # Seed - based scores
    if cfg["seed_scores"]:

        logger.info("""STEP: Seed scores based on draft genome-scale reconstructions.""".center(80, "*"))

        # Get all PATRIC ids corresponding to the relative_genomes
        flat_list = [item for sublist in list(genomes_of_species_nodes.values()) for item in sublist]
        patricIds = get_patric_id_of_gc_accession_list(flat_list)

        # Get Seed scores
        ncbi_id_pairs_with_seed_scores = get_seed_scores_and_complements_for_list_of_pairs_of_ncbiIds(
            species_level_associations,
            genomes_of_species_nodes,
            patricIds
        )

        if not os.path.exists(seed_scores_dir):
            os.mkdir(seed_scores_dir)

        ncbi_id_pairs_with_seed_scores = remove_query(ncbi_id_pairs_with_seed_scores)
        ncbi_id_pairs_with_seed_scores = {tuple_to_str(key): value for key, value in ncbi_id_pairs_with_seed_scores.items()}
        ncbi_id_pairs_with_seed_scores_str = convert_to_json_serializable(ncbi_id_pairs_with_seed_scores)  # convert_tuples_to_strings
        ncbi_id_pairs_witt_seed_scores_js = convert_tuples_to_strings(ncbi_id_pairs_with_seed_scores_str)
        with open(os.path.join(seed_scores_dir, "seed_scores.json"), "w") as f:
            json.dump(ncbi_id_pairs_witt_seed_scores_js, f)

        logger.info("Seed scores have been assigned successfully.". center(80, "*"))

    # Run manta
    if cfg["manta"]:

        import time
        logger.info("""STEP: network clustering using manta and the abundance table""".center(80, "*"))

        # Build base network; no annotations added
        base_network = build_base_graph(edgelist_as_a_list_of_dicts, seqID_taxid_level_repr_genome, cfg)
        base_network = convert_to_json_serializable(base_network)
        base_network_file = os.path.join(out_dir, "basenet.cyjs")
        with open(base_network_file, "w") as f:
            json.dump(base_network, f, indent=4)

        logger.info("Base network has been built and saved.")

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
        m1 = time.time()
        if os.system(manta_command) != 0:
            e = """\
                The manta clustering algorithm failed.
                Most likely this is because clusters could not be grouped based on the provided network and the parameters setup of manta.
                Yet, microbetag will continue to the following steps without considering for clusters.
            """
            logger.warn(e)
            # Changed cfg so the buld_cx_annotated_graph function will not fail.
            cfg["manta"] = False
        else:
            logger.info("""manta ran fine.""")
        m2 = time.time()
        time = " ".join(["Network clustering with manta took:", str(m2 - m1), "sec"])
        logger.info(time)

    # Build output; an annotated graph
    logger.info("""STEP: Constructing the annotated network""".center(80, "*"))
    try:
        if cfg["get_children"]:

            if len(child_genomes_in_net) == 0:
                print("No network node has children genomes.")
                child_genomes_in_net = None

            annotated_network = build_cx_annotated_graph(
                edge_list,
                edgelist_as_a_list_of_dicts,
                seqID_taxid_level_repr_genome,
                env_edges,
                cfg,
                out_dir,
                child_genomes_in_net
            )

        else:
            annotated_network = build_cx_annotated_graph(
                edge_list,
                edgelist_as_a_list_of_dicts,
                seqID_taxid_level_repr_genome,
                env_edges,
                cfg,
                out_dir
            )

    except ValueError:
        e = """\
            microbetag failed building the annotated network.
            This is probably the result of a bug.
            Please consider contacting us at haris.zafeiropoulos@kuleuven.be sharing your abundance table.
        """
        logger.error(e)
        raise ValueError(e)

    # Save annotated network to file
    net_file = "/".join([out_dir, 'annotated.cx'])
    with open(net_file, "w") as f:
        annotated_network2file = convert_to_json_serializable(annotated_network)
        json.dump(annotated_network2file, f)

    return annotated_network
