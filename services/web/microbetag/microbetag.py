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
    # if cfg["taxonomy"] == "GTDB":
    #     gtdb_accession_ncbi_tax_id = os.path.join(
    #         MAPPINGS, "gtdbSpecies2ncbiId2accession.tsv")
    # else:
    #     gtdb_accession_ncbi_tax_id = os.path.join(
    #         MAPPINGS, "species2ncbiId2accession.tsv")

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

    # Load vars
    flashweave_output_dir = os.path.join(out_dir, "flashweave")
    flashweave_tmp_input = os.path.join(flashweave_output_dir, "abundance_table_flashweave_format.tsv")
    flashweave_edgelist = os.path.join(flashweave_output_dir,"network_output.edgelist")

    faprotax_output_dir = os.path.join(out_dir, "faprotax")
    faprotax_funct_table = os.path.join(faprotax_output_dir, "functional_otu_table.tsv")
    faprotax_sub_tables = os.path.join(faprotax_output_dir, "sub_tables")

    phen_output_dir = os.path.join(out_dir, "phen_predictions")
    pathway_complementarity_dir = os.path.join(out_dir, "path_compl")
    seed_scores_dir = os.path.join(out_dir, "seed_scores")

    if cfg["delimiter"]:
        tax_delim = cfg["delimiter"]

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
        try:
            seqID_taxid_level_repr_genome, repr_genomes_present = map_seq_to_ncbi_tax_level_and_id(ext, TAX_COL, SEQ_COL, cfg["taxonomy"], tax_delim)

        except BaseException:
            logging.error("""microbetag was not able to map your table's taxonomies to NCBI and GTDB ids.
                           Check on how you describe your table.
                           Also check whether you have set the edge_list_file parameter in the config file;
                           if there is not an edge list, leave it blank.""")
            sys.exit(0)

        seqID_taxid_level_repr_genome.to_csv(
            os.path.join(
                out_dir,
                "taxa_with_ndbiId_and_repr_gtdbIds.csv"),
            "\t")

    # Get co-occurrence network
    if not edge_list_file:

        logging.info("STEP: Build co-occurrence network".center(80, "*"))

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
        #   'node_b': 'microbetag_21', 'ncbi_tax_id_node_b': 136703, 'gtdb_gen_repr_node_b': nan, 'ncbi_tax_level_node_b': 'species'},.. {..}]
        edgelist_as_a_list_of_dicts = edge_list.applymap(lambda x: str(x) if pd.notna(x) else 'null').to_dict(orient="records")

        edge_list.to_csv("edgelist.csv", sep="\t")

        microb_id_taxonomy = ext[["microbetag_id", "taxonomy"]]

        base_network = build_base_graph(edgelist_as_a_list_of_dicts, microb_id_taxonomy)

        with open(os.path.join(out_dir, "basenet.json"), "w") as f:
            json.dump(base_network, f)

        logging.info("Base network has been built and saved.")

    # Phen annotations

    if cfg["phenDB"]:

        logging.info("STEP: PhenDB ".center(80, "*"))

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
                    to species/strain to species/strain level associations of the network
                    Let's find those pairs!""")

        species_level_associations, genomes_of_species_nodes = export_species_level_associations(edgelist_as_a_list_of_dicts)

    # Pathway complementarity
    if cfg["pathway_complement"]:

        logging.info("STEP: Pathway complementarity".center(80, "*"))

        if edge_list_file:
            edge_list = edge_list_file.copy()

        ncbi_id_pairs_with_complements = get_complements_of_list_of_pair_of_ncbiIds(species_level_associations, genomes_of_species_nodes)
        os.mkdir(pathway_complementarity_dir)

        with open(os.path.join(pathway_complementarity_dir, "complements.json"), "w") as f:
            ncbi_id_pairs_with_complements = {tuple_to_str(key): value for key, value in ncbi_id_pairs_with_complements.items()}
            json.dump(ncbi_id_pairs_with_complements, f)

        logging.info("Pathway complementarity has been completed successfully.".center(80, "*"))

    # Seed - based scores
    if cfg["seed_scores"]:
        logging.info("""STEP: Competition and complementarity scores between associated species/strains based on their drafte genome-scale
                   reconstructions and the Seed appraoch.""".center(80, "*"))

        ncbi_id_pairs_with_seed_scores = get_seed_scores_for_list_of_pair_of_ncbiIds(species_level_associations, genomes_of_species_nodes)
        os.mkdir(seed_scores_dir)

        print(ncbi_id_pairs_with_seed_scores)

        with open(os.path.join(seed_scores_dir, "seed_scores.json"), "w") as f:
            # ncbi_id_pairs_with_seed_scores = {tuple_to_str(key): value for key, value in ncbi_id_pairs_with_seed_scores.items()}
            json.dump(ncbi_id_pairs_with_seed_scores, f)

        logging.info("Seed scores have been assigned successfully.". center(80, "*"))

    if cfg["manta"]:
        if not edge_list_file:
            logging.info("""STEP: network clustering using manta and the abundance table""".center(80, "*"))

            # fix input data  
            os.system("sed  '1,2d'  flashweave/network_detailed_output.edgelist  > edgelist_for_manta.txt")

            # 

            manta_params = [
                "manta",
                "-i", "edgelist_for_manta.txt"
            ]


    # Build output; an annotated graph
    logging.info("""STEP: Constructing the annotated network""")
    annotated_network = annotate_network(base_network, cfg, seqID_taxid_level_repr_genome, out_dir)


    with open(os.path.join(out_dir, "annotated_network.json"), "w") as f:
        json.dump(annotated_network, f)

    return annotated_network
