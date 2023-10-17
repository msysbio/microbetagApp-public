"""
Parse user's input OTU table
"""
from .db_functions import *
from .variables import *
from fuzzywuzzy import fuzz
import pandas as pd
import numpy as np
import logging
import os
import time
import re
import json
import multiprocessing
import networkx as nx


def init_run():
    """
    fix it so we have a route for cytoscape app
    another one for terminal where user gets a zip file
    """


def annotate_network(base_network, config, ids_map, out_dir):
    """
    Takes as input a base network and adds microbetag annotations from the various channels of information.
    ids_map: sequence id (ASV, OTU id etc), NCBI Tax id and level and last GTDB representative genome
    """
    annotated_network = base_network.copy()

    if config["phenDB"]:

        df_traits_table = pd.read_csv(os.path.join(out_dir, "phen_predictions/phen_traits.tsv"), sep="\t")
        for node in annotated_network["elements"]["nodes"]:
            ncbiId = node["data"]["NCBI-Tax-Id"]
            if ncbiId == "null":
                continue
            if int(ncbiId) in df_traits_table['NCBI_ID'].values:
                node_traits = df_traits_table.iloc[df_traits_table["NCBI_ID"].values == int(ncbiId)].to_dict()
                node_traits = {key: list(value.values())[0] for key, value in node_traits.items()}
                node["data"]["phenotypic-traits"] = node_traits

    if config["faprotax"]:

        assignments_per_seqId = seqId_faprotax_functions_assignment(os.path.join(out_dir, "faprotax/sub_tables/"))
        assignments_per_seqId = {key: [value.rstrip('.txt').replace('_', ' ') for value in values] for key, values in assignments_per_seqId.items()}
        for node in annotated_network["elements"]["nodes"]:
            for seqId, assignments in assignments_per_seqId.items():
                if seqId == node["data"]["id"]:
                    node["data"]["faprotax-assignments"] = assignments

    if config["pathway_complement"]:

        complements_dict = json.load(open(os.path.join(out_dir, "path_compl/complements.json")))

        for pair, complements in complements_dict.items():
            ncbi_a, ncbi_b, = pair.split(",")
            for edge in annotated_network["elements"]["edges"]:
                if ncbi_a == edge["data"]["source-ncbi-tax-id"] and ncbi_b == edge["data"]["target-ncbi-tax-id"]:
                    edge["source-to-target-complements"] = complements

                if ncbi_a == edge["data"]["target-ncbi-tax-id"] and ncbi_b == edge["data"]["source-ncbi-tax-id"]:
                    edge["target-to-source-complements"] = complements

    if config["seed_scores"]:
        seed_scores_dict = json.load(open(os.path.join(out_dir, "seed_scores/seed_scores.json")))

        for pair, seed_scores in seed_scores_dict.items():
            ncbi_a, ncbi_b, = pair.split(",")
            for edge in annotated_network["elements"]["edges"]:
                if ncbi_a == edge["data"]["source-ncbi-tax-id"] and ncbi_b == edge["data"]["target-ncbi-tax-id"]:
                    edge["seed-scores"] = seed_scores

    annotated_network = convert_to_float(annotated_network)

    return annotated_network


def export_species_level_associations(edgelist_as_a_list_of_dicts):
    """
    This functions gets as input the edges of the network as a list of dictionaries
    checks the ncbi_tax_level of each node and in case where for both nodes it is species or strain
    gets their corresponding GTDB genomes on microbetagDB
    Returns:
    a. pairs_of_interest: a set with the ncbi tax ids of the linking nodes
    b. related_genomes: a dictionary with the genomes assigned (values) to each ncbi id (key)
    c. parent_children_ncbiIds_present: {}
    """

    set_of_ncbiids_of_interest = set()
    list_of_non_usefules_ids = ["77133"]
    pairs_of_interest = set()

    for pair in edgelist_as_a_list_of_dicts:
        taxon_a = str(pair["ncbi_tax_id_node_a"]).split(".")[0]
        taxon_b = str(pair["ncbi_tax_id_node_b"]).split(".")[0]
        if taxon_a in list_of_non_usefules_ids or taxon_b in list_of_non_usefules_ids:
            continue

        # We keep as a pair both a->b and b->a association
        if pair["ncbi_tax_level_node_a"] == pair["ncbi_tax_level_node_b"] == "mspecies":
            set_of_ncbiids_of_interest.add(taxon_a)
            set_of_ncbiids_of_interest.add(taxon_b)
            pairs_of_interest.add((taxon_a, taxon_b))
            pairs_of_interest.add((taxon_b, taxon_a))

    # Build nodes dict
    """
    ATTENTION! We have noticed that in some cases, we have strain genomes but no genomes of their corresponding genomes.
    On top of that, a species may have several genomes so even if there is a genome for the species, it's good to get
    what its strains suggest.
    """
    ncbi_nodes_dict = get_ncbi_nodes_dict()
    ncbiId_to_genomesIds_queries = []
    parent_children_ncbiIds_present = {}
    for ncbiId in set_of_ncbiids_of_interest:
        nids = []
        if ncbiId in ncbi_nodes_dict:
            nids = ncbi_nodes_dict[ncbiId]
            parent_children_ncbiIds_present[ncbiId] = nids
        nids.append(ncbiId)
        for nid in nids:
            query = "".join(["SELECT genomeId from genome2taxNcbiId where ncbiTaxId = ", str(nid), ";"])
            ncbiId_to_genomesIds_queries.append(query)

    # Execute the queries
    genomes = []
    cnx, cursor = create_cursor()
    for query in ncbiId_to_genomesIds_queries:
        cursor.execute(query)
        query_results = cursor.fetchall()
        genomes.append({query.split(" ")[-1][:-1]: [x[0] for x in query_results if x[0].startswith("GCA_") or x[0].startswith("GCF_")]})

    related_genomes = {key: value for dictionary in genomes for key, value in dictionary.items()}

    return pairs_of_interest, related_genomes, parent_children_ncbiIds_present


def count_comment_lines(my_abundance_table, my_taxonomy_column):
    """
    Get the number of rows of the OTU table that should be skipped
    """
    skip_rows = 0
    with open(my_abundance_table, 'r') as f:
        for line in f:
            if line.startswith('#') and my_taxonomy_column not in line:
                skip_rows += 1
            elif my_taxonomy_column in line:
                line = line.split("\t")
                line[-1] = line[-1][:-1]
            else:
                break
    return skip_rows


def map_seq_to_ncbi_tax_level_and_id(abundance_table, my_taxonomy_column, seqId, taxonomy_scheme, tax_delim, get_chiildren):
    """
    Parse user's OTU table and the Silva database to get to add 2 extra columns in the OTU table:
    1. the lowest taxonomic level of the taxonomy assigned in an OTU, for which an NCBI Taxonomy id exists (e.g., "genus")
    2. the corresponding NCBI Taxonomy Id (e.g., "343")

    Returns:
    a. splitted_taxonomies: a pandas df with the sequence id, the ncbi tax id of species, genus, family level when available
                        and the lowest taxonomic level for which there is an ncbi tax id
                        also a list (as column in th df) with gtdb genomes if available
                        species_ncbi_id genus_ncbi_id family_ncbi_id  tax_level  ncbi_tax_id,
    b. repr_genomes_present: a list with the GTDB genomes of interest
    """
    if taxonomy_scheme == "GTDB":
        taxon_to_ncbiId = os.path.join(MAPPINGS, "gtdbSpecies2ncbiId2accession.tsv")

        gtdb_accession_ids = pd.read_csv(taxon_to_ncbiId, sep="\t")
        gtdb_accession_ids.columns = ["Species", "ncbi_tax_id", "gtdb_gen_repr"]
        gtdb_accession_ids["Species"].str.strip()

        species_ncbi_id = pd.read_csv(SPECIES_NCBI_IDS, sep="\t")
        species_ncbi_id.columns = ['Species', 'ncbi_tax_id']
        species_ncbi_id['Species'].str.strip()

    elif taxonomy_scheme == "Silva":
        # [TODO] add the
        print("remember me!")

    else:
        taxon_to_ncbiId = os.path.join(MAPPINGS, "species2ncbiId2accession.tsv")
        taxon_to_ncbiId = pd.read_csv(taxon_to_ncbiId, sep="\t", names=["name", "ncbi", "gtdb"])

    # Split the taxonomy column and split it based on semicolumn!
    taxonomies = abundance_table.filter([seqId, my_taxonomy_column, "microbetag_id"])
    splitted_taxonomies = taxonomies[my_taxonomy_column].str.split(tax_delim, expand=True)
    splitted_taxonomies.columns = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    splitted_taxonomies[seqId] = taxonomies[seqId]
    splitted_taxonomies["microbetag_id"] = taxonomies["microbetag_id"]

    # Check if "s__" taxonomy-like and whether species in 1 or 2 steps
    s_pattern = "s__"
    if splitted_taxonomies['Species'].str.contains(s_pattern).all():
        pattern_for_complete_name = r"__[\s_]"
        if splitted_taxonomies['Species'].str.contains(pattern_for_complete_name).all():
            splitted_taxonomies["extendedSpecies"] = splitted_taxonomies['Species'].apply(process_underscore_taxonomy)
        else:
            # We need to use the Genus column too
            genus = splitted_taxonomies['Genus'].apply(process_underscore_taxonomy)
            species = splitted_taxonomies['Species'].apply(process_underscore_taxonomy)
            splitted_taxonomies["extendedSpecies"] = pd.DataFrame({'extendedSpecies': np.where(species.isna(), None, genus + ' ' + species)})
    else:
        splitted_taxonomies["extendedSpecies"] = splitted_taxonomies['Species']

    unique_species = splitted_taxonomies['extendedSpecies'].unique()
    unique_species = [item for item in unique_species if item is not None]

    # Run a pool to get NCBI ids at SPECIES level
    chunk_size = round(len(unique_species) / 2)
    chunks = [unique_species[i: i + chunk_size] for i in range(0, len(unique_species), chunk_size)]

    pool = multiprocessing.Pool(2)
    data = [(chunk, taxon_to_ncbiId) for chunk in chunks]
    ncbi_ids_species_level = pool.map(calculate_fuzzy_similarity_chunk, data)

    # Fix ncbi tax ids for species level
    ncbi_ids_species_name_matched_df = pd.DataFrame.from_dict(ncbi_ids_species_level[0], orient="index")
    ncbi_ids_species_name_df = ncbi_ids_species_name_matched_df["ncbi"]
    ncbi_ids_species_name = dict(ncbi_ids_species_name_df)

    splitted_taxonomies["species_ncbi_id"] = splitted_taxonomies["extendedSpecies"].map(ncbi_ids_species_name)

    # Get NCBI ids at GENUS level
    genera_ncbi_id = pd.read_csv(GENERA_NCBI_IDS, sep="\t")
    genera_ncbi_id.columns = ['Genus', 'ncbi_tax_id']
    genera_ncbi_id['Genus'].str.strip()

    # genera_to_get_ids = splitted_taxonomies.loc[splitted_taxonomies['species_ncbi_id'].isna(), 'Genus']
    genera = splitted_taxonomies['Genus'].apply(process_underscore_taxonomy)
    splitted_taxonomies["extendedGenus"] = genera
    genera.name = "Genus"
    genera = genera.to_frame()
    genera = genera.merge(genera_ncbi_id, on='Genus', how='inner').drop_duplicates()
    splitted_taxonomies = pd.merge(splitted_taxonomies, genera, left_on='extendedGenus', right_on='Genus', how='left')
    splitted_taxonomies = splitted_taxonomies.drop('Genus_y', axis=1)
    splitted_taxonomies = splitted_taxonomies.rename(columns={'ncbi_tax_id': 'genus_ncbi_id'})
    splitted_taxonomies['genus_ncbi_id'] = np.where(splitted_taxonomies['species_ncbi_id'].notna(), np.nan, splitted_taxonomies['genus_ncbi_id'])

    # Get NCBI ids at FAMILY level
    family_ncbi_id = pd.read_csv(FAMILIES_NCBI_IDS, sep="\t")
    family_ncbi_id.columns = ['Family', 'ncbi_tax_id']
    family_ncbi_id['Family'].str.strip()
    families = splitted_taxonomies['Family'].apply(process_underscore_taxonomy)
    splitted_taxonomies["extendedFamily"] = families
    families.name = "Family"
    families = families.to_frame()
    families = families.merge(family_ncbi_id, on='Family', how='inner').drop_duplicates()
    splitted_taxonomies = pd.merge(splitted_taxonomies, families, left_on='extendedFamily', right_on='Family', how='left')
    splitted_taxonomies = splitted_taxonomies.drop('Family_y', axis=1)
    splitted_taxonomies = splitted_taxonomies.rename(columns={'ncbi_tax_id': 'family_ncbi_id'})

    # Remove proxy taxon columns
    # splitted_taxonomies = splitted_taxonomies.drop('extendedSpecies', axis=1)
    splitted_taxonomies = splitted_taxonomies.drop('extendedGenus', axis=1)
    splitted_taxonomies = splitted_taxonomies.drop('extendedFamily', axis=1)

    # Get GTDB genomes for taxa available
    unique_species_present_ncbi_ids = list(splitted_taxonomies["species_ncbi_id"].astype(pd.Int64Dtype()).astype(str).unique())
    unique_species_present_ncbi_ids = [x for x in unique_species_present_ncbi_ids if x != "<NA>"]

    # Dictionary with ncbi ids of species level nodes and their gtdb genomes
    species_ncbi_ids_to_gtdb_genomes = {}
    for ncbi_id in unique_species_present_ncbi_ids:
        genomes = get_genomes_for_ncbi_tax_id(ncbi_id)
        gc_genomes = [genome for genome in list(genomes.values())[0] if genome.startswith("GCA_") or genome.startswith("GCF_")]
        if len(gc_genomes) > 0:
            species_ncbi_ids_to_gtdb_genomes[ncbi_id] = gc_genomes

    # A list with the gtdb present on the dataset
    repr_genomes_present = [item for sublist in list(species_ncbi_ids_to_gtdb_genomes.values()) for item in sublist]

    # Fix final df: add GTDB genomes to those with one
    splitted_taxonomies['gtdb_gen_repr'] = splitted_taxonomies['species_ncbi_id'].map(species_ncbi_ids_to_gtdb_genomes)

    # Now map those GTDB ids to their corresponding entries in the splitted_taxonomy df
    splitted_taxonomies['ncbi_tax_id'] = splitted_taxonomies.apply(assign_tax_id_for_node_level, axis=1).astype(str)
    splitted_taxonomies['ncbi_tax_id'] = splitted_taxonomies['ncbi_tax_id'].fillna('').astype(str).str.rstrip('.0').replace('', np.nan).astype(float).astype('Int64')

    # Keep what is the taxonomic level of the node
    splitted_taxonomies['ncbi_tax_level'] = splitted_taxonomies.apply(determine_tax_level, axis=1)

    # Remove any white spaces from the dataframe's columns
    splitted_taxonomies = splitted_taxonomies.map(lambda x: x.strip() if isinstance(x, str) else x)

    if get_chiildren:
        ncbi_parent_to_children = {}
        ncbi_nodes_dict = get_ncbi_nodes_dict()

        species_df = splitted_taxonomies[splitted_taxonomies['ncbi_tax_level'] == "species"]

        for potent_parent in species_df["ncbi_tax_id"]:
            potent_parent = str(int(potent_parent))
            if potent_parent in ncbi_nodes_dict:
                ncbi_parent_to_children[potent_parent] = ncbi_nodes_dict[potent_parent]

        rows = [(key, value) for key, values in ncbi_parent_to_children.items() for value in values]

        # Create a DataFrame from the list of tuples
        children_df = pd.DataFrame(rows, columns=['parent_ncbi_tax_id', 'child_ncbi_tax_id'])
        children_ncbi_ids_to_gtdb_genomes = {}
        for ncbi_id in children_df["child_ncbi_tax_id"]:
            genomes = get_genomes_for_ncbi_tax_id(ncbi_id)
            gc_genomes = [genome for genome in list(genomes.values())[0] if genome.startswith("GCA_") or genome.startswith("GCF_")]
            if len(gc_genomes) > 0:
                children_ncbi_ids_to_gtdb_genomes[ncbi_id] = gc_genomes

        children_df['gtdb_gen_repr'] = children_df['child_ncbi_tax_id'].map(children_ncbi_ids_to_gtdb_genomes)

    splitted_taxonomies = splitted_taxonomies.explode("gtdb_gen_repr")

    children_df = children_df.explode("gtdb_gen_repr")

    return splitted_taxonomies, repr_genomes_present, children_df


def calculate_fuzzy_similarity_chunk(args_set):
    """
    Gets a list of taxa names and checks what is their best hit against a ref taxonomy.
    """
    buzz_taxa = [
        "uncultured bacterium", "uncultured organism", "uncultured beta proteobacterium",
        "uncultured soil bacterium", "uncultured rumen bacterium", "uncultured gamma proteobacterium",
        "D_8__uncultured rumen protozoa", "uncultured delta proteobacterium"
    ]
    chunk, taxon_to_ncbiId_df = args_set[0], args_set[1]
    start = time.time()
    result = []
    _shared_dict = {}
    for taxon in chunk:
        taxon = taxon.strip()
        _shared_dict[taxon] = {}
        # check if it is not a species name at all
        tokens = re.split(r'[ _]', taxon)
        if len(tokens) == 1:
            print("entry '", taxon, "' is single word, thus not a species name")
            _shared_dict[taxon]["ncbi"] = None
        elif taxon in buzz_taxa:
            print("entry '", taxon, "' is not a species name")
            _shared_dict[taxon]["ncbi"] = None
        # check whether name exists as is in the ncbi db
        else:
            is_in_column = taxon_to_ncbiId_df['name'].isin([taxon])
            if is_in_column.any():
                print("!!! species name", taxon, " exactly as in NCBI Taxonomy")
                _shared_dict[taxon]["ncbi"] = taxon_to_ncbiId_df[is_in_column]["ncbi"].values[0]
                _shared_dict[taxon]["matched_name"] = taxon_to_ncbiId_df[is_in_column]["name"].values[0]
            else:
                print(".... ***looking for closest entry in NCBI Taxonomy for taxon '", taxon, "'")
                best_similarity = 0
                best_match = ""
                ncbi = ""
                for row in taxon_to_ncbiId_df.iterrows():
                    accurate_similarity = fuzz.ratio(row[1][0], taxon)
                    if accurate_similarity > best_similarity:
                        best_similarity = accurate_similarity
                        best_match = row[1][0]
                        ncbi = row[1][1]
                result.append((best_match, best_similarity))
                if best_similarity > 80:
                    _shared_dict[taxon]["matched_name"] = best_match
                    _shared_dict[taxon]["ncbi"] = ncbi
                else:
                    print(taxon, best_match, best_similarity)
                    _shared_dict[taxon]["ncbi"] = None
    end = time.time()
    print("\n\nTotal chunk time: ", str(end - start))
    return _shared_dict


def seqId_faprotax_functions_assignment(path_to_subtables):
    """
    Parse the sub tables of the faprotax analysis
    to assign the biological processes related to each sequence id
    """
    seqId_faprotax_assignments = {}

    for process_name in os.listdir(path_to_subtables):

        f = os.path.join(path_to_subtables, process_name)
        table_file = open(f, "r")
        table_file = table_file.readlines()

        for line in table_file[2:]:

            seqId = line.split("\t")[1]

            if seqId not in seqId_faprotax_assignments:

                seqId_faprotax_assignments[seqId] = [process_name]

            else:

                seqId_faprotax_assignments[seqId].append(process_name)

    return seqId_faprotax_assignments


def build_a_base_node(taxon, edge, taxonomy, side):
    """
    Builds a node for the base network
    """
    node = {}
    node["data"] = {}
    node["data"]["id"] = taxon
    node["data"]["selected"] = False
    node["data"]["taxonomy"] = taxonomy
    node["data"]["degree_layout"] = 1
    node["data"]["name"] = taxonomy.split(";")[-1]
    node["data"]["NCBI-Tax-Id"] = edge["ncbi_tax_id_node_" + side]  # ) if not pd.isna(edge["ncbi_tax_id_node_" + side]) else edge["ncbi_tax_id_node_" + side]
    node["data"]["GTDB-representative"] = edge["gtdb_gen_repr_node_" + side]
    node["data"]["taxonomy-level"] = edge["ncbi_tax_level_node_" + side]

    return node


def get_ncbi_nodes_dict():
    """
    Returns all children NCBI Taxonomy IDs of a parent.
    """
    nodes_dict = {}
    with open(NCBI_NODES_TMP, 'r') as nodes:
        for line in nodes:
            fields = line.strip().split('\t|\t')
            tax_id, parent_id = fields[0], fields[1]
            if parent_id not in nodes_dict:
                nodes_dict[parent_id] = [tax_id]
            else:
                nodes_dict[parent_id].append(tax_id)
    return nodes_dict


def build_base_graph(edgelist_as_a_list_of_dicts, microb_id_taxonomy):
    """
    Get a list of dictionaries where each dictionary is an edge and returns
    the basenetwork in a JSON format.
    """
    base_network = {}
    base_network["elements"] = {}

    nodes = []
    edges = []

    processed_nodes = set()

    counter = 1

    for edge in edgelist_as_a_list_of_dicts:

        taxon_a = edge["node_a"]  # microbetag_id
        taxonomy_a = microb_id_taxonomy.loc[microb_id_taxonomy['microbetag_id'] == taxon_a, 'taxonomy'].item()
        if taxon_a not in processed_nodes:
            processed_nodes.add(taxon_a)
            node_a = build_a_base_node(taxon_a, edge, taxonomy_a, "a")
            nodes.append(node_a)

        taxon_b = edge["node_b"]  # microbetag_id
        taxonomy_b = microb_id_taxonomy.loc[microb_id_taxonomy['microbetag_id'] == taxon_b, 'taxonomy'].item()
        if taxon_b not in processed_nodes:
            processed_nodes.add(taxon_b)
            node_b = build_a_base_node(taxon_b, edge, taxonomy_b, "b")
            nodes.append(node_b)

        new_edge = {}
        new_edge["data"] = {}
        new_edge["data"]["id"] = str(counter)
        new_edge["data"]["source"] = taxon_a
        new_edge["data"]["source-ncbi-tax-id"] = edge["ncbi_tax_id_node_a"]
        new_edge["data"]["target"] = taxon_b
        new_edge["data"]["target-ncbi-tax-id"] = edge["ncbi_tax_id_node_b"]
        new_edge["data"]["selected"] = False
        new_edge["data"]["shared_name"] = taxonomy_a.split(";")[-1] + "-" + taxonomy_b.split(";")[-1]
        new_edge["data"]["SUID"] = str(counter)
        new_edge["data"]["name"] = "co-occurrence"
        new_edge["data"]["weight"] = float(edge["score"])
        new_edge["selected"] = False

        edges.append(new_edge)
        counter += 1

    # Ensure .cyjs format
    base_network["elements"]["nodes"] = nodes
    base_network["elements"]["edges"] = edges

    base_network["data"] = {}
    base_network["data"]["title"] = "microbetag annotated microbial co-occurrence network"
    base_network["data"]["tags"] = ["v1.0"]

    return base_network


def is_tab_separated(my_abundance_table, my_taxonomy_column):
    """
    Read the OTU table and make sure it is tab separated and not empty
    Takes as input a .tsv file and returns a pandas dataframe.
    """
    number_of_commented_lines = count_comment_lines(
        my_abundance_table, my_taxonomy_column)
    try:
        abundance_table = pd.read_csv(
            my_abundance_table,
            sep=ABUNDANCE_TABLE_DELIM,
            skiprows=number_of_commented_lines)
    except BaseException:
        logging.error(
            "The OTU table provided is not a tab separated file. Please convert your OTU table to .tsv or .csv format.")

    if abundance_table.shape[1] < 2:
        logging.error("The OTU table you provide has no records.")

    return abundance_table


def ensure_flashweave_format(my_abundance_table, my_taxonomy_column, seqId, outdir):
    """
    Build an OTU table that will be in a FlashWeave-based format.
    """
    flashweave_table = my_abundance_table.drop(my_taxonomy_column, axis=1)
    float_col = flashweave_table.select_dtypes(include=['float64'])

    for col in float_col.columns.values:
        flashweave_table[col] = flashweave_table[col].astype('int64')

    flashweave_table[seqId] = flashweave_table[seqId].astype(str)
    my_abundance_table['microbetag_id'] = flashweave_table[seqId]

    file_to_save = os.path.join(
        outdir,
        "abundance_table_flashweave_format.tsv")
    flashweave_table.to_csv(file_to_save, sep='\t', index=False)

    return my_abundance_table


def edge_list_of_ncbi_ids(edgelist, abundance_table_with_ncbi_ids):
    """
    Read an edge list and build a dataframe with the corresponding ncbi ids for each pair
    if and only if, both OTUs have been mapped to a NCBI tax id
    e.g.
                ncbi_tax_id_a  ncbi_tax_level_a  ncbi_tax_id_b   ncbi_tax_level_b
          0        838              genus           171552          family
          1       186803           family           186807          family
    """
    pd_edgelist = pd.read_csv(edgelist, sep="\t", skiprows=2, header=None)
    pd_edgelist.columns = ["node_a", "node_b", "score"]

    pd_edgelist["joint"] = pd_edgelist['node_a'].astype(str) + ":" + pd_edgelist["node_b"]

    # The pd.explode() function transforms the arrays into separate rows
    abundance_table_with_ncbi_ids_exploded = abundance_table_with_ncbi_ids.explode("gtdb_gen_repr")

    associated_pairs_node_a = pd.merge(
        pd_edgelist[["node_a", "joint"]],
        abundance_table_with_ncbi_ids_exploded[["ncbi_tax_id", "gtdb_gen_repr", "ncbi_tax_level", "microbetag_id"]],
        left_on='node_a', right_on='microbetag_id', how="inner"
    ).drop(["microbetag_id"], axis=1)

    associated_pairs_node_a.rename(columns={
        "ncbi_tax_level": "ncbi_tax_level_node_a",
        "gtdb_gen_repr": "gtdb_gen_repr_node_a",
        "ncbi_tax_id": "ncbi_tax_id_node_a"
    }, inplace=True)

    associated_pairs_node_b = pd.merge(
        pd_edgelist[["node_b", "joint", "score"]],
        abundance_table_with_ncbi_ids_exploded[["ncbi_tax_id", "gtdb_gen_repr", "ncbi_tax_level", "microbetag_id"]],
        left_on='node_b', right_on='microbetag_id', how="inner"
    ).drop(["microbetag_id"], axis=1)

    associated_pairs_node_b.rename(columns={
        "ncbi_tax_level": "ncbi_tax_level_node_b",
        "gtdb_gen_repr": "gtdb_gen_repr_node_b",
        "ncbi_tax_id": "ncbi_tax_id_node_b"
    }, inplace=True)

    associated_pairs = pd.merge(associated_pairs_node_a, associated_pairs_node_b, on="joint").drop(["joint"], axis=1)

    # associated_pairs['ncbi_tax_id_node_a'] = associated_pairs['ncbi_tax_id_node_a'].astype(pd.Int64Dtype())
    # associated_pairs['ncbi_tax_id_node_b'] = associated_pairs['ncbi_tax_id_node_b'].astype(pd.Int64Dtype())

    return associated_pairs


def export_phen_traits_to_file(column_names, rows, filename):
    with open(filename, 'w') as file:
        file.write('\t'.join(column_names) + '\n')
        for row in rows:
            file.write('\t'.join(row) + '\n')
    return True


def process_underscore_taxonomy(entry):
    """
    Convert "s__" like taxa to ncbi like ones
    """
    parts = entry.split("__")
    if len(parts) == 1:
        return entry
    if len(parts[1]) > 0:
        return parts[1].replace("_", " ")
    else:
        return np.nan


def tuple_to_str(t):
    """
    Joins members of a tuple to a comma-separated string
    """
    return ','.join(map(str, t))


def convert_tuples_to_strings(obj):
    """
    Recursive function to convert all keys of tuple structure in any level of a nested dictionary to strings.
    """
    if isinstance(obj, tuple):
        return str(obj)
    if isinstance(obj, dict):
        return {convert_tuples_to_strings(key): convert_tuples_to_strings(value) for key, value in obj.items()}
    return obj


def determine_tax_level(row):
    if isinstance(row['gtdb_gen_repr'], list):
        return 'mspecies'
    elif not pd.isna(row['species_ncbi_id']):
        return 'species'
    elif not pd.isna(row['genus_ncbi_id']):
        return 'genus'
    elif not pd.isna(row['family_ncbi_id']):
        return 'family'
    else:
        return 'not known'


def assign_tax_id_for_node_level(row):
    if isinstance(row['gtdb_gen_repr'], list):
        return row['species_ncbi_id']
    elif not pd.isna(row['species_ncbi_id']):
        return row['species_ncbi_id']
    elif not pd.isna(row['genus_ncbi_id']):
        return row['genus_ncbi_id']
    elif not pd.isna(row['family_ncbi_id']):
        return row['family_ncbi_id']
    else:
        return np.nan


def convert_to_float(data):
    """
    Convert float numbers recursively in a nested dict.
    """
    if isinstance(data, str):
        try:
            return float(data)
        except ValueError:
            return data
    elif isinstance(data, dict):
        return {key: convert_to_float(value) for key, value in data.items()}
    elif isinstance(data, list):
        return [convert_to_float(item) for item in data]
    else:
        return data


def remove_query(dictionary):
    if isinstance(dictionary, dict):
        # Create a copy of the dictionary to avoid modifying it while iterating
        dictionary_copy = dictionary.copy()

        # Iterate over the items in the dictionary
        for key, value in dictionary.items():
            if key == "query":
                # Remove the "query" key
                del dictionary_copy[key]
            elif isinstance(value, dict):
                # Recursively call the function for nested dictionaries
                dictionary_copy[key] = remove_query(value)
        return dictionary_copy
    else:
        return dictionary


def read_cyjson(filename, direction=False):
    """Function based on the corresponding of the manta library: https://github.com/ramellose/manta/blob/master/manta/cyjson.py
    Small utility function for reading Cytoscape json files
    generated with CoNet.
    In our case, it also gets the layout and adds it as part of the node data.
    Parameters
    ----------
    :param filename: Filepath to location of cyjs file.
    :param direction: If true, graph is imported as a NetworkX DiGraph
    :return: NetworkX graph.
    """
    with open(filename) as f:
        data = json.load(f)
    name = 'name'
    ident = 'id'
    if len(set([ident, name])) < 2:
        raise nx.NetworkXError('Attribute names are not unique.')
    if direction:
        graph = nx.DiGraph()
    else:
        graph = nx.Graph()
    graph.graph = dict(data.get('data'))
    i = 0
    for d in data["elements"]["nodes"]:
        # only modification: 'value' key is not included in CoNet output
        # now graph only needs ID and name values
        node_data = d["data"].copy()
        position = d["position"]
        node_data["position"] = position
        try:
            node = d["data"].get(ident)
        except KeyError:
            # if no index is found, one is generated
            node = i
            i += 1
        if d["data"].get(name):
            node_data[name] = d["data"].get(name)

        graph.add_node(node)
        graph.nodes[node].update(node_data)
    for d in data["elements"]["edges"]:
        edge_data = d["data"].copy()
        sour = d["data"].pop("source")
        targ = d["data"].pop("target")
        graph.add_edge(sour, targ)
        graph.edges[sour, targ].update(edge_data)
    if 'interactionType' in graph.edges[list(graph.edges)[0]]:
        # this indicates this is a CoNet import
        for edge in graph.edges:
            if graph.edges[edge]['interactionType'] == 'copresence':
                graph.edges[edge]['weight'] = 1
            elif graph.edges[edge]['interactionType'] == 'mutualExclusion':
                graph.edges[edge]['weight'] = -1
            else:
                graph.edges[edge]['weight'] = 0
    return graph
