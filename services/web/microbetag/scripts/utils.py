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
import ast


def export_species_level_associations(edgelist_as_a_list_of_dicts, seqID_taxid_level_repr_genome, children_df=None):
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
    pairs_of_ncbi_id_of_interest = set()
    pairs_of_seqId_of_interest = set()
    list_of_non_usefules_ids = ["77133", "91750"]
    for pair in edgelist_as_a_list_of_dicts:
        # print(pair)
        # print(pair["node_b"])
        # print(
        #     seqID_taxid_level_repr_genome[seqID_taxid_level_repr_genome["seqId"] == pair["node_b"]]
        # )
        # taxon_a and taxon_b variables are int type; in case there is a <NA> value then it is a pandas missing NAType
        taxon_a = seqID_taxid_level_repr_genome[seqID_taxid_level_repr_genome["seqId"] == pair["node_a"]]["ncbi_tax_id"].values.item()
        taxon_b = seqID_taxid_level_repr_genome[seqID_taxid_level_repr_genome["seqId"] == pair["node_b"]]["ncbi_tax_id"].values.item()
        taxon_a_level = seqID_taxid_level_repr_genome[seqID_taxid_level_repr_genome["seqId"] == pair["node_a"]]["ncbi_tax_level"].values.item()
        taxon_b_level = seqID_taxid_level_repr_genome[seqID_taxid_level_repr_genome["seqId"] == pair["node_b"]]["ncbi_tax_level"].values.item()
        if str(taxon_a) in list_of_non_usefules_ids or str(taxon_b) in list_of_non_usefules_ids:
            continue
        # We keep as a pair both a->b and b->a association
        if taxon_a_level == taxon_b_level == "mspecies":
            set_of_ncbiids_of_interest.add(taxon_a)
            set_of_ncbiids_of_interest.add(taxon_b)
            pairs_of_ncbi_id_of_interest.add((taxon_a, taxon_b))
            pairs_of_ncbi_id_of_interest.add((taxon_b, taxon_a))
            pairs_of_seqId_of_interest.add((pair["node_a"], pair["node_b"]))
    # Start building dics
    related_genomes = {}
    children_genomes = {}
    for ncbiId in set_of_ncbiids_of_interest:
        ncbiId_genomes = seqID_taxid_level_repr_genome[
            seqID_taxid_level_repr_genome["ncbi_tax_id"] == ncbiId
            ]["gtdb_gen_repr"].to_list()[0]
        related_genomes[ncbiId] = ncbiId_genomes

        if children_df is not None:
            children_ncbiId_genomes = children_df[(children_df["parent_ncbi_tax_id"] == ncbiId) & 
                                                  (children_df["gtdb_gen_repr"].notna())
                                                  ][["child_ncbi_tax_id", "gtdb_gen_repr"]].to_dict(orient="records")
            for case in children_ncbiId_genomes:
                c_genomes = case["gtdb_gen_repr"]
                children_genomes[case["child_ncbi_tax_id"]] = c_genomes
                related_genomes[ncbiId].append(c_genomes)

    # for k, v in related_genomes.items():
    #     related_genomes[k] = flatten_list(v)

    return pairs_of_ncbi_id_of_interest, pairs_of_seqId_of_interest, related_genomes, children_genomes


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


def aggregate_genomes(genomes):
    """
    Aggregate rows of a df with the same seqId and different gtdb_repr_gen so 
    you have a single row with a list of the genomes found
    """
    unique_genomes = set()
    for genome_list in genomes:
        unique_genomes.update(genome_list)
    return list(unique_genomes)


def replace_empty_list(value):
    """
    Replace a df's cell that is an empty list with np.nan 
    """
    return np.nan if isinstance(value, list) and not value else value


def map_seq_to_ncbi_tax_level_and_id(abundance_table, my_taxonomy_column, seqId, taxonomy_scheme, tax_delim, get_chiildren):
    """
    Parse user's OTU table and the Silva database to get to add 2 extra columns in the OTU table:
    1. the lowest taxonomic level of the taxonomy assigned in an OTU, for which an NCBI Taxonomy id exists (e.g., "genus")
    2. the corresponding NCBI Taxonomy Id (e.g., "343")

    Returns:
    splitted_taxonomies: a pandas df with the sequence id, the ncbi tax id of species, genus, family level when available
                        and the lowest taxonomic level for which there is an ncbi tax id
                        also a list (as column in th df) with gtdb genomes if available
                        species_ncbi_id genus_ncbi_id family_ncbi_id  tax_level  ncbi_tax_id,
    repr_genomes_present: a list with the GTDB genomes of interest
    children_df: (optional)
    """
    # Split the taxonomy column and split it based on semicolumn!
    taxonomies = abundance_table.filter([seqId, my_taxonomy_column, "microbetag_id"])
    splitted_taxonomies = taxonomies[my_taxonomy_column].str.split(tax_delim, expand=True)
    if splitted_taxonomies[0].str.contains('Root').all():
        splitted_taxonomies = splitted_taxonomies.drop(splitted_taxonomies.columns[0], axis=1)

    if len(list(splitted_taxonomies.columns)) != 7:
        raise ValueError(f"Error: The taxonomy scheme provided is not a 7-level one.")

    splitted_taxonomies.columns = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    splitted_taxonomies[seqId] = taxonomies[seqId]
    splitted_taxonomies["microbetag_id"] = taxonomies["microbetag_id"]

    # Check if "s__" taxonomy-like and whether species in 1 or 2 steps
    s_pattern = "s__"
    if splitted_taxonomies['Species'].str.contains(s_pattern).all():
        underscore = True
        pattern_for_complete_name = r"__[\s_]"
        if splitted_taxonomies['Species'].str.contains(pattern_for_complete_name).all():
            splitted_taxonomies["extendedSpecies"] = splitted_taxonomies['Species'].apply(process_underscore_taxonomy)
        else:
            # We might need to use the Genus column too, check if genus is part of species
            genus = splitted_taxonomies['Genus'].apply(process_underscore_taxonomy)
            species = splitted_taxonomies['Species'].apply(process_underscore_taxonomy)
            genus_list = genus.tolist() ; genus_list = [str(c) for c in genus_list]
            species_list = species.tolist() ; species_list = [str(c) for c in species_list]
            results = [True if genus in species else False for genus, species in zip(genus_list, species_list) if species != 'nan' ]
            # Check if genus is included in species for non-empty species
            if results.count(True) > 0.8*len(results):
                splitted_taxonomies["extendedSpecies"] = pd.DataFrame({'extendedSpecies': np.where(species.isna(), None, species)})
            else: 
                splitted_taxonomies["extendedSpecies"] = pd.DataFrame({'extendedSpecies': np.where(species.isna(), None, genus + ' ' + species)})
    else:
        splitted_taxonomies["extendedSpecies"] = splitted_taxonomies['Species']
        underscore = False

    """
    Get NCBI ids at SPECIES level - Use proper ref file
    NOTE: If the taxon_to_ncbiId file has more than one ncbi ids or genomes for a taxonomy is going to be an issue 
    
    Unique taxonomies: a seqId will hit exactly one genome
    ------------------------------------------------------
    gc_accession_16s_gtdb_ncbid.tsv
    gtdbSpecies2ncbiId2accession.tsv


    fuzzywazzy output
    Silva (gtdb_silvaSpecies2ncbi2accession) > 
    
    """
    if taxonomy_scheme == "GTDB":
        taxon_to_ncbiId = os.path.join(MAPPINGS, "gtdbSpecies2ncbiId2accession.tsv")

    elif taxonomy_scheme == "Silva":
        taxon_to_ncbiId = os.path.join(MAPPINGS, "gtdb_silvaSpecies2ncbi2accession.tsv")

    elif taxonomy_scheme == "microbetag_prep":
        taxon_to_ncbiId = os.path.join(MAPPINGS, "gc_accession_16s_gtdb_ncbid.tsv")

    else:
        taxon_to_ncbiId = os.path.join(MAPPINGS, "species2ncbiId2accession.tsv")


    # Get species NCBI ids
    if taxonomy_scheme == "GTDB" or taxonomy_scheme == "Silva" or taxonomy_scheme == "microbetag_prep":
        gtdb_accession_ids = pd.read_csv(taxon_to_ncbiId, sep="\t")
        gtdb_accession_ids.columns = ["refSpecies", "species_ncbi_id", "gtdb_gen_repr"]
        gtdb_accession_ids["refSpecies"].str.strip()

        splitted_taxonomies =  pd.merge(splitted_taxonomies, gtdb_accession_ids, left_on='extendedSpecies', right_on='refSpecies', how='left')
        repr_genomes_present = [ c for c in splitted_taxonomies["gtdb_gen_repr"].to_list() if isinstance(c, str) ]
        
        if taxonomy_scheme != "Silva":
            splitted_taxonomies['gtdb_gen_repr'] = splitted_taxonomies['gtdb_gen_repr'].apply(lambda x: [x] if pd.notna(x) else np.nan)

        else:
            # Remove cases where for a single seqId there are more than 2 ncbi taxids in case both hit the same genome
            splitted_taxonomies['gtdb_gen_repr'] = splitted_taxonomies['gtdb_gen_repr'].apply(lambda x: tuple([x]) if isinstance(x, str) else tuple())

            # Group by 'seqId' and aggregate 'gtdb_gen_repr' with custom function
            splitted_taxonomies = splitted_taxonomies.groupby('seqId').agg({
                **{col: 'first' for col in splitted_taxonomies.columns if col not in ['seqId', 'gtdb_gen_repr']},
                'gtdb_gen_repr': aggregate_genomes
            }).reset_index()

            splitted_taxonomies = splitted_taxonomies.map(replace_empty_list)

        gtdb_genomes_on = True


    else:
        taxon_to_ncbiId_df = pd.read_csv(taxon_to_ncbiId, sep="\t", names=["name", "ncbi", "gtdb"])
        unique_species = splitted_taxonomies['extendedSpecies'].unique()
        unique_species = [item for item in unique_species if item is not None]

        # Run a pool to get NCBI ids at SPECIES level
        chunk_size = round(len(unique_species) / 2)
        chunks = [unique_species[i: i + chunk_size] for i in range(0, len(unique_species), chunk_size)]
        pool = multiprocessing.Pool(2)
        data = [(chunk, taxon_to_ncbiId_df) for chunk in chunks]
        ncbi_ids_species_level = pool.map(calculate_fuzzy_similarity_chunk, data)

        # Fix ncbi tax ids for species level
        ncbi_ids_species_name_matched_df = pd.DataFrame.from_dict(ncbi_ids_species_level[0], orient="index")
        ncbi_ids_species_name_df = ncbi_ids_species_name_matched_df["ncbi"]
        ncbi_ids_species_name = dict(ncbi_ids_species_name_df)

        splitted_taxonomies["species_ncbi_id"] = splitted_taxonomies["extendedSpecies"].map(ncbi_ids_species_name)

        gtdb_genomes_on = False

    """
    Get NCBI ids at GENUS level
    """
    genera_ncbi_id = pd.read_csv(GENERA_NCBI_IDS, sep="\t")
    genera_ncbi_id.columns = ['Genus', 'ncbi_tax_id']
    genera_ncbi_id['Genus'].str.strip()

    if underscore:
        genera = splitted_taxonomies['Genus'].apply(process_underscore_taxonomy)
    else:
        genera = splitted_taxonomies['Genus']
    splitted_taxonomies["extendedGenus"] = genera
    genera.name = "Genus"
    genera = genera.to_frame()
    genera = genera.merge(genera_ncbi_id, on='Genus', how='inner').drop_duplicates()
    splitted_taxonomies = pd.merge(splitted_taxonomies, genera, left_on='extendedGenus', right_on='Genus', how='left')
    splitted_taxonomies = splitted_taxonomies.drop('Genus_y', axis=1)
    splitted_taxonomies = splitted_taxonomies.rename(columns={'ncbi_tax_id': 'genus_ncbi_id'})
    splitted_taxonomies['genus_ncbi_id'] = np.where(splitted_taxonomies['species_ncbi_id'].notna(), np.nan, splitted_taxonomies['genus_ncbi_id'])

    """
    Get NCBI ids at FAMILY level
    """
    family_ncbi_id = pd.read_csv(FAMILIES_NCBI_IDS, sep="\t")
    family_ncbi_id.columns = ['Family', 'ncbi_tax_id']
    family_ncbi_id['Family'].str.strip()
    if underscore:
        families = splitted_taxonomies['Family'].apply(process_underscore_taxonomy)
    else:
        families = splitted_taxonomies['Family']
    splitted_taxonomies["extendedFamily"] = families
    families.name = "Family"
    families = families.to_frame()
    families = families.merge(family_ncbi_id, on='Family', how='inner').drop_duplicates()
    splitted_taxonomies = pd.merge(splitted_taxonomies, families, left_on='extendedFamily', right_on='Family', how='left')
    splitted_taxonomies = splitted_taxonomies.drop('Family_y', axis=1)
    splitted_taxonomies = splitted_taxonomies.rename(columns={'ncbi_tax_id': 'family_ncbi_id'})

    # Remove proxy taxon columns
    splitted_taxonomies = splitted_taxonomies.drop('extendedGenus', axis=1)
    splitted_taxonomies = splitted_taxonomies.drop('extendedFamily', axis=1)

    """
    Get GTDB genomes for taxa available
    """
    if not gtdb_genomes_on:
        unique_species_present_ncbi_ids = [int(x) for x in list(splitted_taxonomies[splitted_taxonomies["species_ncbi_id"].notna()]["species_ncbi_id"])]

        # Dictionary with ncbi ids of species level nodes and their gtdb genomes
        species_ncbi_ids_to_gtdb_genomes = {}
        for ncbi_id in unique_species_present_ncbi_ids:
            genomes = get_genomes_for_ncbi_tax_id(ncbi_id)
            gc_genomes = [genome for genome in list(genomes.values())[0] if genome.startswith("GCA_") or genome.startswith("GCF_")]
            if len(gc_genomes) > 0:
                species_ncbi_ids_to_gtdb_genomes[float(ncbi_id)] = gc_genomes

        # A list with the gtdb present on the dataset
        repr_genomes_present = [item for sublist in list(species_ncbi_ids_to_gtdb_genomes.values()) for item in sublist]

        # Now map those GTDB ids to their corresponding entries in the splitted_taxonomy df
        splitted_taxonomies['gtdb_gen_repr'] = splitted_taxonomies['species_ncbi_id'].map(species_ncbi_ids_to_gtdb_genomes)
        # {float(k): v[0] for k, v in species_ncbi_ids_to_gtdb_genomes.items()}

    # Keep track of the taxonomic level of the ncbi tax id of the lowest level found
    splitted_taxonomies['ncbi_tax_id'] = splitted_taxonomies.apply(assign_tax_id_for_node_level, axis=1).astype(str)
    splitted_taxonomies['ncbi_tax_id'] = pd.to_numeric(splitted_taxonomies['ncbi_tax_id'], errors='coerce').astype('Int64').astype(str)
    # splitted_taxonomies['ncbi_tax_id'] = splitted_taxonomies['ncbi_tax_id'].fillna('').astype(str) #.str.rstrip('.0').replace('', np.nan).astype(float).astype('Int64')

    # Keep what is the taxonomic level of the node
    splitted_taxonomies['ncbi_tax_level'] = splitted_taxonomies.apply(determine_tax_level, axis=1)

    # Remove any white spaces from the dataframe's columns
    splitted_taxonomies = splitted_taxonomies.map(lambda x: x.strip() if isinstance(x, str) else x)

    # Beautify df
    splitted_taxonomies = splitted_taxonomies.rename(columns = {"Genus_x": "Genus", "Family_x": "Family"})
    desired_order = ["microbetag_id", seqId, "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "extendedSpecies",  
                     "family_ncbi_id", "genus_ncbi_id", "species_ncbi_id", "ncbi_tax_id", "ncbi_tax_level", "gtdb_gen_repr"]
    splitted_taxonomies = splitted_taxonomies[desired_order]

    # Get children NCBI Tax ids and their corresponding genomes
    if get_chiildren:
        ncbi_parent_to_children = {}
        ncbi_nodes_dict = get_ncbi_nodes_dict()
        species_df = splitted_taxonomies[(splitted_taxonomies['ncbi_tax_level'] == "species") | (splitted_taxonomies['ncbi_tax_level'] == "mspecies")]
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
        children_df['parent_ncbi_tax_id'] = children_df['parent_ncbi_tax_id'].astype(int)
        children_df['child_ncbi_tax_id'] = children_df['child_ncbi_tax_id'].astype(int)
        repr_genomes_present = repr_genomes_present + \
            [x for item in children_df['gtdb_gen_repr'].to_list() if not (isinstance(item, float) and np.isnan(item)) for x in item]
        return splitted_taxonomies, repr_genomes_present, children_df

    else:
        return splitted_taxonomies, repr_genomes_present


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
        # Check if it is not a species name at all
        tokens = re.split(r'[ _]', taxon)
        if len(tokens) == 1:
            _shared_dict[taxon]["ncbi"] = None
        elif taxon in buzz_taxa:
            _shared_dict[taxon]["ncbi"] = None
        # Check whether name exists as is in the ncbi db
        else:
            is_in_column = taxon_to_ncbiId_df['name'].isin([taxon])
            # Species name exactly as in NCBI Taxonomy
            if is_in_column.any():
                _shared_dict[taxon]["ncbi"] = taxon_to_ncbiId_df[is_in_column]["ncbi"].values[0]
                _shared_dict[taxon]["matched_name"] = taxon_to_ncbiId_df[is_in_column]["name"].values[0]
            # Run fuzzywuzzy to get closest
            else:
                print(">> Looking for closest entry in NCBI Taxonomy for taxon '", taxon, "'")
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
    for subtable_file in os.listdir(path_to_subtables):
        f = os.path.join(path_to_subtables, subtable_file)
        process_name = subtable_file.split(".")[0].replace("_", " ")
        table_file = open(f, "r")
        table_file = table_file.readlines()
        for line in table_file[2:]:
            seqId = line.split("\t")[1]
            if seqId not in seqId_faprotax_assignments:
                seqId_faprotax_assignments[seqId] = [process_name]
            else:
                seqId_faprotax_assignments[seqId].append(process_name)
    return seqId_faprotax_assignments


def build_a_base_node(taxon, map_seq, cfg):
    """
    Builds a node for the base network.
    [TODO] Remove not necessary entries.
    """
    case = map_seq[map_seq["seqId"] == taxon]
    node = {}
    node["data"] = {}
    node["data"]["id"] = taxon
    node["data"]["selected"] = False
    node["data"]["taxonomy"] = case["taxonomy"].item()
    node["data"]["degree_layout"] = 1
    node["data"]["name"] = case["taxonomy"].item().split(cfg["delimiter"])[-1]
    node["data"]["NCBI-Tax-Id"] = case["ncbi_tax_id"].item()
    node["data"]["GTDB-representative"] = case["gtdb_gen_repr"]
    node["data"]["taxonomy-level"] = case["ncbi_tax_level"].item()
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


def build_base_graph(edgelist_as_a_list_of_dicts, microb_id_taxonomy, cfg):
    """
    Runs if manta has been asked for from the user.
    manta gets a .cyjs input file.
    This function builds an non-annotated graph using only the scores and the taxonomies of the taxa of the network.
    It get a list of dictionaries where each dictionary is an edge and returns the basenetwork in a .cyjs format.
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
            node_a = build_a_base_node(taxon_a, microb_id_taxonomy, cfg)
            nodes.append(node_a)
        taxon_b = edge["node_b"]
        taxonomy_b = microb_id_taxonomy.loc[microb_id_taxonomy['microbetag_id'] == taxon_b, 'taxonomy'].item()
        if taxon_b not in processed_nodes:
            processed_nodes.add(taxon_b)
            node_b = build_a_base_node(taxon_b, microb_id_taxonomy, cfg)
            nodes.append(node_b)
        new_edge = {}
        new_edge["data"] = {}
        new_edge["data"]["id"] = str(counter)
        new_edge["data"]["source"] = taxon_a
        # new_edge["data"]["source-ncbi-tax-id"] = edge["ncbi_tax_id_node_a"]
        new_edge["data"]["source-ncbi-tax-id"] = microb_id_taxonomy[microb_id_taxonomy["seqId"] == taxon_a]["ncbi_tax_id"]
        new_edge["data"]["target"] = taxon_b
        # new_edge["data"]["target-ncbi-tax-id"] = edge["ncbi_tax_id_node_b"]
        new_edge["data"]["target-ncbi-tax-id"] = microb_id_taxonomy[microb_id_taxonomy["seqId"] == taxon_b]["ncbi_tax_id"]
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


def edge_list_of_ncbi_ids(edgelist, metadata_file=None):  # abundance_table_with_ncbi_ids
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

    if metadata_file:
        elements_to_exclude = pd.read_csv(metadata_file, sep="\t", header=None, index_col=0).index.to_list() 
        mask =pd_edgelist.isin(elements_to_exclude).any(axis=1)
        pd_edgelist = pd_edgelist[~mask]

    pd_edgelist["pair-of-taxa"] = pd_edgelist['node_a'].astype(str) + ":" + pd_edgelist["node_b"]

    # # The pd.explode() function transforms the arrays into separate rows
    # abundance_table_with_ncbi_ids_exploded = abundance_table_with_ncbi_ids.explode("gtdb_gen_repr")
    # associated_pairs_node_a = pd.merge(
    #     pd_edgelist[["node_a", "joint"]],
    #     abundance_table_with_ncbi_ids_exploded[["ncbi_tax_id", "gtdb_gen_repr", "ncbi_tax_level", "microbetag_id"]],
    #     left_on='node_a', right_on='microbetag_id', how="inner"
    # ).drop(["microbetag_id"], axis=1)
    # associated_pairs_node_a.rename(columns={
    #     "ncbi_tax_level": "ncbi_tax_level_node_a",
    #     "gtdb_gen_repr": "gtdb_gen_repr_node_a",
    #     "ncbi_tax_id": "ncbi_tax_id_node_a"
    # }, inplace=True)
    # associated_pairs_node_b = pd.merge(
    #     pd_edgelist[["node_b", "joint", "score"]],
    #     abundance_table_with_ncbi_ids_exploded[["ncbi_tax_id", "gtdb_gen_repr", "ncbi_tax_level", "microbetag_id"]],
    #     left_on='node_b', right_on='microbetag_id', how="inner"
    # ).drop(["microbetag_id"], axis=1)
    # associated_pairs_node_b.rename(columns={
    #     "ncbi_tax_level": "ncbi_tax_level_node_b",
    #     "gtdb_gen_repr": "gtdb_gen_repr_node_b",
    #     "ncbi_tax_id": "ncbi_tax_id_node_b"
    # }, inplace=True)
    # associated_pairs = pd.merge(associated_pairs_node_a, associated_pairs_node_b, on="joint").drop(["joint"], axis=1)
    return pd_edgelist


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


def convert_to_json_serializable(obj):
    """
    Recursively serializes entries of an object, i.e. a set is converted to a list, a list is split to its items
    and a dictionary keeps its key and their values get serialized
    """
    if isinstance(obj, (int, float, str, bool, type(None))):
        return obj
    elif isinstance(obj, set):
        return list(obj)
    elif isinstance(obj, list):
        return [convert_to_json_serializable(item) for item in obj]
    elif isinstance(obj, dict):
        return {key: convert_to_json_serializable(value) for key, value in obj.items()}
    else:
        try:
            return json.dumps(obj)
        except TypeError:
            return str(obj)


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
    """ Recursive function to remove the sql quert from a dictionary"""
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
    """
    Function based on the corresponding of the manta library: 
    https://github.com/ramellose/manta/blob/master/manta/cyjson.py
    
    Small utility function for reading Cytoscape json files generated with CoNet.
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
    return graph


def flatten_list(lista, flat_list=[]):
    """
    Recursive function taking as input a nested list and returning a flatten one.
    E.g. ['GCF_003252755.1', 'GCF_900638025.1', 'GCF_003252725.1', 'GCF_003253005.1', 'GCF_003252795.1', ['GCF_000210895.1'], ['GCF_000191405.1']]
    becomes ['GCF_003252755.1', 'GCF_900638025.1', 'GCF_003252725.1', 'GCF_003253005.1', 'GCF_003252795.1'].
    """
    for i in lista:
        if isinstance(i, list):
            flatten_list(i, flat_list)
        else:
            flat_list.append(i)
    return set(flat_list)


def build_cx_annotated_graph(edgelist_as_df, edgelist_as_a_list_of_dicts, seq_map, cfg, out_dir, children_df=None):
    """
    Builds the annotated network object in the .cx format
    (https://cytoscape.org/cx/cx2/specification/cytoscape-exchange-format-specification-(version-2)/)
    """
    base_cx = []
    init = {}
    init["numberVerification"] = [{"longNumber": 281474976710655}]
    base_cx.append(init)

    metadata = {}
    metadata["metaData"] = [
        {"name": "cyTableColumn", "version": "1.0"},
        {"name": "nodes", "version": "1.0"},
        {"name": "edges", "version": "1.0"},
        {"name": "nodeAttributes", "version": "1.0"},
        {"name": "edgeAttributes", "version": "1.0"},
        {"name": "networkAttributes", "version": "1.0"},
        {"name": "cartesianLayout", "version": "1.0"},
        {"name": "cyVisualProperties", "version": "1.0"}
        # {"name": "cyGroups", "version": "1.0"},
        # {"name": "cyHiddenAttributes", "version": "1.0"},
        # {"name": "cyNetworkRelations", "version": "1.0"},
        # {"name": "cySubNetworks", "version": "1.0"}
    ]
    base_cx.append(metadata)

    # GET ALL COLUMNS OF ALL TABLES
    table_columns = {}
    table_columns["cyTableColumn"] = []
    cyTableColumns = [
        # for node table mandatory
        {"applies_to": "node_table", "n": "@id", "d": "string"},
        {"applies_to": "node_table", "n": "name", "d": "string"},
        {"applies_to": "node_table", "n": "shared name", "d": "string"},
        {"applies_to": "node_table", "n": "display name", "d": "string"},
        {"applies_to": "node_table", "n": "microbetag::taxon name", "d": "string"},
        {"applies_to": "node_table", "n": "microbetag::namespace"},
        {"applies_to": "node_table", "n": "microbetag::taxonomy", "d": "string"},
        {"applies_to": "node_table", "n": "microbetag::ncbi-tax-id", "d": "string"},
        {"applies_to": "node_table", "n": "microbetag::ncbi-tax-level", "d": "string"},
        {"applies_to": "node_table", "n": "microbetag::gtdb-genomes", "d": "list_of_string"},
        # for edge table mandatory
        {"applies_to": "edge_table", "n": "shared name"},
        # {"applies_to": "edge_table", "n": "shared interaction"},
        {"applies_to": "edge_table", "n": "name"},
        {"applies_to": "edge_table", "n": "interaction type"},
        {"applies_to": "edge_table", "n": "weight::weight", "d": "double"}  # flashweave score
    ]

    """CHILDREN GENOMES AND NCBI IDS"""
    if cfg["get_children"]:
        cyTableColumns.append({"applies_to": "node_table", "n": "microbetag::children-ncbi-ids", "d": "list_of_string"})
        cyTableColumns.append({"applies_to": "node_table", "n": "microbetag::children-gtdb-genomes", "d": "list_of_string"})

    """Phen traits"""
    if cfg["phenDB"]:
        df_traits_table = pd.read_csv(os.path.join(out_dir, "phen_predictions/phen_traits.tsv"), sep="\t")
        phen_traits = [c for c in df_traits_table.columns.to_list() if "Score" not in c]
        phen_traits.remove("NCBI_ID")
        phen_traits.remove("gtdb_id")
        for term in phen_traits:
            cyTableColumns.append({"applies_to": "node_table", "n": "::".join(["phendb", term]), "d": "boolean"})
            cyTableColumns.append({"applies_to": "node_table", "n": "::".join(["phendbScore", "".join([term, "Score"])]), "d": "double"})
            cyTableColumns.append({"applies_to": "node_table", "n": "::".join(["phendbScore", "".join([term, "Std"])]), "d": "double"})
        df_traits_table["NCBI_ID"] = pd.to_numeric(df_traits_table["NCBI_ID"], errors='coerce')

    """FAPROTAX traits"""
    if cfg["faprotax"]:
        assignments_per_seqId = seqId_faprotax_functions_assignment(os.path.join(out_dir, "faprotax/sub_tables/"))
        faprotax_traits = list(flatten_list(assignments_per_seqId.values()))
        for term in faprotax_traits:
            column = {"applies_to": "node_table", "n": "::".join(["faprotax", term]), "d": "boolean"}
            cyTableColumns.append(column)

    """COMPLEMENTS"""
    if cfg["pathway_complement"]:
        complements_dict = json.load(open(os.path.join(out_dir, "path_compl/complements.json")))

        descrps = pd.read_csv(os.path.join(MAPPINGS, "module_descriptions"), sep="\t", header=None)
        descrps.columns= ["category", "moduleId", "description"]
        column_order = ["moduleId", "description", "category"]
        descrps = descrps[column_order]
        for n_pair, g_pair in complements_dict.items():
            for genomes, compls in g_pair.items():
                pairs_cats = set()  # NOTE use it if needed as last entry for each genome pair so cytoscapeApp can build filtering
                for index, compl in enumerate(compls):
                    try:
                        triplet = descrps[descrps["moduleId"] == compl[0][0]].values.tolist()[0]
                        pairs_cats.add(triplet[2])
                        extend = triplet + compl[0][1:]
                        complements_dict[n_pair][genomes][index] = [extend]
                    except:
                        logging.warn("Module with id:", compl[0][0], " not in current the module_descriptions file.")
                        extend = [compl[0][0]] + ["N/A", "N/A",] + compl[0][1:]
                        complements_dict[n_pair][genomes][index] = [extend]


        # ast.literal_eval function to evaluate the key as a Python literal expression and converts it into an object.
        # to safely evaluate a string containing a Python literal or container display (e.g., a dictionary or list) without running arbitrary code. 
        # It's a safer alternative to eval when dealing with untrusted input.
        complements_keys = [ast.literal_eval(x) for x in list(complements_dict.keys())]
        complements_dict = dict(zip(complements_keys, complements_dict.values()))

        for edge in edgelist_as_a_list_of_dicts:
            if (
                (seq_map[seq_map["seqId"] == edge["node_a"]]["gtdb_gen_repr"].notna()).any() and
                (seq_map[seq_map["seqId"] == edge["node_b"]]["gtdb_gen_repr"].notna()).any()
            ):
                genomes_a = seq_map[seq_map["seqId"] == edge["node_a"]]["gtdb_gen_repr"].to_list()[0]
                genomes_b = seq_map[seq_map["seqId"] == edge["node_b"]]["gtdb_gen_repr"].to_list()[0]
                for genome_a in genomes_a:
                    for genome_b in genomes_b:
                        edge_col = {"applies_to": "edge_table", "n": "".join(["compl::", genome_a, ":", genome_b]), "d": "list_of_string"}
                        cyTableColumns.append(edge_col)
                        # example of edge_col:  {"applies_to": "edge_table", "n": "compl::GCF_000174815.1:GCF_003253005.1", "d": "list_of_string"}

    """SEED SCORES"""
    if cfg["seed_scores"]:
        seed_scores_dict = json.load(open(os.path.join(out_dir, "seed_scores/seed_scores.json")))
        seed_scores_keys = [ast.literal_eval(x) for x in list(seed_scores_dict.keys())]
        seed_scores_dict = dict(zip(seed_scores_keys, seed_scores_dict.values()))
        cyTableColumns.append({"applies_to": "edge_table", "n": "::".join(["seed", "competition"]), "d": "double"})
        cyTableColumns.append({"applies_to": "edge_table", "n": "::".join(["seed", "competition-std"]), "d": "double"})
        cyTableColumns.append({"applies_to": "edge_table", "n": "::".join(["seed", "cooperation"]), "d": "double"})
        cyTableColumns.append({"applies_to": "edge_table", "n": "::".join(["seed", "cooperation-std"]), "d": "double"})

    """MANTA CLUSTERS"""
    if cfg["manta"]:
        cartesianLayout = {}; cartesianLayout["cartesianLayout"] = []
        m1 = {"applies_to": "node_table", "n": "::".join(["manta", "cluster"]), "d": "double"}
        m2 = {"applies_to": "node_table", "n": "::".join(["manta", "assignment"]), "d": "string"}
        cyTableColumns.append(m1)
        cyTableColumns.append(m2)
        manta_output_file = "/".join([out_dir, 'manta_annotated.cyjs'])
        manta_net = read_cyjson(manta_output_file)
        clusters = list(manta_net.nodes(data="cluster"))
        assignments = list(manta_net.nodes(data="assignment"))
        positions = list(manta_net.nodes(data="position"))
        manta_annotations = {}
        for pair in clusters:
            manta_annotations[pair[0]] = {}
            manta_annotations[pair[0]]["cluster"] = pair[1]
        for pair in assignments:
            manta_annotations[pair[0]]["assignment"] = pair[1]
        for pair in positions:
            manta_annotations[pair[0]]["position"] = pair[1]

    table_columns["cyTableColumn"] = cyTableColumns
    base_cx.append(table_columns)

    # NETWORK TABLE
    networkAttributes = {}
    networkAttributes["networkAttributes"] = []
    networkAttributes["networkAttributes"].append({"n": "database", "v": "microbetagDB"})
    networkAttributes["networkAttributes"].append({"n": "network type", "v": "annotated network of microbial co-occurrences"})
    networkAttributes["networkAttributes"].append({"n": "name", "v": "microbetag network"})
    networkAttributes["networkAttributes"].append({"n": "uri", "v": "https://hariszaf.github.io/microbetag/"})
    networkAttributes["networkAttributes"].append({"n": "version", "v": "1.0"})
    base_cx.append(networkAttributes)

    # NODES TABLE
    nodes = {}; nodes["nodes"] = []
    nodeAttributes = {}; nodeAttributes["nodeAttributes"] = []

    set_of_nodes = set(edgelist_as_df["node_a"].to_list() + edgelist_as_df["node_b"].to_list())
    node_counter = 1000
    seq_to_nodeID = {}

    for seq in set_of_nodes:
        node_counter += 1
        case = seq_map[seq_map["seqId"] == seq]
        node = {"@id": node_counter, "n": seq}
        seq_to_nodeID[seq] = node_counter
        nodes["nodes"].append(node)

        # NCBI attribute
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "@id", "v": seq, "d": "string"})
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "shared name", "v": seq, "d": "string"})
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "display name", "v": seq, "d": "string"})
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::taxon name", "v": case["extendedSpecies"].item(), "d": "string"})
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::taxonomy", "v": case["taxonomy"].item(), "d": "string"})

        # Clusters
        if cfg["manta"]:
            if seq in manta_annotations.keys():
                nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "manta::cluster", "v": str(manta_annotations[seq]["cluster"]), "d": "double"})
                nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "manta::assignment", "v": manta_annotations[seq]["assignment"], "d": "string"})
                cartesianLayout["cartesianLayout"].append({"node": node_counter, "x": manta_annotations[seq]["position"]["x"], "y": manta_annotations[seq]["position"]["y"]})

        if isinstance(case["ncbi_tax_id"].item(), pd._libs.missing.NAType):
            continue

        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::ncbi-tax-id", "v": str(case["ncbi_tax_id"].item()), "d": "string"})
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::ncbi-tax-level", "v": case["ncbi_tax_level"].item(), "d": "string"})

        if not isinstance(case["gtdb_gen_repr"].item(), float):
            nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::gtdb-genomes", "v": case["gtdb_gen_repr"].item(), "d": "list_of_string"})

        if cfg["get_children"]:
            ch_case = children_df[children_df["parent_ncbi_tax_id"] == case["ncbi_tax_id"].item()].dropna()
            if len(ch_case["child_ncbi_tax_id"].to_list()) > 0:
               children_ncbi_ids = [str(x) for x in list(flatten_list(ch_case["child_ncbi_tax_id"].to_list(), []))]
               nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::children-ncbi-ids", "v": children_ncbi_ids, "d": "list_of_string"})
            if len(ch_case["gtdb_gen_repr"].to_list()):
                children_genomes = list(flatten_list(ch_case["gtdb_gen_repr"].to_list(), []))
                nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::children-gtdb-genomes", "v": children_genomes, "d": "list_of_string"})

        if cfg["phenDB"]:
            df_traits_table["NCBI_ID_str"] = df_traits_table["NCBI_ID"].astype(str) 
            traits_case = df_traits_table[ df_traits_table["NCBI_ID_str"] ==  case["ncbi_tax_id"].item() ]

            if traits_case.empty:
                continue
            for trait in phen_traits:
                trait_value = list(set(traits_case[trait].to_list()))
                trait_score = traits_case["".join([trait, "Score"])].mean()
                trait_std = traits_case["".join([trait, "Score"])].std()
                if len(trait_value) == 1:
                    trait_value = True if trait_value[0] == "YES" else False
                    nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "::".join(["phendb", trait]), "v": trait_value, "d": "boolean"})
                    nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "::".join(["phendbScore", "".join([trait, "Score"])]), "v": str(trait_score), "d": "double"})
                    if traits_case.shape[0] == 1:
                        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "::".join(["phendbScore", "".join([trait, "Std"])]), "v": "0.0", "d": "double"})                    
                    else:
                        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "::".join(["phendbScore", "".join([trait, "Std"])]), "v": str(trait_std), "d": "double"})

        if cfg["faprotax"]:
            if seq in assignments_per_seqId:
                for term in assignments_per_seqId[seq]:
                    nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "::".join(["faprotax", term]), "v": True, "d": "boolean"})

    base_cx.append(nodes); base_cx.append(nodeAttributes)
    if cfg["manta"]:
        base_cx.append(cartesianLayout)

    # ADD EDGES AND THEIR ATTRIBUTES
    # [NOTE] Each seqId corresponds to a single node id, however each co-occurrence association from Flashweave or any other tool leads up to 3 edges in
    # the microbetag network; that is because in the case of an mspecies-mspecies association, we build 2 extra edges to describe this pairwised (A-B) association:
    # one where we have A as source and B as target and another that is vice-versa; we do that, so it is clear which complements have A as beneficiary species and B as donor and
    # the other way around. The same applies for the seed scores.
    # [NOTE] With respect to the path complementarities, we need to keep in mind that species A may benefit from species B but not the other way around.
    # Thus, we may have only one of the two extra edges.
    edges = {}; edges["edges"] = []
    edgeAttributes = {}; edgeAttributes["edgeAttributes"] = []
    edge_counter = node_counter + 1000
    for case in edgelist_as_a_list_of_dicts:
        id_a = seq_to_nodeID[case["node_a"]]
        id_b = seq_to_nodeID[case["node_b"]]
        edge = {"@id": edge_counter, "s": id_a, "t": id_b, "i": "cooccurs"}
        edges["edges"].append(edge)

        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "weight::weight", "v": str(case["score"]), "d": "double"})
        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "shared name", "v": " ".join([case["node_a"], "(cooccurss with)", case["node_b"]]), "d": "string"})
        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "interaction type", "v": "cooccurs with", "d": "string"})

        ncbi_a = seq_map[seq_map["seqId"] == case["node_a"]]["ncbi_tax_id"].item()
        ncbi_b = seq_map[seq_map["seqId"] == case["node_b"]]["ncbi_tax_id"].item()
        if isinstance(ncbi_a, pd._libs.missing.NAType) or isinstance(ncbi_b, pd._libs.missing.NAType):
            continue

        if ncbi_a == "<NA>" or ncbi_b == "<NA>":
            continue
        ncbi_pair_as_tuple_a_b = (int(ncbi_a), int(ncbi_b))
        ncbi_pair_as_tuple_b_a = (int(ncbi_b), int(ncbi_a))
        edge_counter += 1

        if cfg["pathway_complement"] or cfg["seed_scores"]:

            check = False

            """Edge for A -> B"""
            pot_edge = {"@id": (edge_counter), "s": id_a, "t": id_b, "i": "comp_coop"}

            # Path complements A -> B
            if cfg["pathway_complement"]:
                if ncbi_pair_as_tuple_a_b in complements_keys:
                    edges["edges"].append(pot_edge)
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "shared name", "v": " ".join([case["node_a"], "(completes/competes with)", case["node_b"]]), "d": "string"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "interaction type", "v": "completes/competes with", "d": "string"})
                    check = True
                    for genome_combo, complements in complements_dict[ncbi_pair_as_tuple_a_b].items():
                        genome_combo = ast.literal_eval(genome_combo)
                        attr = "".join(["compl::", genome_combo[0], ":", genome_combo[1]])
                        merged_compl = ["^".join(gcompl[0]) for gcompl in complements]
                        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": attr, "v": merged_compl, "d": "list_of_string"})

            # Seed scores A -> B
            if cfg["seed_scores"]:
                if ncbi_pair_as_tuple_a_b in seed_scores_keys:
                    if pot_edge not in edges["edges"]:
                        edges["edges"].append(pot_edge)                    
                        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "shared name", "v": " ".join([case["node_a"], "(completes/competes with)", case["node_b"]]), "d": "string"})
                        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "interaction type", "v": "completes/competes with", "d": "string"})
                        check = True

                    comp = []; coop = []
                    for index in seed_scores_dict[ncbi_pair_as_tuple_a_b].values():
                        if "competition" in index.keys():
                            comp.append(index["competition"])
                            coop.append(index["cooperation"])
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "competition"]), "v": str(np.mean(comp)), "d": "double"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "competition-std"]), "v": str(np.std(comp)), "d": "double"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "cooperation"]), "v": str(np.mean(coop)), "d": "double"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "cooperation-std"]), "v": str(np.std(coop)), "d": "double"})

            if check:
                edge_counter += 1
                check = False

            """Edge for B -> A"""
            pot_edge = {"@id": (edge_counter), "s": id_b, "t": id_a, "i": "comp_coop"}

            # Path complements B -> A
            if cfg["pathway_complement"]:
                if ncbi_pair_as_tuple_b_a in complements_keys:
                    edges["edges"].append(pot_edge)
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "shared name", "v": " ".join([case["node_b"], "(completes/competes with)", case["node_a"]]), "d": "string"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "interaction type", "v": "completes/competes with", "d": "string"})
                    check = True
                    for genome_combo, complements in complements_dict[ncbi_pair_as_tuple_b_a].items():
                        genome_combo = ast.literal_eval(genome_combo)
                        merged_compl = ["^".join(gcompl[0]) for gcompl in complements]
                        attr = "".join(["compl::", genome_combo[0], ":", genome_combo[1]])
                        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": attr, "v": merged_compl, "d": "list_of_string"})

            # Seed scores B -> A
            if cfg["seed_scores"]:
                if ncbi_pair_as_tuple_b_a in seed_scores_keys:
                    if pot_edge not in edges["edges"]:
                        edges["edges"].append(pot_edge)
                        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "shared name", "v": " ".join([case["node_b"], "(completes/competes with)", case["node_a"]]), "d": "string"})
                        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "interaction type", "v": "completes/competes with", "d": "string"})
                        check = True
                    comp = []; coop = []
                    for index in seed_scores_dict[ncbi_pair_as_tuple_b_a].values():
                        if "competition" in index.keys():
                            comp.append(index["competition"])
                            coop.append(index["cooperation"])
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "competition"]), "v": str(np.mean(comp)), "d": "double"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "competition-std"]), "v": str(np.std(comp)), "d": "double"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "cooperation"]), "v": str(np.mean(coop)), "d": "double"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "cooperation-std"]), "v": str(np.std(coop)), "d": "double"})

            if check:
                edge_counter += 1

    base_cx.append(edges); base_cx.append(edgeAttributes)

    visProp = {}
    visProp["cyVisualProperties"] = [{
        "properties_of": "network",
        "properties": {"NETWORK_ANNOTATION_SELECTION": "false", "NETWORK_BACKGROUND_PAINT": "#FFFFFF", "NETWORK_CENTER_X_LOCATION": "0.26885154030537706", "NETWORK_CENTER_Y_LOCATION": "4.089405406605124", "NETWORK_CENTER_Z_LOCATION": "0.0", "NETWORK_DEPTH": "0.0", "NETWORK_EDGE_SELECTION": "true", "NETWORK_FORCE_HIGH_DETAIL": "false", "NETWORK_HEIGHT": "856.0", "NETWORK_NODE_LABEL_SELECTION": "false", "NETWORK_NODE_SELECTION": "true", "NETWORK_SCALE_FACTOR": "2.532920631819882", "NETWORK_SIZE": "550.0", "NETWORK_WIDTH": "1879.0"}},
        {"properties_of": "nodes:default",
         "properties": {"COMPOUND_NODE_PADDING": "10.0", "COMPOUND_NODE_SHAPE": "ROUND_RECTANGLE", "NODE_BORDER_PAINT": "#CCCCCC", "NODE_BORDER_STROKE": "SOLID", "NODE_BORDER_TRANSPARENCY": "255", "NODE_BORDER_WIDTH": "0.0", "NODE_CUSTOMGRAPHICS_1": "org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove Graphics ], ", "NODE_CUSTOMGRAPHICS_2": "org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove Graphics ], ", "NODE_CUSTOMGRAPHICS_3": "org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove Graphics ], ", "NODE_CUSTOMGRAPHICS_4": "org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove Graphics ], ", "NODE_CUSTOMGRAPHICS_5": "org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove Graphics ], ", "NODE_CUSTOMGRAPHICS_6": "org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove Graphics ], ", "NODE_CUSTOMGRAPHICS_7": "org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove Graphics ], ", "NODE_CUSTOMGRAPHICS_8": "org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove Graphics ], ", "NODE_CUSTOMGRAPHICS_9": "org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove Graphics ], ", "NODE_CUSTOMGRAPHICS_POSITION_1": "C,C,c,0.00,0.00", "NODE_CUSTOMGRAPHICS_POSITION_2": "C,C,c,0.00,0.00", "NODE_CUSTOMGRAPHICS_POSITION_3": "NE,C,c,0.00,0.00", "NODE_CUSTOMGRAPHICS_POSITION_4": "C,C,c,0.00,0.00", "NODE_CUSTOMGRAPHICS_POSITION_5": "C,C,c,0.00,0.00", "NODE_CUSTOMGRAPHICS_POSITION_6": "C,C,c,0.00,0.00", "NODE_CUSTOMGRAPHICS_POSITION_7": "C,C,c,0.00,0.00", "NODE_CUSTOMGRAPHICS_POSITION_8": "C,C,c,0.00,0.00", "NODE_CUSTOMGRAPHICS_POSITION_9": "C,C,c,0.00,0.00", "NODE_CUSTOMGRAPHICS_SIZE_1": "50.0", "NODE_CUSTOMGRAPHICS_SIZE_2": "50.0", "NODE_CUSTOMGRAPHICS_SIZE_3": "50.0", "NODE_CUSTOMGRAPHICS_SIZE_4": "50.0", "NODE_CUSTOMGRAPHICS_SIZE_5": "50.0", "NODE_CUSTOMGRAPHICS_SIZE_6": "50.0", "NODE_CUSTOMGRAPHICS_SIZE_7": "50.0", "NODE_CUSTOMGRAPHICS_SIZE_8": "50.0", "NODE_CUSTOMGRAPHICS_SIZE_9": "50.0", "NODE_CUSTOMPAINT_1": "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_1, name=Node Custom Paint 1)", "NODE_CUSTOMPAINT_2": "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_2, name=Node Custom Paint 2)", "NODE_CUSTOMPAINT_3": "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_3, name=Node Custom Paint 3)", "NODE_CUSTOMPAINT_4": "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_4, name=Node Custom Paint 4)", "NODE_CUSTOMPAINT_5": "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_5, name=Node Custom Paint 5)", "NODE_CUSTOMPAINT_6": "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_6, name=Node Custom Paint 6)", "NODE_CUSTOMPAINT_7": "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_7, name=Node Custom Paint 7)", "NODE_CUSTOMPAINT_8": "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_8, name=Node Custom Paint 8)", "NODE_CUSTOMPAINT_9": "DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_9, name=Node Custom Paint 9)", "NODE_DEPTH": "0.0", "NODE_FILL_COLOR": "#C0C0C0", "NODE_HEIGHT": "45.0", "NODE_LABEL_COLOR": "#000000", "NODE_LABEL_FONT_FACE": "SansSerif.plain,plain,12", "NODE_LABEL_FONT_SIZE": "12", "NODE_LABEL_POSITION": "C,C,c,0.00,0.00", "NODE_LABEL_ROTATION": "0.0", "NODE_LABEL_TRANSPARENCY": "255", "NODE_LABEL_WIDTH": "200.0", "NODE_NESTED_NETWORK_IMAGE_VISIBLE": "true", "NODE_PAINT": "#1E90FF", "NODE_SELECTED": "false", "NODE_SELECTED_PAINT": "#FFFF00", "NODE_SHAPE": "ELLIPSE", "NODE_SIZE": "35.0", "NODE_TRANSPARENCY": "255", "NODE_VISIBLE": "true", "NODE_WIDTH": "45.0", "NODE_X_LOCATION": "0.0", "NODE_Y_LOCATION": "0.0", "NODE_Z_LOCATION": "0.0"},
         "dependencies": {"nodeCustomGraphicsSizeSync": "true", "nodeSizeLocked": "false"},
         "mappings": {"NODE_CUSTOMGRAPHICS_1": {"type": "PASSTHROUGH", "definition": "COL=stringdb::STRING style,T=string"}, "NODE_CUSTOMGRAPHICS_2": {"type": "PASSTHROUGH", "definition": "COL=stringdb::chemViz Passthrough,T=string"}, "NODE_CUSTOMGRAPHICS_3": {"type": "PASSTHROUGH", "definition": "COL=stringdb::enhancedLabel Passthrough,T=string"}, "NODE_CUSTOMGRAPHICS_POSITION_3": {"type": "DISCRETE", "definition": "COL=stringdb::node type,T=string,K=0=protein,V=0=NE,,C,,c,,0.00,,0.00,K=1=compound,V=1=N,,C,,c,,0.00,,-5.00"}, "NODE_FILL_COLOR": {"type": "DISCRETE", "definition": "COL=name,T=string,K=0=4932.YPR119W,V=0=#FF008B,K=1=4932.YLR079W,V=1=#F380FF,K=2=4932.YJL076W,V=2=#5D00FF,K=3=4932.YOR195W,V=3=#00B9FF,K=4=4932.YFR028C,V=4=#8097FF,K=5=4932.YMR055C,V=5=#80FFDC,K=6=4932.YHR152W,V=6=#00FF2E,K=7=4932.YIL106W,V=7=#E8FF00,K=8=4932.YGR092W,V=8=#AEFF80,K=9=4932.YAR019C,V=9=#FFC580,K=10=4932.YML064C,V=10=#FF0000"}, "NODE_HEIGHT": {"type": "DISCRETE", "definition": "COL=stringdb::node type,T=string,K=0=protein,V=0=50.0,K=1=compound,V=1=40.0"}, "NODE_SHAPE": {"type": "DISCRETE", "definition": "COL=stringdb::node type,T=string,K=0=protein,V=0=ELLIPSE,K=1=compound,V=1=ROUND_RECTANGLE"}, "NODE_TRANSPARENCY": {"type": "DISCRETE", "definition": "COL=stringdb::node type,T=string,K=0=protein,V=0=255,K=1=compound,V=1=0"}, "NODE_WIDTH": {"type": "DISCRETE", "definition": "COL=stringdb::node type,T=string,K=0=protein,V=0=50.0,K=1=compound,V=1=100.0"}}},
        {"properties_of": "edges:default",
         "properties": {"EDGE_CURVED": "true", "EDGE_LABEL_COLOR": "#000000", "EDGE_LABEL_FONT_FACE": "Dialog.plain,plain,10", "EDGE_LABEL_FONT_SIZE": "10", "EDGE_LABEL_ROTATION": "0.0", "EDGE_LABEL_TRANSPARENCY": "255", "EDGE_LABEL_WIDTH": "200.0", "EDGE_LINE_TYPE": "SOLID", "EDGE_PAINT": "#323232", "EDGE_SELECTED": "false", "EDGE_SELECTED_PAINT": "#FF0000", "EDGE_SOURCE_ARROW_SELECTED_PAINT": "#FFFF00", "EDGE_SOURCE_ARROW_SHAPE": "NONE", "EDGE_SOURCE_ARROW_SIZE": "6.0", "EDGE_SOURCE_ARROW_UNSELECTED_PAINT": "#000000", "EDGE_STACKING": "AUTO_BEND", "EDGE_STACKING_DENSITY": "0.5", "EDGE_STROKE_SELECTED_PAINT": "#FF0000", "EDGE_STROKE_UNSELECTED_PAINT": "#1F293D", "EDGE_TARGET_ARROW_SELECTED_PAINT": "#FFFF00", "EDGE_TARGET_ARROW_SHAPE": "NONE", "EDGE_TARGET_ARROW_SIZE": "6.0", "EDGE_TARGET_ARROW_UNSELECTED_PAINT": "#000000", "EDGE_TRANSPARENCY": "255", "EDGE_UNSELECTED_PAINT": "#404040", "EDGE_VISIBLE": "true", "EDGE_WIDTH": "2.0", "EDGE_Z_ORDER": "0.0"}, "dependencies": {"arrowColorMatchesEdge": "false"}, "mappings": {"EDGE_TRANSPARENCY": {"type": "CONTINUOUS", "definition": "COL=stringdb::score,T=double,L=0=34,E=0=34,G=0=34,OV=0=0.2,L=1=85,E=1=85,G=1=85,OV=1=0.5,L=2=170,E=2=170,G=2=170,OV=2=1.0"}, "EDGE_WIDTH": {"type": "CONTINUOUS", "definition": "COL=stringdb::score,T=double,L=0=0.8,E=0=0.8,G=0=0.8,OV=0=0.2,L=1=2.0,E=1=2.0,G=1=2.0,OV=1=0.5,L=2=4.0,E=2=4.0,G=2=4.0,OV=2=1.0"}}}
    ]
    base_cx.append(visProp)

    # POST-metadata
    post_metadata = {}
    post_metadata["metaData"] = []
    post_metadata["metaData"].append({"name": "nodeAttributes", "elementCount": len(nodeAttributes["nodeAttributes"]), "version": 1.0})
    post_metadata["metaData"].append({"name": "edgeAttributes", "elementCount": len(edgeAttributes["edgeAttributes"]), "version": 1.0})
    post_metadata["metaData"].append({"name": "cyTableColumn", "elementCount": len(table_columns["cyTableColumn"]), "version": 1.0})
    post_metadata["metaData"].append({"name": "edges", "elementCount": len(edges["edges"]), "idCounter": node_counter + 1000, "version": 1.0})
    post_metadata["metaData"].append({"name": "nodes", "elementCount": len(nodes["nodes"]), "idCounter": 1001, "version": 1.0})
    post_metadata["metaData"].append({"name": "cyVisualProperties","elementCount":3,"version":"1.0"})
    post_metadata["metaData"].append({"name": "networkPropernetworkAttributesties","elementCount": len(networkAttributes["networkAttributes"]), "version":"1.0"})

    if cfg["manta"]:
        post_metadata["metaData"].append({"name": "cartesianLayout", "elementCount": len(cartesianLayout["cartesianLayout"]), "version": 1.0})

    base_cx.append(post_metadata)

    # Status
    status = {}; status["status"] = []
    status["status"].append({"error":"","success": True})
    base_cx.append(status)


    return base_cx
