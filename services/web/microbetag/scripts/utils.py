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
        # taxon_a and taxon_b variables are int type; in case there is a <NA> value then it is a pandas missing NAType
        taxon_a = seqID_taxid_level_repr_genome[seqID_taxid_level_repr_genome["seqId"] == pair["node_a"]]["ncbi_tax_id"].item()
        taxon_b = seqID_taxid_level_repr_genome[seqID_taxid_level_repr_genome["seqId"] == pair["node_b"]]["ncbi_tax_id"].item()
        taxon_a_level = seqID_taxid_level_repr_genome[seqID_taxid_level_repr_genome["seqId"] == pair["node_a"]]["ncbi_tax_level"].item()
        taxon_b_level = seqID_taxid_level_repr_genome[seqID_taxid_level_repr_genome["seqId"] == pair["node_b"]]["ncbi_tax_level"].item()
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
    return pairs_of_ncbi_id_of_interest, related_genomes, children_genomes


def count_comment_lines(my_abd_tab, tax_col):
    """
    Get the number of rows of the OTU table that should be skipped
    """
    skip_rows = 0
    with open(my_abd_tab, 'r') as f:
        for line in f:
            if line.startswith('#') and tax_col not in line:
                skip_rows += 1
            elif tax_col in line:
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


def custom_agg(series):
    """
    Aggregate rows that share the same value in a column of the df.
    """
    result = []
    for value in series.dropna().unique():
        if isinstance(value, tuple):
            result.extend(value)
        else:
            result.append(value)
    return result or None


def replace_empty_list(value):
    """
    Replace a df's cell that is an empty list with np.nan
    """
    return np.nan if isinstance(value, list) and not value else value


def map_seq_to_ncbi_tax_level_and_id(abd_tab, tax_col, seqId, tax_scheme, tax_delim, get_chiildren):
    """
    Parse user's OTU table and the Silva database to get to add 2 extra columns in the OTU table:
    1. the lowest taxonomic level of the taxonomy assigned in an OTU, for which an NCBI Taxonomy id exists (e.g., "genus")
    2. the corresponding NCBI Taxonomy Id (e.g., "343")

    Returns:
    splt_tax: a pandas df with the sequence id, the ncbi tax id of species, genus, family level when available
                        and the lowest taxonomic level for which there is an ncbi tax id
                        also a list (as column in th df) with gtdb genomes if available
                        species_ncbi_id genus_ncbi_id family_ncbi_id  tax_level  ncbi_tax_id,
    repr_genomes_present: a list with the GTDB genomes of interest
    children_df: (optional)
    """
    # Split the taxonomy column and split it based on semicolumn!
    taxonomies = abd_tab.filter([seqId, tax_col, "microbetag_id"])
    splt_tax = taxonomies[tax_col].str.split(tax_delim, expand=True)
    if splt_tax[0].str.contains('Root').all():
        splt_tax = splt_tax.drop(splt_tax.columns[0], axis=1)

    if len(list(splt_tax.columns)) != 7:
        raise ValueError(f'{"Error: The taxonomy scheme provided is not a 7-level one."}')

    splt_tax.columns = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    splt_tax[seqId] = taxonomies[seqId]
    splt_tax["microbetag_id"] = taxonomies["microbetag_id"]

    # Check if "s__" taxonomy-like and whether species in 1 or 2 steps
    s_pattern = "s__"
    if splt_tax['Species'].str.contains(s_pattern).all():
        underscore = True
        pattern_for_complete_name = r"__[\s_]"
        if splt_tax['Species'].str.contains(pattern_for_complete_name).all():
            splt_tax["extendedSpecies"] = splt_tax['Species'].apply(process_underscore_taxonomy)
        else:
            # We might need to use the Genus column too, check if genus is part of species
            genus = splt_tax['Genus'].apply(process_underscore_taxonomy)
            species = splt_tax['Species'].apply(process_underscore_taxonomy)
            genus_list = genus.tolist(); genus_list = [str(c) for c in genus_list]
            species_list = species.tolist(); species_list = [str(c) for c in species_list]
            results = [True if genus in species else False for genus, species in zip(genus_list, species_list) if species != 'nan']
            # Check if genus is included in species for non-empty species
            if results.count(True) > 0.8 * len(results):
                splt_tax["extendedSpecies"] = pd.DataFrame({'extendedSpecies': np.where(species.isna(), None, species)})
            else:
                splt_tax["extendedSpecies"] = pd.DataFrame({'extendedSpecies': np.where(species.isna(), None, genus + ' ' + species)})
    else:
        splt_tax["extendedSpecies"] = splt_tax['Species']
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
    if tax_scheme == "GTDB":
        taxon_to_ncbiId = os.path.join(MAPPINGS, "gtdbSpecies2ncbiId2accession.tsv")

    elif tax_scheme == "Silva":
        taxon_to_ncbiId = os.path.join(MAPPINGS, "gtdb_silvaSpecies2ncbi2accession.tsv")

    elif tax_scheme == "microbetag_prep":
        taxon_to_ncbiId = os.path.join(MAPPINGS, "gc_accession_16s_gtdb_ncbid.tsv")

    else:
        taxon_to_ncbiId = os.path.join(MAPPINGS, "species2ncbiId2accession.tsv")

    # Get species NCBI ids if you have one of the 3 schemes supported
    if tax_scheme == "GTDB" or tax_scheme == "Silva" or tax_scheme == "microbetag_prep":
        gtdb_accession_ids = pd.read_csv(taxon_to_ncbiId, sep="\t")
        gtdb_accession_ids.columns = ["refSpecies", "species_ncbi_id", "gtdb_gen_repr"]
        gtdb_accession_ids["refSpecies"].str.strip()

        splt_tax = pd.merge(splt_tax, gtdb_accession_ids, left_on='extendedSpecies', right_on='refSpecies', how='left')
        repr_genomes_present = [c for c in splt_tax["gtdb_gen_repr"].to_list() if isinstance(c, str)]

        if tax_scheme != "Silva":
            splt_tax['gtdb_gen_repr'] = splt_tax['gtdb_gen_repr'].apply(lambda x: [x] if pd.notna(x) else np.nan)

        else:
            # Remove cases where for a single seqId there are more than 2 ncbi taxids in case both hit the same genome
            splt_tax['gtdb_gen_repr'] = splt_tax['gtdb_gen_repr'].apply(lambda x: tuple([x]) if isinstance(x, str) else tuple())

            # Group by 'seqId' and aggregate 'gtdb_gen_repr' with custom function
            splt_tax = splt_tax.groupby('seqId').agg({
                **{col: 'first' for col in splt_tax.columns if col not in ['seqId', 'gtdb_gen_repr']},
                'gtdb_gen_repr': aggregate_genomes
            }).reset_index()

            splt_tax = splt_tax.map(replace_empty_list)

        gtdb_genomes_on = True

    # Get species NCBI ids if you have any other taxonomy scheme
    else:

        taxon_to_ncbiId_df = pd.read_csv(taxon_to_ncbiId, sep="\t", names=["name", "ncbi", "gtdb"])
        unique_species = splt_tax['extendedSpecies'].unique()
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

        splt_tax["species_ncbi_id"] = splt_tax["extendedSpecies"].map(ncbi_ids_species_name)

        gtdb_genomes_on = False

    """
    Get NCBI ids at GENUS level
    """
    genera_ncbi_id = pd.read_csv(GENERA_NCBI_IDS, sep="\t")
    genera_ncbi_id.columns = ['Genus', 'ncbi_tax_id']
    genera_ncbi_id['Genus'].str.strip()

    if underscore:
        genera = splt_tax['Genus'].apply(process_underscore_taxonomy)
    else:
        genera = splt_tax['Genus']
    splt_tax["extendedGenus"] = genera
    genera.name = "Genus"
    genera = genera.to_frame()
    genera = genera.merge(genera_ncbi_id, on='Genus', how='inner').drop_duplicates()
    splt_tax = pd.merge(splt_tax, genera, left_on='extendedGenus', right_on='Genus', how='left')
    splt_tax = splt_tax.drop('Genus_y', axis=1)
    splt_tax = splt_tax.rename(columns={'ncbi_tax_id': 'genus_ncbi_id'})
    splt_tax['genus_ncbi_id'] = np.where(splt_tax['species_ncbi_id'].notna(), np.nan, splt_tax['genus_ncbi_id'])

    """
    Get NCBI ids at FAMILY level
    """
    family_ncbi_id = pd.read_csv(FAMILIES_NCBI_IDS, sep="\t")
    family_ncbi_id.columns = ['Family', 'ncbi_tax_id']
    family_ncbi_id['Family'].str.strip()
    if underscore:
        families = splt_tax['Family'].apply(process_underscore_taxonomy)
    else:
        families = splt_tax['Family']
    splt_tax["extendedFamily"] = families
    families.name = "Family"
    families = families.to_frame()
    families = families.merge(family_ncbi_id, on='Family', how='inner').drop_duplicates()
    splt_tax = pd.merge(splt_tax, families, left_on='extendedFamily', right_on='Family', how='left')
    splt_tax = splt_tax.drop('Family_y', axis=1)
    splt_tax = splt_tax.rename(columns={'ncbi_tax_id': 'family_ncbi_id'})

    # Remove proxy taxon columns
    splt_tax = splt_tax.drop('extendedGenus', axis=1)
    splt_tax = splt_tax.drop('extendedFamily', axis=1)

    """
    Get GTDB genomes for taxa available in case where taxonomy scheme is "other"
    """
    if not gtdb_genomes_on:
        unique_species_present_ncbi_ids = [int(x) for x in list(splt_tax[splt_tax["species_ncbi_id"].notna()]["species_ncbi_id"])]

        # Dictionary with ncbi ids of species level nodes and their gtdb genomes
        species_ncbi_ids_to_gtdb_genomes = {}
        for ncbi_id in unique_species_present_ncbi_ids:
            genomes = get_genomes_for_ncbi_tax_id(ncbi_id)
            gc_genomes = [genome for genome in list(genomes.values())[0] if genome.startswith("GCA_") or genome.startswith("GCF_")]
            if len(gc_genomes) > 0:
                # species_ncbi_ids_to_gtdb_genomes[float(ncbi_id)] = gc_genomes
                species_ncbi_ids_to_gtdb_genomes[str(ncbi_id)] = gc_genomes

        # A list with the gtdb present on the dataset
        repr_genomes_present = [item for sublist in list(species_ncbi_ids_to_gtdb_genomes.values()) for item in sublist]

        # Now map those GTDB ids to their corresponding entries in the splitted_taxonomy df
        splt_tax['gtdb_gen_repr'] = splt_tax['species_ncbi_id'].map(species_ncbi_ids_to_gtdb_genomes)
        # {float(k): v[0] for k, v in species_ncbi_ids_to_gtdb_genomes.items()}

    # Keep track of the taxonomic level of the ncbi tax id of the lowest level found
    splt_tax['ncbi_tax_id'] = splt_tax.apply(assign_tax_id_for_node_level, axis=1).astype(str)
    splt_tax['ncbi_tax_id'] = pd.to_numeric(splt_tax['ncbi_tax_id'], errors='coerce').astype('Int64').astype(str)
    # splt_tax['ncbi_tax_id'] = splt_tax['ncbi_tax_id'].fillna('').astype(str) #.str.rstrip('.0').replace('', np.nan).astype(float).astype('Int64')

    # Keep what is the taxonomic level of the node
    splt_tax['ncbi_tax_level'] = splt_tax.apply(determine_tax_level, axis=1)

    # Remove any white spaces from the dataframe's columns
    splt_tax = splt_tax.map(lambda x: x.strip() if isinstance(x, str) else x)

    # Beautify df
    splt_tax = splt_tax.rename(columns={"Genus_x": "Genus", "Family_x": "Family"})
    desired_order = ["microbetag_id", seqId, "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "extendedSpecies",
                     "family_ncbi_id", "genus_ncbi_id", "species_ncbi_id", "ncbi_tax_id", "ncbi_tax_level", "gtdb_gen_repr"]
    splt_tax = splt_tax[desired_order]

    tmp_df = splt_tax; tmp_df2 = tmp_df 
    tmp_df["gtdb_gen_repr"] = tmp_df["gtdb_gen_repr"].apply(lambda x: tuple(x) if isinstance(x, list) else x)
    result_df = splt_tax.groupby('seqId').agg(
        {
            # 'ncbi_tax_id': custom_agg, # by now, we should always have a single ncbi id for each taxonomy. we have made sure of it on the mapping files
            'gtdb_gen_repr': custom_agg
        }
        ).reset_index()

    mapping_dict = dict(zip(result_df['seqId'], result_df['gtdb_gen_repr']))
    splt_tax['gtdb_gen_repr'] = splt_tax['seqId'].map(mapping_dict)

    # Get children NCBI Tax ids and their corresponding genomes
    if get_chiildren:
        ncbi_parent_to_children = {}
        ncbi_nodes_dict = get_ncbi_nodes_dict()
        species_df = splt_tax[(splt_tax['ncbi_tax_level'] == "species") | (splt_tax['ncbi_tax_level'] == "mspecies")]
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
        return splt_tax, repr_genomes_present, children_df

    else:
        return splt_tax, repr_genomes_present  #, tmp_df2


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


def is_tab_separated(my_abd_tab, tax_col):
    """
    Read the OTU table and make sure it is tab separated and not empty
    Takes as input a .tsv file and returns a pandas dataframe.
    """
    number_of_commented_lines = count_comment_lines(
        my_abd_tab, tax_col)
    try:
        abd_tab = pd.read_csv(
            my_abd_tab,
            # sep=ABUNDANCE_TABLE_DELIM,
            skiprows=number_of_commented_lines)
    except BaseException:
        logging.error(
            "The OTU table provided is not a tab separated file. Please convert your OTU table to .tsv or .csv format.")

    if abd_tab.shape[1] < 2:
        logging.error("The OTU table you provide has no records.")

    return abd_tab


def ensure_flashweave_format(my_abd_tab, tax_col, seqId, outdir):
    """
    Build an OTU table that will be in a FlashWeave-based format.
    """
    flashweave_table = my_abd_tab.drop(tax_col, axis=1)
    float_col = flashweave_table.select_dtypes(include=['float64'])

    for col in float_col.columns.values:
        flashweave_table[col] = flashweave_table[col].astype('int64')

    flashweave_table[seqId] = flashweave_table[seqId].astype(str)
    my_abd_tab['microbetag_id'] = flashweave_table[seqId]

    file_to_save = os.path.join(
        outdir,
        "abd_tab_flashweave_format.tsv")
    flashweave_table.to_csv(file_to_save, sep='\t', index=False)

    return my_abd_tab


def edge_list_of_ncbi_ids(edgelist, metadata_file=None):  # abd_tab_with_ncbi_ids
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
        mask = pd_edgelist.isin(elements_to_exclude).any(axis=1)
        pd_edgelist = pd_edgelist[~mask]

    pd_edgelist["pair-of-taxa"] = pd_edgelist['node_a'].astype(str) + ":" + pd_edgelist["node_b"]

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
        descrps.columns = ["category", "moduleId", "description"]
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
                    except ValueError:
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

        if not isinstance(case["gtdb_gen_repr"].item(), type(None)):
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
            traits_case = df_traits_table[df_traits_table["NCBI_ID_str"] == case["ncbi_tax_id"].item()]

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

    # POST-metadata
    post_metadata = {}
    post_metadata["metaData"] = []
    post_metadata["metaData"].append({"name": "nodeAttributes", "elementCount": len(nodeAttributes["nodeAttributes"]), "version": 1.0})
    post_metadata["metaData"].append({"name": "edgeAttributes", "elementCount": len(edgeAttributes["edgeAttributes"]), "version": 1.0})
    post_metadata["metaData"].append({"name": "cyTableColumn", "elementCount": len(table_columns["cyTableColumn"]), "version": 1.0})
    post_metadata["metaData"].append({"name": "edges", "elementCount": len(edges["edges"]), "idCounter": node_counter + 1000, "version": 1.0})
    post_metadata["metaData"].append({"name": "nodes", "elementCount": len(nodes["nodes"]), "idCounter": 1001, "version": 1.0})
    post_metadata["metaData"].append({"name": "networkPropernetworkAttributesties", "elementCount": len(networkAttributes["networkAttributes"]), "version": "1.0"})
    if cfg["manta"]:
        post_metadata["metaData"].append({"name": "cartesianLayout", "elementCount": len(cartesianLayout["cartesianLayout"]), "version": 1.0})

    base_cx.append(post_metadata)

    # Status
    status = {}; status["status"] = []
    status["status"].append({"error": "", "success": True})
    base_cx.append(status)

    return base_cx




"""
df = pd.read_csv("microbetag/mappings/species2ncbiId.tsv", sep="\t")
df.columns = ["taxon", "id"]
duplicates_mask = df.duplicated(subset="taxon", keep=False)
double_species = df[duplicates_mask]
doubles_dict = double_species.to_dict(orient="records")


counter = 0
taxon_id_genome_presence = {}
taxa_with_no_genomes = set()

doubles = {}
for entry in doubles_dict:
    counter += 1
    taxon, gid  = entry["taxon"], entry["id"]
    try:
        check = microbetag.get_genomes_for_ncbi_tax_id(gid)
    except:
        pass
    if isinstance(check, dict):
        for case in check.values():
            for scenario in case:
                if "GC" in scenario:
                    if taxon not in taxon_id_genome_presence:
                        taxon_id_genome_presence[taxon] = gid
                        doubles[taxon] = [gid]
                    else:
                        doubles[taxon].append(gid)
    if counter % 100 == 0:
        print(counter)


for entry in doubles_dict:
    taxon, gid  = entry["taxon"], entry["id"]
    if taxon not in doubles.keys():
        taxa_with_no_genomes.add(entry["taxon"])

g = open("FROMNET.tsv", "r") 
s = g.readlines()
pairs = {}
for case in s:
    gid = case.split("\t")[0]
    taxon = case.split("\t")[1][:-1]
    pairs[gid] = taxon


actually_doubles = {}
for taxon, genomes in doubles.items():
    if len(genomes) > 1:
        actually_doubles[taxon] = genomes




original_file = open("microbetag/mappings/species2ncbiId.tsv", "r")
with open("CLEANSPECIESNCBIT.TSV", "w") as f:
    for line in original_file:
        taxon, nid = line.split("\t")[0], line.split("\t")[1]
        taxon = taxon.strip()
        if taxon in taxa_with_no_genomes:
            continue
        if taxon not in taxon_id_genome_presence:
            f.write(line)
        elif taxon in actually_doubles:
            for genome in actually_doubles[taxon]:
                try:
                    print(pairs[str(int(genome))], genome)
                except:
                    print("not working")
                    pass
        else:
            print(">>", taxon)
            f.write(taxon + "\t" + str(int(taxon_id_genome_presence[taxon])) +"\n")

"""
