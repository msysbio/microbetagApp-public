"""
Parse user's input OTU table
"""
from .db_functions import *
import pandas as pd
import numpy as np
import logging
import os
import sys
from .variables import *
import json
from fuzzywuzzy import fuzz
from multiprocessing import Pool
# import networkx as nx
# import difflib


def init_run():
    """
    fix it so we have a route for cytoscape app
    another one for terminal where user gets a zip file
    """


def annotate_network(base_network, config, ids_map, out_dir):
    """
    [TODO] add seed scores
    """
    annotated_network = base_network.copy()

    if config["phenDB"]:

        df_traits_table = pd.read_csv(os.path.join(out_dir, "phen_predictions/phen_traits.tsv"), sep="\t")
        for node in annotated_network["elements"]["nodes"]:
            ncbiId = node["data"]["NCBI-Tax-Id"]
            if ncbiId == "<NA>":
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
 
        seed_scores_dict = json.load(open(os.path.join(out_dir, "")))


    return annotated_network


def export_species_level_associations(edgelist_as_a_list_of_dicts):
    """
    This functions gets as input the edges of the network as a list of dictionaries
    checks the ncbi_tax_level of each node and in case where for both nodes it is species or strain
    gets their corresponding GTDB genomes on microbetagDB
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

    ncbiId_to_genomesIds_queries = []
    for ncbiId in set_of_ncbiids_of_interest:
        query = "".join(["SELECT genomeId from genome2taxNcbiId where ncbiTaxId = ", str(ncbiId), ";"])
        ncbiId_to_genomesIds_queries.append(query)

    # Execute the queries
    genomes = []
    cnx, cursor = create_cursor()
    for query in ncbiId_to_genomesIds_queries:
        cursor.execute(query)
        query_results = cursor.fetchall()
        genomes.append({query.split(" ")[-1][:-1]: [x[0] for x in query_results]})

    related_genomes = {key: value for dictionary in genomes for key, value in dictionary.items()}

    return pairs_of_interest, related_genomes


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


def map_seq_to_ncbi_tax_level_and_id(
        abundance_table, my_taxonomy_column, seqId, taxonomy_scheme, tax_delim):
    """
    Parse user's OTU table and the Silva database to get to add 2 extra columns in the OTU table:
    1. the lowest taxonomic level of the taxonomy assigned in an OTU, for which an NCBI Taxonomy id exists (e.g., "genus")
    2. the corresponding NCBI Taxonomy Id (e.g., "343")

    The SPECIES_NCBI_ID, GENERA_NCBI_IDS and FAMILIES_NCBI_IDS files are like this:
    Species  ncbi_tax_id
    D_6__Abiotrophia sp. oral clone OH2A       319434

    Returns:
    a. seq_id_taxid_level_repr_genome: a pandas df with the sequence id, the ncbi tax id of the lowest taxonomic level in the taxonomy
         and their corresponding GTDB genome if available
    b. repr_genomes_present: a list of GTDB genomes
    """
    if taxonomy_scheme == "GTDB":
        taxon_to_ncbiId = os.path.join(MAPPINGS, "gtdbSpecies2ncbiId2accession.tsv")

        gtdb_accession_ids = pd.read_csv(taxon_to_ncbiId, sep="\t")
        gtdb_accession_ids.columns = ["Species", "ncbi_tax_id", "gtdb_gen_repr"]
        gtdb_accession_ids["Species"].str.strip()

        species_ncbi_id = pd.read_csv(SPECIES_NCBI_IDS, sep="\t")
        species_ncbi_id.columns = ['Species', 'ncbi_tax_id']
        species_ncbi_id['Species'].str.strip()

    else:
        taxon_to_ncbiId = os.path.join(MAPPINGS, "species2ncbiId2accession.tsv")
        taxon_to_ncbiId = pd.read_csv(taxon_to_ncbiId, sep="\t", names = ["name", "ncbi", "gtdb"])

    # Split the taxonomy column and split it based on semicolumn!
    taxonomies = abundance_table.filter([seqId, my_taxonomy_column, "microbetag_id"])
    splitted_taxonomies = taxonomies[my_taxonomy_column].str.split(tax_delim, expand=True)
    splitted_taxonomies.columns = ["Domain", "Phylum", "Class", "Order", "Family","Genus", "Species"]
    splitted_taxonomies[seqId] = taxonomies[seqId]
    splitted_taxonomies["microbetag_id"] = taxonomies["microbetag_id"]

    splitted_taxonomies_c = splitted_taxonomies.copy()

    # Check if "s__" taxonomy like and whether species in 1 or 2 steps
    if splitted_taxonomies['Species'].str.contains("s__").all():
        pattern_for_complete_name = r"__[\s_]" 
        if splitted_taxonomies['Species'].str.contains(pattern_for_complete_name).all():
            splitted_taxonomies["extendedSpecies"] = splitted_taxonomies['Species'].apply(process_species)
        else:
            # We need to use the Genus column too
            genus = splitted_taxonomies['Genus'].apply(process_species)
            species = splitted_taxonomies['Species'].apply(process_species)
            splitted_taxonomies["extendedSpecies"] = pd.DataFrame({'extendedSpecies': np.where(species.isna(), None, genus + ' ' + species)})
        unique_species = splitted_taxonomies['extendedSpecies'].unique()            
    else:
        unique_species = splitted_taxonomies['Species'].unique()

    unique_species = [item for item in unique_species if item is not None]

    # Run a pool to get NCBI ids at SPECIES level
    chunk_size = round(len(unique_species)/2)
    chunks = [unique_species[i:i+chunk_size] for i in range(0, len(unique_species), chunk_size)]

    pool = multiprocessing.Pool(2) 
    data = [ (chunk, taxon_to_ncbiId) for chunk in chunks]
    ncbi_ids_species_level = pool.map(calculate_fuzzy_similarity_chunk, data)

    ncbi_ids_species_name_matched_df = pd.DataFrame.from_dict(ncbi_ids_species_level[0], orient="index" )
    ncbi_ids_species_name_df = ncbi_ids_species_name_matched_df["ncbi"]
    ncbi_ids_species_name = dict(ncbi_ids_species_name_df)

    splitted_taxonomies["species_ncbi_id"]=splitted_taxonomies["extendedSpecies"].map(ncbi_ids_species_name)


    # Get NCBI ids at GENUS level
    genera_to_get_ids = splitted_taxonomies.loc[splitted_taxonomies['species_ncbi_id'].isna(), 'Genus']





# =============================================================



    genera_ncbi_id = pd.read_csv(GENERA_NCBI_IDS, sep="\t")
    genera_ncbi_id.columns = ['Genus', 'ncbi_tax_id']
    genera_ncbi_id['Genus'].str.strip()

    family_ncbi_id = pd.read_csv(FAMILIES_NCBI_IDS, sep="\t")
    family_ncbi_id.columns = ['Family', 'ncbi_tax_id']
    family_ncbi_id['Family'].str.strip()

    """
    Build a dataframe for the Species, Genus and Family taxonomies present
    on the OTU table along with their corresponding NCBI Tax IDs
    """



    """
    GTDB Ids will be added in another function
    """
    # # GTDB accession ids
    # genomes_present = gtdb_accession_ids.merge(
    #     splitted_taxonomies,
    #     on=["Species"]
    # )[["ncbi_tax_id", "gtdb_gen_repr", seqId]]

    # splitted_taxonomies = pd.merge(genomes_present, splitted_taxonomies,
    #                                on=seqId, how='outer')
    # splitted_taxonomies.loc[splitted_taxonomies["ncbi_tax_id"].notnull(
    # ), "ncbi_tax_level"] = "mspecies"

    # mspecies = splitted_taxonomies[[seqId, "gtdb_gen_repr"]]



    # # Species
    # species_present = species_ncbi_id.merge(splitted_taxonomies.query('ncbi_tax_level != "mspecies"'),
    #                                         on=['Species'],
    #                                         suffixes=('_species', '_mspecies')
    #                                         )[
    #     [seqId, "ncbi_tax_id_species", "ncbi_tax_level"]
    # ]

    species_present.loc[species_present["ncbi_tax_level"] !=
                        "mspecies", "ncbi_tax_level"] = "species"

    splitted_taxonomies = pd.merge(splitted_taxonomies, species_present,
                                   on=seqId,
                                   how="outer",
                                   suffixes=('_species', '_mspecies')
                                   )

    splitted_taxonomies["ncbi_tax_id"] = splitted_taxonomies["ncbi_tax_id"].fillna(
        splitted_taxonomies["ncbi_tax_id_species"])
    splitted_taxonomies["ncbi_tax_level"] = splitted_taxonomies["ncbi_tax_level_mspecies"].fillna(
        splitted_taxonomies["ncbi_tax_level_species"])

    splitted_taxonomies = splitted_taxonomies.drop(
        ["ncbi_tax_id_species", "ncbi_tax_level_species", "ncbi_tax_level_mspecies"], axis=1)

    """
    REMEMBER!
    There is no one-to-one relationship between a NCBI Taxonomy Id and a representative genome
    We do not have a hit for ALL NCBI TAX IDs from the species found in our experiment, aka there s
    no GTDB representative genome for all NCBI Taxonomy Ids, e.g.  Bradyrhizobium sp. J81, NCBI Tax Id: 656743
    """

    # Genera
    pd_genera = splitted_taxonomies[[
        seqId, "ncbi_tax_level", "ncbi_tax_id", "gtdb_gen_repr", "Genus"]]
    genera_present = genera_ncbi_id.merge(
        pd_genera, on=['Genus'], suffixes=(
            "_gen", "_over"), how="right")

    genera_present.loc[genera_present["ncbi_tax_id_gen"].notnull(
    ), "ncbi_tax_level_gen"] = "genus"

    genera_present["ncbi_tax_id"] = genera_present['ncbi_tax_id_over'].combine_first(
        genera_present['ncbi_tax_id_gen'])
    genera_present["ncbi_tax_level"] = genera_present['ncbi_tax_level'].combine_first(
        genera_present['ncbi_tax_level_gen'])

    genera_present = genera_present.drop(
        ["ncbi_tax_id_gen", "ncbi_tax_id_over", "Genus"], axis=1)

    # Families
    pd_families = splitted_taxonomies[[seqId, "Family"]]
    families_present = family_ncbi_id.merge(
        pd_families, on=['Family'], how="right")
    families_present.loc[families_present["ncbi_tax_id"].notnull(
    ), "ncbi_tax_level"] = "family"

    families_present["ncbi_tax_id"] = genera_present["ncbi_tax_id"].combine_first(
        families_present["ncbi_tax_id"])
    families_present["ncbi_tax_level"] = genera_present["ncbi_tax_level"].combine_first(
        families_present["ncbi_tax_level"])
    families_present = families_present.drop(["Family"], axis=1)

    # Build a unified data frame for all levels and the accession ids when
    # available
    otu_taxid_level_repr_genome = pd.merge(
        splitted_taxonomies_c,
        families_present,
        on=seqId)
    otu_taxid_level_repr_genome = pd.merge(
        otu_taxid_level_repr_genome,
        mspecies,
        on=seqId,
        how="outer")

    otu_taxid_level_repr_genome['ncbi_tax_id'] = otu_taxid_level_repr_genome['ncbi_tax_id'].astype(pd.Int64Dtype()).astype(str)

    repr_genomes_present = list(mspecies["gtdb_gen_repr"].dropna())

    return otu_taxid_level_repr_genome, repr_genomes_present


def calculate_fuzzy_similarity_chunk(args_set):
    """
    Gets a list of taxa names and checks what is their best hit against the 
    """
    import time, re, sys
    buzz_taxa = ["uncultured bacterium", "uncultured organism", "uncultured beta proteobacterium",
        "uncultured soil bacterium", "uncultured rumen bacterium", "uncultured gamma proteobacterium",
        "D_8__uncultured rumen protozoa", "uncultured delta proteobacterium"]
    chunk, taxon_to_ncbiId_df = args_set[0], args_set[1]
    start = time.time()
    result = []
    _shared_dict = {}
    for taxon in chunk:
        taxon = taxon.strip()
        _shared_dict[taxon] = {}
         # check if it is not a species name at all
        tokens =  re.split(r'[ _]', taxon)
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
                start2 = time.time()
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
                    print(taxon)
                    print(best_match)
                    print(best_similarity)
                    _shared_dict[taxon]["ncbi"] = None
                end2 = time.time()
    end = time.time()
    print("\n\nTotal chunk time: ", str(end-start))
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


def get_child_taxa(parent_tax_id, nodes_file):
    """
    Returns all children NCBI Taxonomy IDs of a parent.
    """
    child_taxa = []
    with open(nodes_file, 'r') as nodes:
        for line in nodes:
            fields = line.strip().split('\t|\t')
            tax_id, parent_id, rank = fields[0], fields[1], fields[2]
            if parent_id == parent_tax_id:
                child_taxa.append(tax_id)
    return child_taxa


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
        new_edge["data"]["flashweave-score"] = edge["score"]
        new_edge["selected"] = False

        edges.append(new_edge)
        counter += 1

    base_network["elements"]["nodes"] = nodes
    base_network["elements"]["edges"] = edges

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

    associated_pairs_node_a = pd.merge(pd_edgelist[["node_a", "joint"]], abundance_table_with_ncbi_ids[["ncbi_tax_id", "gtdb_gen_repr", "ncbi_tax_level", "microbetag_id"]],
                                       left_on='node_a', right_on='microbetag_id', how="inner").drop(["microbetag_id"], axis=1)

    associated_pairs_node_a.rename(columns={
        "ncbi_tax_level": "ncbi_tax_level_node_a",
        "gtdb_gen_repr": "gtdb_gen_repr_node_a",
        "ncbi_tax_id": "ncbi_tax_id_node_a"
    }, inplace=True)

    associated_pairs_node_b = pd.merge(pd_edgelist[["node_b", "joint", "score"]], abundance_table_with_ncbi_ids[["ncbi_tax_id", "gtdb_gen_repr", "ncbi_tax_level", "microbetag_id"]],
                                       left_on='node_b', right_on='microbetag_id', how="inner").drop(["microbetag_id"], axis=1)

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


def process_species(entry):
    """
    Convert "s__" like taxa to ncbi like ones
    """
    parts = entry.split("__")
    if len(parts[1]) > 0:
        return parts[1].replace("_", " ")
    else:
        return np.nan


def tuple_to_str(t):
    return ','.join(map(str, t))





if __name__ == "__main__":

from microbetag.scripts.utils import *
import pandas as pd
from fuzzywuzzy import fuzz
import multiprocessing
import os
from multiprocessing import Pool
abundance_table = is_tab_separated("static/input_test/k__taxonomy.tsv", "taxonomy")  #dada2_use_case
ext = ensure_flashweave_format(abundance_table, "taxonomy", "ASV_ID", ".")
taxonomies = abundance_table.filter(["ASV_ID","taxonomy", "microbetag_id"])
splitted_taxonomies = taxonomies["taxonomy"].str.split(";", expand=True)
splitted_taxonomies.columns = ["Domain","Phylum", "Class", "Order","Family","Genus", "Species"]
splitted_taxonomies["ASV_ID"] = taxonomies["ASV_ID"]
taxon_to_ncbiId = os.path.join("microbetag/mappings/", "species2ncbiId.tsv")
taxon_to_ncbiId = pd.read_csv(taxon_to_ncbiId, sep="\t", names = ["name", "ncbi"])
# Check if "s__" taxonomy like and whether species in 1 or 2 steps
if splitted_taxonomies['Species'].str.contains("s__").all():
    pattern_for_complete_name = r"__[\s_]" 
    if splitted_taxonomies['Species'].str.contains(pattern_for_complete_name).all():
        splitted_taxonomies["extendedSpecies"] = splitted_taxonomies['Species'].apply(process_species)
    else:
        # We need to use the Genus column too
        genus = splitted_taxonomies['Genus'].apply(process_species)
        species = splitted_taxonomies['Species'].apply(process_species)
        splitted_taxonomies["extendedSpecies"] = pd.DataFrame({'extendedSpecies': np.where(species.isna(), None, genus + ' ' + species)})
    unique_species = splitted_taxonomies['extendedSpecies'].unique()            
else:
    unique_species = splitted_taxonomies['Species'].unique()

unique_species = [item for item in unique_species if item is not None]


chunk_size = round(len(unique_species)/2)
chunks = [unique_species[i:i+chunk_size] for i in range(0, len(unique_species), chunk_size)]

pool = multiprocessing.Pool(2) 
data = [ (chunk, taxon_to_ncbiId) for chunk in chunks]
ncbi_ids_species_level = pool.map(calculate_fuzzy_similarity_chunk, data)

ncbi_ids_species_name_matched_df = pd.DataFrame.from_dict(ncbi_ids_species_level[0], orient="index" )
ncbi_ids_species_name_df = ncbi_ids_species_name_matched_df["ncbi"]
ncbi_ids_species_name = dict(ncbi_ids_species_name_df)


splitted_taxonomies["species_ncbi_id"]=splitted_taxonomies["extendedSpecies"].map(ncbi_ids_species_name)



