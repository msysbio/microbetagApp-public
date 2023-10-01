import sys
import os
import decimal
import pandas as pd
import mysql.connector
from .variables import *


# General
def create_cursor():
    cnx = mysql.connector.connect(user=USER_NAME, password=PASSWORD, host=HOST, database=DB_NAME)
    cnx.get_warnings = True
    cursor = cnx.cursor()
    return cnx, cursor


def query_to_microbetagDB(phrase):
    """
    Function to get functional traits as assigned in a genome from the phenotrex software
    phenotrex.readthedocs.io/ and stored in the microbetagDB.
    """
    try:
        cnx = mysql.connector.connect(user=USER_NAME, password=PASSWORD, host=HOST, database=DB_NAME)
        cnx.get_warnings = True
        cursor = cnx.cursor()
        cursor.execute(phrase)

        res = []
        for row in cursor:
            res.append(row)

        warnings = cursor.fetchwarnings()
        if warnings:
            for i in range(len(warnings)):
                print("\t ** Warning - " + warnings[i][2])
        cursor.close()
        cnx.commit()
        cnx.close()
        return res

    except mysql.connector.Error as err:
        print("Something went wrong: {}".format(err))
        print(phrase)


def get_column_names():
    """
    Get the column names of a database table
    """
    phrase = "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='" +\
             DB_NAME + "' AND `TABLE_NAME`='phenDB';"
    colnames = query_to_microbetagDB(phrase)
    colnames = [x[0] for x in colnames]
    return colnames


def execute(phrase):
    """
    Establish a database connection and perform an action
    """
    # Database connection configuration
    config = {
        'user': 'msysbio',
        'password': 'pass',
        'host': 'db',  # Replace with the actual database service name or IP address
        'port': 3306,  # Replace with the actual port
        'database': 'microbetagDB'
    }
    cnx = mysql.connector.connect(**config)
    cursor = cnx.cursor()
    cursor.execute(phrase)
    rows = cursor.fetchall()
    cursor.close()
    cnx.commit()
    cnx.close()
    return rows


# Mapping
def get_genomes_for_ncbi_tax_id(ncbi_tax_id=1281578):
    """
    Get the genomes IDs that correspond to a NCBI Taxonomy Id and are present in the microbetagDB
    """
    query = "".join(
        ["SELECT genomeId from genome2taxNcbiId where ncbiTaxId = ", str(ncbi_tax_id), ";"])
    genome_ids = execute(query)

    return {ncbi_tax_id: [x[0] for x in genome_ids]}


# Phen related
def get_phendb_traits(gtdb_genome_id="GCA_018819265.1"):
    """
    Get phenotypical traits based on phenDB classes based on its GTDB representative genome
    """
    query = "".join(
        ["SELECT * FROM phenDB WHERE gtdbId = '", gtdb_genome_id, "';"])

    rows = execute(query)

    if len(rows) == 0:
        return 0

    query_colnames = "SHOW COLUMNS FROM phenDB;"
    colnames = [list(x)[0] for x in execute(query_colnames)]
    genomes_traits = {i: j for i, j in zip(colnames, rows[0])}

    return genomes_traits


# Pathway complementarity related
def get_complements_for_pair_of_genomes(
        beneficiarys_genome_id="GCA_003184265.1", donors_genome_id="GCA_000015645.1"):
    """
    For a certain beneficiary genome and a certain donor genome, retrieve all complements
    available in the database.
    """
    q = query_for_getting_compl_ids(beneficiarys_genome_id, donors_genome_id)
    # complements_ids_list_query = "".join(['SELECT complmentId FROM pathwayComplementarity WHERE beneficiaryGenome = "',
    #                                       beneficiarys_genome_id,
    #                                       '" AND donorGenome = "',
    #                                       donors_genome_id,
    #                                       '";'
    #                                       ])
    export = execute(q)
    if len(export) > 0:
        complements_ids_list = export[0][0].split(",")
    else:
        return
    complements = []
    for complementId in complements_ids_list:
        query = "".join(
            [
                "SELECT KoModuleId, complement, pathway FROM uniqueComplements WHERE complementId = '",
                complementId,
                "';"])
        compl = execute(query)
        complements.append(compl)

    return complements


def get_complements_for_pair_of_ncbiIds(
        beneficiarys_ndbi_id=1281578, donors_ncbi_id=146891):
    """
    Gets a pair of NCBI Taxonomy Ids and returns and all the potential complements,
    based on all the corresponding genomes.
    """
    beneficiarys_genomes = get_genomes_for_ncbi_tax_id(beneficiarys_ndbi_id)
    donors_genomes = get_genomes_for_ncbi_tax_id(donors_ncbi_id)
    complements = {}
    for beneficiary_genome in list(beneficiarys_genomes.values())[0]:
        complements[beneficiary_genome] = {}
        for donor_genome in list(donors_genomes.values())[0]:
            pair_compl = get_complements_for_pair_of_genomes(
                beneficiary_genome, donor_genome)
            if pair_compl is not None:
                colored_pair_compl = build_kegg_urls(pair_compl)
                complements[beneficiary_genome][donor_genome] = colored_pair_compl

    return complements


def get_complements_of_list_of_pair_of_ncbiIds(pairs_of_interest, relative_genomes):
    """
    Gets a list of dictionaries as input that describes the edges of a network and returns the complementarities
    as a dictionary where a pair of ncbi ids is the key and a list with complementarities the value
    e.g.: ()
    This function is running as part of the microbetag pipeline while the others mostly support the microbetag API.
    By running chunks of queries using the same cursor, we save quite some time.
    """
    # set_of_ncbiids_of_interest = set()
    # list_of_non_usefules_ids = ["77133"]
    # pairs_of_interest = set()
    # for pair in ncbi_id_pairs:
    #     taxon_a = str(pair["ncbi_tax_id_node_a"]).split(".")[0]
    #     taxon_b = str(pair["ncbi_tax_id_node_b"]).split(".")[0]
    #     if taxon_a in list_of_non_usefules_ids or taxon_b in list_of_non_usefules_ids:
    #         continue
    #     # we keep as a pair both a->b and b->a association
    #     if pair["ncbi_tax_level_node_a"] == pair["ncbi_tax_level_node_b"] == "mspecies":
    #         set_of_ncbiids_of_interest.add(taxon_a)
    #         set_of_ncbiids_of_interest.add(taxon_b)
    #         pairs_of_interest.add((taxon_a, taxon_b))
    #         pairs_of_interest.add((taxon_b, taxon_a))

    # ncbiId_to_genomesIds_queries = []
    # for ncbiId in set_of_ncbiids_of_interest:
    #     query = "".join(["SELECT genomeId from genome2taxNcbiId where ncbiTaxId = ", str(ncbiId), ";"])
    #     ncbiId_to_genomesIds_queries.append(query)

    # # Execute the queries
    # genomes = []
    # cnx, cursor = create_cursor()
    # for query in ncbiId_to_genomesIds_queries:
    #     cursor.execute(query)
    #     query_results = cursor.fetchall()
    #     genomes.append({query.split(" ")[-1][:-1]: [x[0] for x in query_results]})

    # relative_genomes = {key: value for dictionary in genomes for key, value in dictionary.items()}
    # # relative_genomes = {}
    # # for dictionary in genomes:
    # #     for key, value in dictionary.items():
    # #         relative_genomes[key] = value

    # Build the queries
    complements_ids_queries = {}
    for pair in list(pairs_of_interest):
        ncbi_a = pair[0]
        ncbi_b = pair[1]
        for ncbi_a_genome in relative_genomes[ncbi_a]:
            for ncbi_b_genome in relative_genomes[ncbi_b]:
                q = query_for_getting_compl_ids(ncbi_a_genome, ncbi_b_genome)
                if pair in complements_ids_queries:
                    complements_ids_queries[pair].append(q)
                else:
                    complements_ids_queries[pair] = [q]

    # Queries to the database to get complementarities
    cnx, cursor = create_cursor()
    pairs_to_compl_ids = {}
    for pair, queries in complements_ids_queries.items():
        for query in queries:
            cursor.execute(query)
            query_result = cursor.fetchall()
            if len(query_result) > 0:
                complements_ids_list = query_result[0][0].split(",")
                """
                [TODO] Check if correct; probably only last genome pair is returned for each ncbi id
                """
                pairs_to_compl_ids[pair] = complements_ids_list

    # Queries to map unique complements ids to their extended. human readable version
    cnx, cursor = create_cursor()
    pairs_complements = {}
    for pair, complementIds in pairs_to_compl_ids.items():
        for complementId in complementIds:
            query = "".join(
                ["SELECT KoModuleId, complement, pathway FROM uniqueComplements WHERE complementId = '", complementId, "';"])
            cursor.execute(query)
            query_result = cursor.fetchall()
            colored_pair_compl = build_kegg_urls([query_result])[0][0]
            if pair in pairs_complements:
                pairs_complements[pair].append(colored_pair_compl)
            else:
                pairs_complements[pair] = [colored_pair_compl]

    return pairs_complements


def query_for_getting_compl_ids(beneficiary, donor):
    return "".join(['SELECT complmentId FROM pathwayComplementarity WHERE beneficiaryGenome = "', beneficiary,
                    '" AND donorGenome = "', donor, '";'])


def build_kegg_urls(complements_for_a_pair_of_genomes):
    """
    Takes as input the complements list between two genomes and
    build urls to colorify the related to the module kegg map based on the KO terms of the beneficiary (pink)
    and those it gets from the donor (green).
    """
    # Load the dictionary with the kegg modules and their corresponding maps
    color_mapp_base_url = "https://www.kegg.jp/kegg-bin/show_pathway?"
    present_kos_color = "%09%23EAD1DC/"
    complemet_kos_color = "%09%2300A898/"
    maps = open("/home/app/web/static/module_map_pairs.tsv", "r")
    module_map = {}
    for line in maps:
        module, mmap = line.split("\t")
        module_map[module[3:-1]] = mmap[1:-1]

    # Make a url pointing at a colored kegg map based on what's on the
    # beneficiary's genome and what it gets as complement from the donor
    tmp_compl = complements_for_a_pair_of_genomes.copy()
    for index, compl in enumerate(tmp_compl):
        compl = list(compl[0])
        md = compl[0]
        actual_compl = compl[1].split(";")
        path = compl[2].split(";")

        beneficiarys_kos = ""
        complements_kos = ""
        for ko_term in path:
            if ko_term not in actual_compl:
                beneficiarys_kos = "".join(
                    [beneficiarys_kos, ko_term, present_kos_color])
            else:
                complements_kos = "".join(
                    [complements_kos, ko_term, complemet_kos_color])

        # In rare cases, the module might not have a related map
        try:
            url_ko_map_colored = "".join(
                [color_mapp_base_url, module_map[md], "/", beneficiarys_kos, complements_kos])
        except KeyError:
            url_ko_map_colored = "N/A"

        compl.append(url_ko_map_colored)

        complements_for_a_pair_of_genomes[index][0] = compl

    return complements_for_a_pair_of_genomes


# Seed scores related
def get_seed_scores_for_list_of_pair_of_ncbiIds(pairs_of_interest, relative_genomes):

    # Build the queries
    seed_scores_queries = {}
    for pair in list(pairs_of_interest):
        ncbi_a = pair[0]
        ncbi_b = pair[1]
        for ncbi_a_genome in relative_genomes[ncbi_a]:
            for ncbi_b_genome in relative_genomes[ncbi_b]:
                q = query_for_getting_seed_scores(ncbi_a_genome, ncbi_b_genome)
                if pair in seed_scores_queries:
                    seed_scores_queries[pair].append(q)
                else:
                    seed_scores_queries[pair] = [q]

    # Query to the database
    cnx, cursor = create_cursor()
    pairs_to_compl_ids = {}
    for pair, queries in seed_scores_queries.items():
        for query in queries:
            cursor.execute(query)
            query_result = cursor.fetchall()
            if len(query_result) > 0:
                complements_ids_list = query_result[0][0].split(",")
                pairs_to_compl_ids[pair] = complements_ids_list
    """
    [TODO] add once db ready
    """

    return


def get_seed_scores_for_pair_of_ncbiIds(ncbiId_A, ncbiId_B):
    genomes_for_species_A = get_genomes_for_ncbi_tax_id(ncbiId_A)
    genomes_for_species_B = get_genomes_for_ncbi_tax_id(ncbiId_B)
    seeds_A_B = []
    seeds_B_A = []
    for genome_A in list(genomes_for_species_A.values())[0]:
        for genome_B in list(genomes_for_species_B.values())[0]:
            case_seeds_A_B, case_seeds_B_A = get_seed_scores_for_pair_of_genomes(genome_A, genome_B)
            seeds_A_B.append(case_seeds_A_B)
            seeds_B_A.append(case_seeds_B_A)

    return seeds_A_B, seeds_B_A


def get_seed_scores_for_pair_of_genomes(genome_A="GCA_003184265.1", genome_B="GCA_000015645.1"):

    query_from_A_to_B, query_from_B_to_A = query_for_getting_seed_scores(genome_A, genome_B)

    export = execute(query_from_A_to_B)
    competition_score_A_B, complementarity_score_A_B = export[0][0], export[0][1]

    export = execute(query_from_B_to_A)
    competition_score_B_A, complementarity_score_B_A = export[0][0], export[0][1]

    seeds_A_B = genome_A, genome_B, competition_score_A_B, complementarity_score_A_B
    seeds_B_A = genome_B, genome_A, competition_score_B_A, complementarity_score_B_A

    return seeds_A_B, seeds_B_A


def query_for_getting_seed_scores(genome_A, genome_B):

    query_from_A_to_B = "".join([
        'SELECT Competition, Complementarity FROM seedScores WHERE genomeA = "',
        genome_A,
        '" AND genomeB = "',
        genome_B,
        '";'
    ])

    query_from_B_to_A = "".join([
        'SELECT Competition, Complementarity FROM seedScores WHERE genomeA = "',
        genome_B,
        '" AND genomeB = "',
        genome_A,
        '";'
    ])

    return query_from_A_to_B, query_from_B_to_A
