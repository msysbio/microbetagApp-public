import sys
import os
import decimal
import pandas as pd
import mysql.connector
from .variables import *


# General
def create_cursor():
    """
    Creates a connection to the microbetagDB; gets values from the variables.py.
    """
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


def get_column_names(db_table):
    """
    Get the column names of a database table
    """
    phrase = "".join([
        "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='",
        DB_NAME,
        "' AND `TABLE_NAME`='",
        db_table, "';"
    ])
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
    Returns a dictionary.
    """
    query = "".join(
        ["SELECT genomeId from genome2taxNcbiId where ncbiTaxId = ", str(ncbi_tax_id), ";"])
    genome_ids = execute(query)
    return {ncbi_tax_id: gc_unify([x[0] for x in genome_ids])}


# Phen related
def get_phendb_traits(gtdb_genome_id="GCA_018819265.1"):
    """
    Get phenotypical traits based on phenDB classes based on its GTDB representative genome
    """
    gtdb_genome_id = gtdb_genome_id.strip()
    query = "".join(["SELECT * FROM phenDB WHERE gtdbId = '", gtdb_genome_id, "';"])
    rows = execute(query)
    if len(rows) == 0:
        if gtdb_genome_id.startswith("GCA_"):
            gtdb_genome_id = gtdb_genome_id.replace("GCA_", "GCF_")
        elif gtdb_genome_id.startswith("GCF_"):
            gtdb_genome_id = gtdb_genome_id.replace("GCF_", "GCA_")
        else:
            # Handle cases where the input doesn't start with either prefix
            return 0
        query = "".join(["SELECT * FROM phenDB WHERE gtdbId = '", gtdb_genome_id, "';"])
        rows = execute(query)
        if len(rows) == 0:
            print("This is not a NCBI accession id. It could be a MGnify or a KEGG one.")
            return 0
    query_colnames = "SHOW COLUMNS FROM phenDB;"
    colnames = [list(x)[0] for x in execute(query_colnames)]
    genomes_traits = {i: j for i, j in zip(colnames, rows[0])}
    return genomes_traits  # gtdb_genome_id


# Pathway complementarity related
def get_complements_for_pair_of_genomes(beneficiarys_genome_id="GCA_003184265.1", donors_genome_id="GCA_000015645.1"):
    """
    For a certain beneficiary genome and a certain donor genome, retrieve all complements
    available in the database.
    Genome A is the beneficiary genome while genome B the one that provides the complement.
    """
    if not beneficiarys_genome_id.startswith("GC") or not donors_genome_id.startswith("GC"):
        return None
    q = query_for_getting_compl_ids(beneficiarys_genome_id, donors_genome_id)
    export = execute(q)
    if len(export) > 0:
        complements_ids_list = export[0][0].split(",")
    else:
        return
    complements = []
    for complementId in complements_ids_list:
        query = "".join([
            "SELECT KoModuleId, complement, pathway FROM uniqueComplements WHERE complementId = '",
            complementId,
            "';"
        ])
        compl = execute(query)
        complements.append(compl)
    coloured_complements = build_kegg_urls(complements)
    return coloured_complements


def get_complements_for_pair_of_ncbiIds(beneficiarys_ndbi_id=1281578, donors_ncbi_id=146891):
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
            pair_compl = get_complements_for_pair_of_genomes(str(beneficiary_genome), str(donor_genome))
            # if pair_compl is not None:
            #     colored_pair_compl = build_kegg_urls(pair_compl)
            #     complements[beneficiary_genome][donor_genome] = colored_pair_compl

    return complements


def get_complements_of_list_of_pair_of_ncbiIds(pairs_of_interest, relative_genomes):
    """
    Gets a list of dictionaries as input that describes the edges of a network and returns the complementarities
    as a dictionary where a pair of ncbi ids is the key and a list with complementarities the value
    e.g.: ()
    This function is running as part of the microbetag pipeline while the others mostly support the microbetag API.
    By running chunks of queries using the same cursor, we save quite some time.
    """
    # Build the queries
    complements_ids_queries = {}
    for ncbi_pair in list(pairs_of_interest):
        ncbi_a = ncbi_pair[0]
        ncbi_b = ncbi_pair[1]
        for ncbi_a_genome in relative_genomes[ncbi_a]:
            for ncbi_b_genome in relative_genomes[ncbi_b]:
                genomes_pair = (ncbi_a_genome, ncbi_b_genome)
                q = query_for_getting_compl_ids(ncbi_a_genome, ncbi_b_genome)
                if ncbi_pair in complements_ids_queries:
                    complements_ids_queries[ncbi_pair][genomes_pair] = q
                else:
                    complements_ids_queries[ncbi_pair] = {}
                    complements_ids_queries[ncbi_pair][genomes_pair] = q

    # Queries to the database to get complementarities
    cnx, cursor = create_cursor()
    pairs_to_compl_ids = {}
    for ncbi_pair, genomes_pair_query in complements_ids_queries.items():
        for genome_pair, query in genomes_pair_query.items():  # this should always be 1 case
            cursor.execute(query)
            query_result = cursor.fetchall()
            if len(query_result) > 0:
                complements_ids_list = query_result[0][0].split(",")
                if ncbi_pair not in pairs_to_compl_ids:
                    pairs_to_compl_ids[ncbi_pair] = {}
                elif genome_pair not in pairs_to_compl_ids:
                    pairs_to_compl_ids[ncbi_pair][genome_pair] = []

                pairs_to_compl_ids[ncbi_pair][genome_pair] = complements_ids_list

    # Queries to map unique complements ids to their extended. human readable version
    cnx, cursor = create_cursor()
    pairs_complements = {}
    for ncbi_pair, genomes_pair_complement in pairs_to_compl_ids.items():
        for pair_of_genomes, complementId_list in genomes_pair_complement.items():
            for complementId in complementId_list:
                query = "".join(["SELECT KoModuleId, complement, pathway FROM uniqueComplements WHERE complementId = '", complementId, "';"])
                cursor.execute(query)
                query_result = cursor.fetchall()
                colored_pair_compl = build_kegg_urls([query_result])[0][0]
                if ncbi_pair not in pairs_complements:
                    pairs_complements[ncbi_pair] = {}
                if pair_of_genomes not in pairs_complements[ncbi_pair]:
                    pairs_complements[ncbi_pair][pair_of_genomes] = []
                pairs_complements[ncbi_pair][pair_of_genomes].append(colored_pair_compl)

    return pairs_complements


def query_for_getting_compl_ids(beneficiary, donor):
    """
    Gets 2 gc accession ids and returns a query for their pathway complementarities
    """
    beneficiary_alt = alt_genome_prefix(beneficiary)
    donor_alt = alt_genome_prefix(donor)
    return "".join(['SELECT complmentId FROM pathwayComplementarity WHERE (beneficiaryGenome = "', str(beneficiary),
                    '" or beneficiaryGenome = "', beneficiary_alt,
                    '") AND (donorGenome = "', str(donor), '" OR donorGenome = "', donor_alt, '");'])


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
    """
    """
    # Build the queries
    seed_scores_queries = {}
    # Get all PATRIC ids corresponding to the relative_genomes
    flat_list = [item for sublist in list(relative_genomes.values()) for item in sublist]
    patricIds = get_patric_id_of_gc_accession_list(flat_list)
    for pair in list(pairs_of_interest):
        ncbi_a = pair[0]
        ncbi_b = pair[1]
        for ncbi_a_genome in gc_unify(relative_genomes[ncbi_a]):
            for ncbi_b_genome in gc_unify(relative_genomes[ncbi_b]):
                if patricIds[ncbi_a_genome] is None or patricIds[ncbi_b_genome] is None:
                    continue
                # [NOTE] We do not use query_from_B_to_A as all associations are double (both A->B and B->A) in the pairs_of_interest
                query_from_A_to_B, query_from_B_to_A = query_for_getting_seed_scores(patricIds[ncbi_a_genome], patricIds[ncbi_b_genome])
                if len(query_from_A_to_B) > 1:
                    if pair in seed_scores_queries:
                        seed_scores_queries[pair][len(seed_scores_queries[pair])] = {}
                        seed_scores_queries[pair][len(seed_scores_queries[pair]) - 1]["organism-A"] = ncbi_a
                        seed_scores_queries[pair][len(seed_scores_queries[pair]) - 1]["organism-B"] = ncbi_b
                        seed_scores_queries[pair][len(seed_scores_queries[pair]) - 1]["genome-A"] = ncbi_a_genome
                        seed_scores_queries[pair][len(seed_scores_queries[pair]) - 1]["patric-genome-A"] = patricIds[ncbi_a_genome]
                        seed_scores_queries[pair][len(seed_scores_queries[pair]) - 1]["genome-B"] = ncbi_b_genome
                        seed_scores_queries[pair][len(seed_scores_queries[pair]) - 1]["patric-genome-B"] = patricIds[ncbi_b_genome]
                        seed_scores_queries[pair][len(seed_scores_queries[pair]) - 1]["query"] = query_from_A_to_B
                    else:
                        seed_scores_queries[pair] = {}
                        seed_scores_queries[pair][0] = {}
                        seed_scores_queries[pair][0]["organism-A"] = ncbi_a
                        seed_scores_queries[pair][0]["organism-B"] = ncbi_b
                        seed_scores_queries[pair][0]["genome-A"] = ncbi_a_genome
                        seed_scores_queries[pair][0]["patric-genome-A"] = patricIds[ncbi_a_genome]
                        seed_scores_queries[pair][0]["genome-B"] = ncbi_b_genome
                        seed_scores_queries[pair][0]["patric-genome-B"] = patricIds[ncbi_b_genome]
                        seed_scores_queries[pair][0]["query"] = query_from_A_to_B
    # Query to the database
    cnx, cursor = create_cursor()
    for pair, cases in seed_scores_queries.items():
        for number_case, case in cases.items():
            cursor.execute(case["query"])
            seed_score = cursor.fetchall()

            if len(seed_score) > 0:
                competition, cooperation = float(seed_score[0][0]), float(seed_score[0][1])
                seed_scores_queries[pair][number_case]["competition"] = competition
                seed_scores_queries[pair][number_case]["cooperation"] = cooperation

    seed_scores = {}
    for k, v in seed_scores_queries.items():
        seed_scores[k] = {}
        # filtered_data = {l: m for l, m in v.items() if 'competition' in list(m.keys())}
        filtered_data = {rindex: rdata for rindex, rdata in v.items() if 'competition' in list(rdata.keys())}
        seed_scores[k] = filtered_data

    return seed_scores


def get_seed_scores_for_pair_of_ncbiIds(ncbiId_A, ncbiId_B):
    """
    Metabolic Complementarity Index is calculated as the fraction of A's seed set that is found within B's metabolic network
    but not part of B's seed set, normalized by the number of A's seed set in B's entire metabolic network.

    Metabolic Competition Index is calculated as the fraction of A's seed set that is also in B's seed set,
    normalized by the weighted sum of the confidence score.

    Gets a pair of NCBI Taxonomy Ids and returns all the seed scores (competition and cooperation) having ncbiId_A as A
    and ncbiId_B as B but also the other way around.
    It does so for the combinations of all genomes related to each of the 2 NCBI Ids provided.
    Thus, it returns a tuple of 2 lists; list 1 is a list of tuples having genomes of ncbiId_A as A and list 2 another list of
    tuples with genomes of ncbiId_B as A.
    """
    genomes_for_species_A = get_genomes_for_ncbi_tax_id(ncbiId_A)
    genomes_for_species_B = get_genomes_for_ncbi_tax_id(ncbiId_B)
    genomes_for_species_A = [x for x in list(genomes_for_species_A.values())[0] if x.startswith("GCA_") or x.startswith("GCF_")]
    genomes_for_species_B = [x for x in list(genomes_for_species_B.values())[0] if x.startswith("GCA_") or x.startswith("GCF_")]
    seeds_A_B = []
    seeds_B_A = []
    for genome_A in genomes_for_species_A:
        for genome_B in genomes_for_species_B:
            case_seeds_A_B, case_seeds_B_A = get_seed_scores_for_pair_of_genomes(genome_A, genome_B)
            seeds_A_B.append(case_seeds_A_B)
            seeds_B_A.append(case_seeds_B_A)
    return seeds_A_B, seeds_B_A


def get_seed_scores_for_pair_of_genomes(genome_A="GCA_003184265.1", genome_B="GCA_000015645.1"):
    """
    Function to get the seed scores of 2 GC accession ids.
    Returns 2 tuples with the GC ids in the order used (source-target genomes), their competition and their cooperation score.
    """
    patricIds = get_patric_id_of_gc_accession_list([genome_A, genome_B])
    if None in patricIds.values():
        return 0, 0
    query_from_A_to_B, query_from_B_to_A = query_for_getting_seed_scores(patricIds[genome_A], patricIds[genome_B])
    export = execute(query_from_A_to_B)
    try:
        competition_score_A_B, cooperation_score_A_B = export[0][0], export[0][1]
    except ValueError:
        return "Competition and cooperation scores based on the seed sets derived from genome-scale reconstructions of those 2 genomes have not been calculated."
    export = execute(query_from_B_to_A)
    competition_score_B_A, cooperation_score_B_A = export[0][0], export[0][1]
    seeds_A_B = genome_A, genome_B, competition_score_A_B, cooperation_score_A_B
    seeds_B_A = genome_B, genome_A, competition_score_B_A, cooperation_score_B_A
    return seeds_A_B, seeds_B_A


def query_for_getting_seed_scores(patric_genome_A, patric_genome_B):
    """
    Queries for getting from MySQL seed scores of a pair of genomes using both as sources and targets
    Returns 2 strings
    """
    query_from_A_to_B = "".join([
        'SELECT competitionScore, complementaritScore FROM seedScores WHERE patricGenomeA = "',
        patric_genome_A,
        '" AND patricGenomeB = "',
        patric_genome_B,
        '";'
    ])
    query_from_B_to_A = "".join([
        'SELECT competitionScore, complementaritScore FROM seedScores WHERE patricGenomeA = "',
        patric_genome_B,
        '" AND patricGenomeB = "',
        patric_genome_A,
        '";'
    ])
    return query_from_A_to_B, query_from_B_to_A


def get_patric_id_of_gc_accession_list(gc_accesion_list=["GCA_003184265.1"]):
    """
    Gets as input a list of GC accession ids
    Returns a dictionary where the gc ids are the keys and their corresponding PATRIC ids their values
    """
    gc_to_patric_dict = {}
    for gc in gc_accesion_list:
        gc_alt = alt_genome_prefix(gc)
        if gc_alt != 0:
            query = "".join(["SELECT patricId FROM patricId2genomeId where gtdbGenomeAccession = '", gc, "' OR gtdbGenomeAccession = '", gc_alt, "';"])
            patricId = execute(query)
            if len(patricId) > 0 and len(patricId[0]):
                gc_to_patric_dict[gc] = patricId[0][0]
            else:
                gc_to_patric_dict[gc] = None
    return gc_to_patric_dict


def gc_unify(gc_list):
    """
    Removes duplicates of a genome that has entries both as GCA and GCF in the db.
    """
    return list(set([x.replace("GCA", "GCF") for x in gc_list]))


def alt_genome_prefix(gc):
    """
    Switches GCA prefic of a genome accession id to GCF and vice-versa.
    """
    if gc.startswith("GCA_"):
        gc_alt = gc.replace("GCA_", "GCF_")
    elif gc.startswith("GCF_"):
        gc_alt = gc.replace("GCF_", "GCA_")
    else:
        return 0
    return gc_alt


def scores_to_dict(list_of_tuples):
    tmp = {}
    counter = 0
    for case in list_of_tuples:
        if isinstance(case, tuple):
            tmp[counter] = {}
            tmp[counter]["genome_A"] = case[0]
            tmp[counter]["genome_B"] = case[1]
            tmp[counter]["competition"] = case[2]
            tmp[counter]["cooperatiom"] = case[3]
    return tmp
