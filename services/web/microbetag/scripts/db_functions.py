import sys
import os
import decimal
import pandas as pd
import mysql.connector
from .variables import *


def check_connection_to_db(db_user, db_password, db_host, db_database):
    """
    Set a connection to the database
    """
    try:
        db = mysql.connector.connect(
            user=db_user,
            password=db_password,
            host=db_host,
            database=db_database)
        cursor = db.cursor()
        cursor.execute("SELECT VERSION()")
        results = cursor.fetchone()

        if results:
            return True
        else:
            return False

    except mysql.connector.Error:
        return False


def query_to_microbetagDB(phrase):
    """
    Function to get functional traits as assigned in a genome from the phenotrex software
    phenotrex.readthedocs.io/ and stored in the microbetagDB.
    """
    try:
        cnx = mysql.connector.connect(
            user=USER_NAME,
            password=PASSWORD,
            host=HOST,
            database=DB_NAME)
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


def get_phenDB_traits(genome_id):
    """
    Query to get the phenotypic traits for a genome id
    """
    phrase = "SELECT * FROM phenDB WHERE gtdbId = '" + genome_id + "';"
    phendb_traits = list(query_to_microbetagDB(phrase)[0])
    phendb_traits = [
        str(x) if isinstance(
            x, decimal.Decimal) else x for x in phendb_traits]
    return phendb_traits


def get_column_names():
    """
    Get the column names of a database table
    """
    phrase = "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='" +\
             DB_NAME + "' AND `TABLE_NAME`='phenDB';"
    colnames = query_to_microbetagDB(phrase)
    colnames = [x[0] for x in colnames]
    return colnames


# HERE IS THE UTILS
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


def get_phendb_traits(gtdb_genome_id="GCA_018819265.1"):
    """
    Get phenotypical traits based on phenDB classes based on its GTDB representative genome
    """
    query = "".join(
        ["SELECT * FROM phenDB WHERE gtdbId = '", gtdb_genome_id, "';"])
    rows = execute(query)

    query_colnames = "SHOW COLUMNS FROM phenDB;"
    colnames = [list(x)[0] for x in execute(query_colnames)]
    genomes_traits = {i: j for i, j in zip(colnames, rows[0])}

    return genomes_traits


def get_genomes_for_ncbi_tax_id(ncbi_tax_id=1281578):
    """
    Get the genomes IDs that correspond to a NCBI Taxonomy Id and are present in the microbetagDB
    """
    query = "".join(
        ["SELECT genomeId from genome2taxNcbiId where ncbiTaxId = ", str(ncbi_tax_id), ";"])
    genome_ids = execute(query)

    return {ncbi_tax_id: [x[0] for x in genome_ids]}


def get_complements_for_pair_of_genomes(
        beneficiarys_genome_id="GCA_003184265.1", donors_genome_id="GCA_000015645.1"):
    """
    For a certain beneficiary genome and a certain donor genome, retrieve all complements
    available in the database.
    """
    complements_ids_list_query = "".join(['SELECT complmentId FROM pathwayComplementarity WHERE beneficiaryGenome = "',
                                          beneficiarys_genome_id,
                                          '" AND donorGenome = "',
                                          donors_genome_id,
                                          '";'
                                          ])
    export = execute(complements_ids_list_query)
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
            # tmp = {"beneficiary_genome": beneficiary_genome, "donor_genome": donor_genome, "potential complementarities": get_complements_for_pair_of_genomes(beneficiary_genome, donor_genome)}

    return complements


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
