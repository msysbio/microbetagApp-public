
from typing import List, Dict
from flask import Flask, jsonify, flash, request, session
from werkzeug.utils import secure_filename
import mysql.connector
import json
import os

app = Flask(__name__)
UPLOAD_FOLDER = '/tmp'
ALLOWED_EXTENSIONS = {'txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif'}

app = Flask(__name__)
app.debug = True
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 16 * 1000 * 1000
app.secret_key = "\xf0+\x91\xe0c\xb8ov\xf3c\xc7b#mS\xbc\xf1\xdd\xb2\x06\x8an*\xc0"

# Database connection configuration
config = {
    'user': 'msysbio',
    'password': 'pass',
    'host': 'db',  # Replace with the actual database service name or IP address
    'port': 3306,  # Replace with the actual port
    'database': 'microbetagDB'
}


def execute(phrase):
    """
    Establish a database connection and perform an action
    """
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
    query = "".join(["SELECT * FROM phenDB WHERE gtdbId = '", gtdb_genome_id, "';"])
    rows = execute(query)

    query_colnames = "SHOW COLUMNS FROM phenDB;"
    colnames = [list(x)[0] for x in execute(query_colnames)]
    genomes_traits = {i: j for i, j in zip(colnames, rows[0])}

    return genomes_traits


def get_genomes_for_ncbi_tax_id(ncbi_tax_id=1281578):
    """
    Get the genomes IDs that correspond to a NCBI Taxonomy Id and are present in the microbetagDB
    """
    query = "".join(["SELECT genomeId from genome2taxNcbiId where ncbiTaxId = ", str(ncbi_tax_id), ";"])
    genome_ids = execute(query)
    return {ncbi_tax_id: [x[0] for x in genome_ids]}


def get_complements_for_pair_of_genomes(beneficiarys_genome_id="GCA_003184265.1", donors_genome_id="GCA_000015645.1"):
    """
    For a certain beneficiary genome and a certain donor genome, retrieve all complements
    available in the database
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
        query = "".join(["SELECT KoModuleId, complement, pathway FROM uniqueComplements WHERE complementId = '", complementId, "';"])
        compl = execute(query)
        complements.append(compl)
    return complements


def get_complements_for_pair_of_ncbiIds(beneficiarys_ndbi_id=1281578, donors_ncbi_id=146891):
    """Gets a pair of NCBI Taxonomy Ids and returns and all the potential complements,
    based on all the corresponding genomes.
    """
    beneficiarys_genomes = get_genomes_for_ncbi_tax_id(beneficiarys_ndbi_id)
    donors_genomes = get_genomes_for_ncbi_tax_id(donors_ncbi_id)
    complements = {}
    for beneficiary_genome in list(beneficiarys_genomes.values())[0]:
        complements[beneficiary_genome] = {}
        for donor_genome in list(donors_genomes.values())[0]:
            complements[beneficiary_genome][donor_genome] = get_complements_for_pair_of_genomes(beneficiary_genome, donor_genome)
            # tmp = {"beneficiary_genome": beneficiary_genome, "donor_genome": donor_genome, "potential complementarities": get_complements_for_pair_of_genomes(beneficiary_genome, donor_genome)}
    return complements


@app.route('/', methods=['GET', 'POST'])
def index() -> str:
    import random
    message = " ".join(["Hello friend!\nThis is microbetag interface."])
    return message


@app.route('/phen-traits/<string:gtdb_genome_id>', methods=['GET'])
def phen_traits(gtdb_genome_id):
    """ Route to get phenotypic traits for a certain genome id based on phenDB classes """
    return jsonify(get_phendb_traits(gtdb_genome_id))


@app.route('/ncbiTaxId-to-genomeId/<int:ncbiTaxId>', methods=['GET'])
def get_related_genomes(ncbiTaxId):
    return jsonify(get_genomes_for_ncbi_tax_id(ncbiTaxId))


@app.route('/genome-complements/<string:beneficiary_genome>/<string:donor_genome>', methods=['GET'])
def genomes_complements(beneficiary_genome, donor_genome):
    return jsonify(get_complements_for_pair_of_genomes(beneficiary_genome, donor_genome))


@app.route('/complements/<int:beneficiary_ncbi_tax_id>/<int:donor_ncbi_tax_id>', methods=['GET'])
def ncbi_ids_complements(beneficiary_ncbi_tax_id, donor_ncbi_tax_id):
    return jsonify(get_complements_for_pair_of_ncbiIds(beneficiary_ncbi_tax_id, donor_ncbi_tax_id))


@app.route('/upload-abundance-table', methods=['GET', 'POST'])
def upload():

    if request.method == 'POST':
        if 'file' not in request.files:
            flash('No file part')

        elif 'file' in request.files:
            file = request.files['file']
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            return 'File copied to the server successfully'

        json_array = request.get_json()
        ifile = os.path.join(app.config['UPLOAD_FOLDER'], "test_from_get_json.tsv")
        json_string = json.dumps(json_array)
        with open(ifile, "w") as outfile:
            outfile.write(json_string)
            """this is until microbetag runs as expected"""
            f = open("/var/tmp/dev_output.json", "r")
            g = json.load(f)
            return g

        """ Now you need to initiate a microbetag run.
        In the end, the file just uploaded will be deleted
        """


if __name__ == '__main__':

    app.run()
