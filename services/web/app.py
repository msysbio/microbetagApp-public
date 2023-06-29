
from typing import List, Dict
from flask import Flask, jsonify, flash, send_from_directory, request, redirect, url_for
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
    genomes_traits = {i: j for i, j in zip(colnames[0], rows)}

    return genomes_traits


def get_complements_for_pair_of_genomes(beneficiarys_genome_id, donors_genome_id):
    """
    For a certain beneficiary genome and a certain donor genome, retrieve all complements
    available in the database
    """
    query = "".join([
            "SELECT uniqueComplements.KoModuleId, uniqueComplements.complement, uniqueComplements.pathway ",
            "FROM uniqueComplements ",
            "INNER JOIN  pathwayComplementarity ",
            "ON pathwayComplementarity.complmentId = uniqueComplements.complementId ",
            "WHERE pathwayComplementarity.beneficiaryGenome = '", beneficiarys_genome_id,
            "' AND pathwayComplementarity.donorGenome = '", donors_genome_id, "';"])

    complements = execute(query)

    return complements


@app.route('/', methods=['GET', 'POST'])
def index() -> str:
    import random
    noise = str(random.randint(1, 100))
    message = " ".join(['Hello friend : ', noise])
    return message


@app.route('/phen-traits/<string:gtdb_genome_id>', methods=['GET'])
def phen_traits(gtdb_genome_id):
    return jsonify(get_phendb_traits(gtdb_genome_id))


@app.route('/complements/<string:beneficiary_genome>/<string:donor_genome>', methods=['GET'])
def complements(beneficiary_genome, donor_genome):
    return jsonify(get_complements_for_pair_of_genomes(beneficiary_genome, donor_genome))


@app.route('/upload-abundance-table', methods=['GET', 'POST'])
def upload():

    if request.method == 'POST':
        if 'file' not in request.files:
            flash('No file part')
            return 'No file added'

        file = request.files['file']
        filename = secure_filename(file.filename)
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

        """ Now you need to initiate a microbetag run. 
        In the end, the file just uploaded will be deleted""" 

        return 'File copied to the server successfully'


if __name__ == '__main__':

    app.run()
