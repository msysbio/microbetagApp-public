
from typing import List, Dict
from flask import Flask, jsonify, send_from_directory, request
import mysql.connector
import json

app = Flask(__name__)

# Database connection configuration
config = {
    'user': 'msysbio',
    'password': 'pass',
    'host': 'db',  # Replace with the actual database service name or IP address
    'port': 3306,  # Replace with the actual port
    'database': 'microbetagDB'
}


def get_phendb_traits(gtdb_genome_id="GCA_018819265.1") -> List[Dict]:

    # Establish a database connection
    cnx = mysql.connector.connect(**config)
    cursor = cnx.cursor()

    query = "".join(["SELECT * FROM phenDB WHERE gtdbId = '", gtdb_genome_id, "';"])
    cursor.execute(query)
    rows = cursor.fetchall()

    return rows


@app.route('/')
def index() -> str:
    return 'Hello friend'


@app.route('/phen-traits-db/<string:gtdb_genome_id>', methods=['GET'])
def phen_traits(gtdb_genome_id):
    return jsonify(get_phendb_traits(gtdb_genome_id))


@app.route('/upload-abundance-table')
def upload():

    print("hello friend")

    if 'file' not in request.files:
        return 'No file uploaded', 400

    file = request.files['file']
    file.save(file.filename)

    # process_file(file.filename)

    return "hello friend"
    # return 'File uploaded and processed successfully', 200


if __name__ == '__main__':

    import os

    app.run()
