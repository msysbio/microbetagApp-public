from flask import Flask, jsonify, flash, request, session
from flask_caching import Cache
from werkzeug.utils import secure_filename
import json
import os
import logging
import pandas as pd
import microbetag.microbetag as microbetag
import microbetag.scripts.db_functions as db_functions
from datetime import datetime
import time

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    handlers=[
                        logging.FileHandler('app.log'),
                        logging.StreamHandler()
                    ]
                    )
# Interface for user-server
UPLOAD_FOLDER = '/tmp'

app = Flask(__name__)
app.debug = True
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 16 * 1000 * 1000
app.secret_key = "\xf0+\x91\xe0c\xb8ov\xf3c\xc7b#mS\xbc\xf1\xdd\xb2\x06\x8an*\xc0"

cache = Cache(app, config={"CACHE_TYPE": "simple", "CACHE_DEFAULT_TIMEOUT": 20})


@app.route('/', methods=['GET', 'POST'])
def index() -> str:
    import random
    message = " ".join(["Hello friend!\nThis is microbetag interface."])
    return message


@app.route('/phen-traits/<string:gtdb_genome_id>', methods=['GET'])
def phen_traits(gtdb_genome_id):
    """ Route to get phenotypic traits for a certain genome id based on phenDB classes """
    phen_traits = db_functions.get_phendb_traits(gtdb_genome_id)
    if phen_traits == 0:
        message = """No Phen traits for the genome id asked. Try replacing GCA with GCF or vice versa.\
        If still no hits, then there is no relative info currently on microbetag DB."""
        return message
    else:
        return jsonify(phen_traits)


@app.route('/ncbiTaxId-to-genomeId/<int:ncbiTaxId>', methods=['GET'])
def get_related_genomes(ncbiTaxId):
    return jsonify(db_functions.get_genomes_for_ncbi_tax_id(ncbiTaxId))


@app.route('/genome-complements/<string:beneficiary_genome>/<string:donor_genome>', methods=['GET'])
def genomes_complements(beneficiary_genome, donor_genome):
    return jsonify(db_functions.get_complements_for_pair_of_genomes(beneficiary_genome, donor_genome))


@app.route('/complements/<int:beneficiary_ncbi_tax_id>/<int:donor_ncbi_tax_id>', methods=['GET'])
def ncbi_ids_complements(beneficiary_ncbi_tax_id, donor_ncbi_tax_id):
    return jsonify(db_functions.get_complements_for_pair_of_ncbiIds(beneficiary_ncbi_tax_id, donor_ncbi_tax_id))


@app.route('/upload-abundance-table', methods=['GET', 'POST'])
def upload():
    if request.method == 'POST':

        if 'file' in request.files:
            tfile = request.files['file']
            filename = secure_filename(tfile.filename)
            tfile.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            return 'File copied to the server successfully'

        f = open("/var/tmp/dev_output.json", "r")
        g = json.load(f)

        return g


@app.route('/upload-abundance-table-dev', methods=['GET', 'POST'])
def upload_dev():

    if request.method == 'POST':

        # Generate a unique identifier for the request based on its content; an integer
        request_id = hash(frozenset(request.form.items()))

        # Check if the request has already been processed recently (within 20 seconds; set in the cache above)
        if cache.get(request_id):
            # Return the cached response for the duplicate request
            response = cache.get(request_id)
            return jsonify({"status": "success", "data": response})

        # Create output dir for the session
        now = datetime.now()
        dt_string = now.strftime("%d-%m-%Y--%H.%M.%S")
        out_dir = os.path.join(app.config['UPLOAD_FOLDER'], dt_string)
        os.mkdir(out_dir)

        # Parse user's input
        json_array = request.get_json()

        data = json_array["data"]
        data = pd.DataFrame(data)  # pd.read_json(data)

        metadata = json_array["metadata"]

        data.iloc[0, 0] = "seqId"
        data.iloc[0, -1] = "taxonomy"

        # CytoApp will provide a list of args while API is better to give a dictionary in terms of user friendliness
        if isinstance(json_array["inputParameters"], list):
            arguments_list = ["input_category", "taxonomy", "delimiter", "sensitive", "heterogeneous", 
                "phenDB", "faprotax", "pathway_complement", "seed_scores", "manta"]
            args = {}
            for i, j in zip(arguments_list, json_array["inputParameters"]):
                if j == "true":
                    j = True
                args[i] = j
        else:
            args = json_array["inputParameters"]

        # Write the user's abundance table as a .tsv file
        abundance_table = os.path.join(out_dir, "abundance_table.tsv")
        data.to_csv(abundance_table, sep='\t', header=False, index=False)

        if metadata:
            # [TODO]: check if an empty metadata entry is a nonetype or a string of 0 length 
            metadata_table = os.path.join(out_dir, "metadata.tsv")
            metadata.to_csv(metadata_table, sep="\t", header=False, index=False)

        # run microbetag
        annotated_network = microbetag.main(out_dir, args, abundance_table)

        # Cache the response and set the request_id as the cache key
        # cache.set(request_id, annotated_network)


        # Zip file Initialization and you can change the compression type
        zipfolder = zipfile.ZipFile('microbetag.zip','w', compression = zipfile.ZIP_STORED)

        # zip all the files which are inside in the folder
        for root,dirs, files in os.walk(out_dir):
            for file in files:
                zipfolder.write(os.path.join(out_dir,file))
        zipfolder.close()

        return send_file('microbetag.zip',
                mimetype = 'zip',
                attachment_filename= 'microbetag.zip',
                as_attachment = True)

        # Delete the zip file if not needed
        os.remove("microbetag.zip")



        return annotated_network


if __name__ == '__main__':

    app.run()
