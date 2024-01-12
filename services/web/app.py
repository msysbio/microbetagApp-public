from flask import Flask, jsonify, flash, request, session, send_file
from flask_caching import Cache
# from werkzeug.utils import secure_filename
import json
import os
import logging
import pandas as pd
import microbetag.microbetag as microbetag
import microbetag.scripts.db_functions as db_functions
from datetime import datetime
import zipfile

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
def phen_traits(gtdb_genome_id="GCF_002550015.1"):
    """ Route to get phenotypic traits for a certain genome id based on phenDB classes """
    phen_traits = db_functions.get_phendb_traits(gtdb_genome_id)
    if phen_traits == 0:
        message = """No Phen traits for the genome id asked. \
        Make sure you are asking for a GTDB v202 representative genome."""
        return message
    else:
        return jsonify(phen_traits)
        # return phen_traits


@app.route('/ncbiTaxId-to-genomeId/<int:ncbiTaxId>', methods=['GET'])
def get_related_genomes(ncbiTaxId=1598):
    q = db_functions.get_genomes_for_ncbi_tax_id(ncbiTaxId)
    return jsonify(q)
    # return q


@app.route('/genome-complements/<string:beneficiary_genome>/<string:donor_genome>', methods=['GET'])
def genomes_complements(beneficiary_genome="GCA_011364525.1", donor_genome="GCA_002980625.1"):
    """
    Gets complements that a donor can provide to a baneficiary genome.

    """
    r = db_functions.get_complements_for_pair_of_genomes(beneficiary_genome, donor_genome)
    d = {}
    d["beneficiary-genome"] = beneficiary_genome
    d["donor-genome"] = donor_genome
    d["complements"] = {}
    counter = 0
    for case in r:
        d["complements"][counter] = {}
        d["complements"][counter]["module"] = case[0][0]
        d["complements"][counter]["complement"] = case[0][1]
        d["complements"][counter]["complete-alternative"] = case[0][2]
        d["complements"][counter]["coloured-map"] = case[0][3]
        counter += 1
    # return d
    return jsonify(d)


@app.route('/complements/<string:beneficiary_ncbi_tax_id>/<string:donor_ncbi_tax_id>', methods=['GET'])
def ncbi_ids_complements(beneficiary_ncbi_tax_id="2184738", donor_ncbi_tax_id="86180"):
    q = db_functions.get_complements_for_pair_of_ncbiIds(beneficiary_ncbi_tax_id, donor_ncbi_tax_id)
    # return q
    return jsonify(q)


@app.route('/seed-complements/<string:beneficiaryId>/<string:donorId>/<string:typeCategory>', methods=['GET'])
def seed_complements(beneficiaryId="2184738", donorId="86180", typeCategory="ncbiTaxonomyIds"):
    q = db_functions.get_paired_seed_complements(beneficiary=beneficiaryId, donor=donorId, type=typeCategory)
    # return jsonify(Q)
    return q

@app.route('/seed-scores/<string:ncbiId_A>/<string:ncbiId_B>', methods=['GET'])
def seed_scores(ncbiId_A="86180", ncbiId_B="2184738"):
    """
    Gets genomes related to the ncbi ids provided and returns their seed scores using the all the genomes retrieved as
    the target seed set.
    """
    r = db_functions.get_seed_scores_for_pair_of_ncbiIds(ncbiId_A, ncbiId_B)
    a2b, b2a = r[0], r[1]
    d = {}
    d[0] = {}
    d[1] = {}
    d[0]["A"] = ncbiId_A
    d[0]["B"] = ncbiId_B
    d[0]["scores"] = db_functions.scores_to_dict(a2b)
    d[1]["A"] = ncbiId_B
    d[1]["B"] = ncbiId_A
    d[1]["scores"] = db_functions.scores_to_dict(b2a)
    # return d
    return jsonify(d)


@app.route('/genomes-seed-scores/<string:gc_seed_set_in>/<string:target_gc>', methods=['GET'])
def genomes_seed_scores(gc_seed_set_in="GCA_011364525.1", target_gc="GCA_002980625.1"):
    """
    Returns seed scores using both gc ids as the target seed set.
    """
    p = db_functions.get_seed_scores_for_pair_of_genomes(gc_seed_set_in, target_gc)
    # return p
    return jsonify(p)


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

        if data.shape[0] > 1000 and "network" not in json_array:
            return ValueError("""
                              Your abundance table consists of more than 1000 sequences. 
                              Please build a co-occurrence network either using the microbetag_prep Docker image or in any way of your choice
                              and submit again your job providing it too.
                              """)

        data.iloc[0, 0] = "seqId"
        data.iloc[0, -1] = "taxonomy"

        # CytoApp will provide a list of args while API is better to give a dictionary in terms of user friendliness
        arguments_list = [
            "input_category", "taxonomy", "delimiter",
            "phenDB", "faprotax", "pathway_complement", "seed_scores",
            "get_children", "manta", "sensitive", "heterogeneous",
        ]
        if isinstance(json_array["inputParameters"], list):
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

        if "metadata" in json_array:
            # [TODO]: check if an empty metadata entry is a nonetype or a string of 0 length
            metadata_table = os.path.join(out_dir, "metadata.tsv")
            metadata_table.to_csv(metadata_table, sep="\t", header=False, index=False)

        # run microbetag
        logging.info(json_array["inputParameters"])

        annotated_network = microbetag.main(out_dir, args, abundance_table)

        # # Zip file Initialization and you can change the compression type
        # zipfolder = zipfile.ZipFile('microbetag.zip', 'w', compression=zipfile.ZIP_STORED)

        # # zip all the files which are inside in the folder
        # for root, dirs, files in os.walk(out_dir):
        #     for file in files:
        #         zipfolder.write(os.path.join(out_dir, file))
        # zipfolder.close()

        # return send_file(
        #     'microbetag.zip',
        #     mimetype='zip',
        #     attachment_filename='microbetag.zip',
        #     as_attachment=True
        # )

        # # Delete the zip file if not needed
        # os.remove("microbetag.zip")

        return annotated_network


@app.route('/upload-abundance-table', methods=['GET', 'POST'])
def test_output_dev():

    if request.method == 'POST':

        # Generate a unique identifier for the request based on its content; an integer
        request_id = hash(frozenset(request.form.items()))

        # Check if the request has already been processed recently (within 20 seconds; set in the cache above)
        if cache.get(request_id):
            # Return the cached response for the duplicate request
            response = cache.get(request_id)
            return jsonify({"status": "success", "data": response})

    d = json.load(open("static/eth_test2.cx", "r"))
    return d


if __name__ == '__main__':

    app.run()
