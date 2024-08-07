from flask import Flask, jsonify, flash, request, session, send_file
from flask_caching import Cache
import json
import os
from microbetag.scripts.logger import setup_logger
import pandas as pd
import microbetag.microbetag as microbetag
import microbetag.scripts.db_functions as db_functions
from datetime import datetime


__version__ = "1.0.1"


# Interface for user-server
UPLOAD_FOLDER = '/tmp'

app = Flask(__name__)
app.debug = True
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 16 * 1000 * 1000
app.secret_key = "\xf0+\x91\xe0c\xb8ov\xf3c\xc7b#mS\xbc\xf1\xdd\xb2\x06\x8an*\xc0"

cache = Cache(app, config={"CACHE_TYPE": "simple", "CACHE_DEFAULT_TIMEOUT": 20})

# Create the log file path
log_dir = "/tmp"
logfile = os.path.join(log_dir, "log.txt")

# Setup logger using the function from logger.py
logger = setup_logger(logfile)


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


@app.route('/ncbiTaxId-to-genomeId/<int:ncbiTaxId>', methods=['GET'])
def get_related_genomes(ncbiTaxId=1598):
    q = db_functions.get_genomes_for_ncbi_tax_id(ncbiTaxId)
    return jsonify(q)


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
    return jsonify(d)


@app.route('/complements/<string:beneficiary_ncbi_tax_id>/<string:donor_ncbi_tax_id>', methods=['GET'])
def ncbi_ids_complements(beneficiary_ncbi_tax_id="2184738", donor_ncbi_tax_id="86180"):
    q = db_functions.get_complements_for_pair_of_ncbiIds(beneficiary_ncbi_tax_id, donor_ncbi_tax_id)
    return jsonify(q)


@app.route('/seed-complements/<string:beneficiaryId>/<string:donorId>/<string:typeCategory>', methods=['GET'])
def seed_complements(beneficiaryId="2184738", donorId="86180", typeCategory="ncbiTaxonomyIds"):
    q = db_functions.get_paired_seed_complements(beneficiary=beneficiaryId, donor=donorId, type=typeCategory)
    return jsonify(q)


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
        users_input = request.get_json()
        logger.info("MGG sends following params: %s", users_input["inputParameters"])

        data = users_input["data"]
        data = pd.DataFrame(data)  # pd.read_json(data)

        # CytoApp will provide a list of args while API is better to give a dictionary in terms of user friendliness
        if isinstance(users_input["inputParameters"], list):
            # If user does not set an argument MGG will set it as false
            args = {}
            for arg in users_input["inputParameters"]:
                argument, value = arg.split(":")
                logger.info(argument)
                logger.info(value)
                if value == "true":
                    evalue = True
                elif value == "false":
                    evalue = False
                else:
                    evalue = value
                args[argument] = evalue
        else:
            args = users_input["inputParameters"]

        logger.info("haris>", args)

        # Check whether the analysis is possible through the server
        if data.shape[0] > 1000:
            accepted_tax_schemes = ["GTDB", "Silva", "microbetag_prep"]
            if "network" not in users_input or args["taxonomy"] not in accepted_tax_schemes:
                return ValueError("""
                                Your abundance table consists of more than 1000 sequences.
                                Please build a co-occurrence network either using the microbetag_prep Docker image or in any way of your choice
                                and submit again your job providing it too.
                                """)

        if "network" in users_input:
            # Skip header line
            network = pd.DataFrame(users_input["network"][1:])
            network_path = os.path.join(out_dir, "edgelist.tsv")
            network.to_csv(network_path, sep='\t', header=False, index=False)
            edge_list_file = network_path
        else:
            edge_list_file = None

        if "metadata" in users_input:
            metadata = pd.DataFrame(users_input["metadata"])
            metadata_path = os.path.join(out_dir, "metadata.tsv")
            metadata.to_csv(metadata_path, sep='\t', header=False, index=False)
            metadata_file = metadata_path
        else:
            metadata_file = None

        # Write the user's abundance table as a .tsv file
        abundance_table = os.path.join(out_dir, "abundance_table.tsv")
        data.iloc[0, 0] = "seqId"
        data.iloc[0, -1] = "taxonomy"
        data.to_csv(abundance_table, sep='\t', header=False, index=False)

        # run microbetag
        try:
            annotated_network = microbetag.main(out_dir=out_dir,
                                                cfg=args,
                                                abundance_table_file=abundance_table,
                                                edge_list_file=edge_list_file,
                                                metadata_file=metadata_file,
                                                logger=logger)
            return annotated_network
        except Exception as e:
            logger.error(e)
            return jsonify({'status': 'error', 'message': str(e)})


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
