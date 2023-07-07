from flask import Flask, jsonify, flash, request, session
from werkzeug.utils import secure_filename
import json
import os
import logging
import pandas
import utils
import microbetag.scripts.db_functions as db_functions

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


@app.route('/', methods=['GET', 'POST'])
def index() -> str:
    import random
    message = " ".join(["Hello friend!\nThis is microbetag interface."])
    return message


@app.route('/phen-traits/<string:gtdb_genome_id>', methods=['GET'])
def phen_traits(gtdb_genome_id):
    """ Route to get phenotypic traits for a certain genome id based on phenDB classes """
    return jsonify(db_functions.get_phendb_traits(gtdb_genome_id))


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

        import microbetag.microbetag as microbetag

        if 'file' not in request.files:
            flash('No file part')

        elif 'file' in request.files:
            file = request.files['file']
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            return 'File copied to the server successfully'

        ifile = os.path.join(app.config['UPLOAD_FOLDER'], "users_abundance_table.tsv")

        json_array = request.get_json()
        json_string = json.dumps(json_array)

        df = pandas.read_json(json_string)
        first_col = df.columns[0]
        last_col = df.columns[-1]
        new_first = "seqId"
        new_last = "taxonomy"
        new = {first_col: new_first, last_col: new_last}
        df = df.rename(columns=new)
        df.to_csv(ifile, sep='\t', header=False, index=False)

        # load arguments

        # run microbetag
        microbetag.main()

        """this is until microbetag runs as expected"""
        f = open("/var/tmp/dev_output.json", "r")
        g = json.load(f)

        return g

        """ Now you need to initiate a microbetag run.
        In the end, the file just uploaded will be deleted
        """


if __name__ == '__main__':

    app.run()
