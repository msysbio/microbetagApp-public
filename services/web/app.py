
from typing import List, Dict
from flask import Flask, jsonify, send_from_directory
from flask_sqlalchemy import SQLAlchemy
import mysql.connector
import json

app = Flask(__name__)

app.config.from_object("config.Config")
db = SQLAlchemy(app)

# Database connection configuration
config = {
    'user': 'root',
    'password': 'pass',
    'host': 'db',  # Replace with the actual database service name or IP address
    'port': 3306,  # Replace with the actual port
    'database': 'microbetagDB'
}


class User(db.Model):
    __tablename__ = "users"

    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(128), unique=True, nullable=False)
    active = db.Column(db.Boolean(), default=True, nullable=False)

    def __init__(self, email):
        self.email = email


def favorite_colors() -> List[Dict]:

    # Establish a database connection
    cnx = mysql.connector.connect(**config)
    cursor = cnx.cursor()

    query = "SELECT * FROM test"
    cursor.execute(query)

    results = []
    for row in cursor:
        results.append(row)

    return results


@app.route('/')
def index() -> str:
    return json.dumps({'favorite_colors': favorite_colors()})


if __name__ == '__main__':
    import os
    # print(os.environ)

    app.run()
    # app.run()
