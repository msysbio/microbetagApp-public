
from typing import List, Dict
from flask import Flask, jsonify, send_from_directory
from flask_sqlalchemy import SQLAlchemy
import mysql.connector
import json

app = Flask(__name__)

app.config.from_object("config.Config")
db = SQLAlchemy(app)


# cnx = mysql.connector.connect(user=os.getenv("MYSQL_USER"), password=os.getenv("MYSQL_PASSWORD"), host='localhost', database=os.getenv("MYSQL_DATABASE"))


class User(db.Model):
    __tablename__ = "users"

    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(128), unique=True, nullable=False)
    active = db.Column(db.Boolean(), default=True, nullable=False)

    def __init__(self, email):
        self.email = email


def favorite_colors() -> List[Dict]:
    # config = {
    #     'user': 'root',
    #     'password': 'root',
    #     'host': 'db',
    #     'port': '3306',
    #     'database': 'knights'
    # }
    # connection = mysql.connector.connect(**config)
    # cursor = connection.cursor()
    # cursor.execute('SELECT * FROM favorite_colors')
    # results = [{name: color} for (name, color) in cursor]
    # cursor.close()
    # connection.close()

    results = [1, 2, 3, 4, 5, 5, 6, 67]

    return results


@app.route('/')
def index() -> str:
    return json.dumps({'favorite_colors': favorite_colors()})


if __name__ == '__main__':
    import os
    # print(os.environ)

    app.run()
    # app.run()
