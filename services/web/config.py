import os


basedir = os.path.abspath(os.path.dirname(__file__))


class Config(object):
    'username:password@localhost/db_name'
    # Based on https://docs.sqlalchemy.org/en/20/dialects/mysql.html#module-sqlalchemy.dialects.mysql.mysqlconnector
    # mysql+mysqlconnector://<user>:<password>@<host>[:<port>]/<dbname>

    SQLALCHEMY_DATABASE_URI = os.getenv("DATABASE_URL", "mysql+mysqlconnector://")
    SQLALCHEMY_TRACK_MODIFICATIONS = False
