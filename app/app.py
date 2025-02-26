#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The main file of the app
"""

from flask import Flask, flash, url_for, redirect  # imports an instance of Flask
from flask_login import LoginManager
from pony.orm import Database
from pony.flask import Pony
from pony.orm import Required, StrArray, Set
from flask_login import UserMixin
from werkzeug.security import generate_password_hash, check_password_hash


def define_database():
    class User(db.Entity, UserMixin):
        """The User class inherits from db.entity which is used as a base class
        for all entities stored in the database and UserMixin which provides
        default implementations for all of these properties and methods
        This class defines attributes as class variables where the kind
        and type of attribute and attribute options are defined."""
        _table_ = "User"
        username = Required(str, unique=True)
        email = Required(str, unique=True)
        fullname = Required(str)
        password_hash = Required(str)

        """Password hashing is implemented 
        by the two following methods""" # TODO RA: Should be a comment. 

        def check_password(self, password):
            return check_password_hash(self.password_hash, password)

        @classmethod
        def set_password(cls, password):
            return generate_password_hash(password)

    class Calculation(db.Entity):
        _table_ = "Results"
        email = Required(str)
        name = Required(StrArray)
        smiles = Required(StrArray)
        radical = Required(StrArray)
        mode = Required(str)
        date_time = Required(str)
        random_string = Required(str)

    @login.user_loader
    def load_user(user_id):
        """Flask-Login keeps track of the logged in user
        by storing its unique identifier in Flask's user
        session. Each time the logged-in user navigates
        to a new page, Flask-Login retrieves the ID of
        the user from the session, and then loads that
        user into memory by the user loader function"""
        return db.User.get(id=user_id)


pycharm = 0

app = Flask(__name__)  # creates the app object

# setting up the database
login = LoginManager(app)  # initializes the login manager
db = Database()  # Creates an instance of a Pony database
define_database()
PONY_TEST = {
        'provider': 'sqlite',
        'filename': 'test_database.db',
        'create_db': True
    }
db.bind(**PONY_TEST)  # Binds the database using the config app
db.generate_mapping(create_tables=True)
Pony(app)  # links the app to the user database
login.login_view = 'auth.auth'

# Now we need to import the app modules as blueprints
from main import main_bp
app.register_blueprint(main_bp)
from auth import auth_bp
app.register_blueprint(auth_bp)
from info import info_bp
app.register_blueprint(info_bp)
from calculations import calculations_bp
app.register_blueprint(calculations_bp)


# error handler for 500
@app.errorhandler(500)
def internal_error(error):
    flash("Something went wrong, please try again later.")
    return redirect(url_for("main.home"))


# set secret key
app.config['SECRET_KEY'] = 'the random string'

if __name__ == "__main__":
    app.run(host='0.0.0.0', port=80)  # to run the app in PyCharm
