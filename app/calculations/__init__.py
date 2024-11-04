from flask import Blueprint

calculations_bp = Blueprint('calculations', __name__)

from calculations import routes