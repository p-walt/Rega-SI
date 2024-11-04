from flask import Blueprint

info_bp = Blueprint('info', __name__)

from info import routes