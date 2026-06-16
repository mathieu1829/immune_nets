from flask import Blueprint
from flask import render_template

datasets_bp = Blueprint("datasets", __name__)

@datasets_bp.route("/datasets")
def index():

    datasets = [
        {"id": 1, "name": "Dataset A"},
        {"id": 2, "name": "Dataset B"},
    ]

    return render_template(
        "datasets.html",
        datasets=datasets
    )
