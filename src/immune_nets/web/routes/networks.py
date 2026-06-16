from flask import Blueprint
from flask import render_template
from flask import request

networks_bp = Blueprint("networks", __name__)

@networks_bp.route("/networks/build", methods=["GET", "POST"])
def index():

    datasets = [
        {"id": 1, "name": "Dataset A"},
        {"id": 2, "name": "Dataset B"},
    ]

    if request.method == "POST":

        dataset_id = request.form["dataset_id"]
        threshold = float(
            request.form["threshold"]
        )

        # network_service.build(...)

    return render_template(
        "networks.html",
        datasets=datasets
    )
