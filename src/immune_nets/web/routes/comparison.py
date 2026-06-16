from flask import Blueprint
from flask import render_template
from flask import request

comparison_bp = Blueprint("comparison", __name__)

@comparison_bp.route("/compare", methods=["GET", "POST"])
def index():

    networks = [
        {"id": 1, "name": "Network A"},
        {"id": 2, "name": "Network B"},
    ]

    result = None

    if request.method == "POST":

        net_a = request.form["network_a"]
        net_b = request.form["network_b"]

        # result = comparison_service.compare(...)

        result = 0.87

    return render_template(
        "comparison.html",
        networks=networks,
        result=result
    )
