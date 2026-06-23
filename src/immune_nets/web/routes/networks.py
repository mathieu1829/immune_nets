from flask import (
    Blueprint,
    render_template,
    request,
    redirect,
    url_for,
    flash,
    Response
)

import io
from matplotlib import pyplot as plt

from immune_nets.services.networks_service import get_all_networks, delete_network, add_network, get_network
from immune_nets.services.repertoire_service import get_all_repertoires, get_repertoire
from immune_nets.mappers import RepertoireMapper, NetworkMapper
from immune_nets.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from immune_nets.creation.distance.alignment import sequenceAligner
from immune_nets.analysis.visualization.graphVisualization import graphVisualization

import matplotlib
matplotlib.use("Agg")

networks_bp = Blueprint("networks", __name__)

@networks_bp.route("/networks", methods=["GET", "POST"])
def index():

    networks = get_all_networks() 
    repertoires = get_all_repertoires()

    query = request.args.get("q")

    if query:
        networks = [
            n for n in networks
            if query.lower() in str(n.network_id).lower()
        ]

    if request.method == "POST":
        repertoire_id = request.form.get(
            "repertoire_id"
        )
        if not repertoire_id:
            flash(
                "Please select a repertoire.",
                "error"
            )
            return render_template(
                "networks.html",
                networks=networks,
                repertoires=repertoires
            )
        repertoire = get_repertoire(repertoire_id)
        immuneRepertoire = RepertoireMapper.toImmuneRepertoire(repertoire)
        
        distanceFun = sequenceAligner("BLOSUM62")
        algorithm = simpleVectorBetaDistance
        threshold = 0.22

        immuneNetwork = algorithm(repertoire=immuneRepertoire, 
                                  distance=distanceFun, 
                                  threshold=threshold)
        network_id = add_network(immuneNetwork)
        
        flash(
            "Network created.",
            "success"
        )
        return redirect(
            url_for("networks.view_network", network_id=network_id)
        )

        

    return render_template(
        "networks.html",
        networks=networks,
        repertoires=repertoires
    )

@networks_bp.route(
    "/networks/<uuid:network_id>/delete",
    methods=["POST"]
)
def delete(network_id):

    success = delete_network(network_id)
    
    if not success:

        flash(
            "Network not found.",
            "error"
        )

        return redirect(
            url_for("datasets.index")
        )

    flash(
        "Network deleted.",
        "success"
    )

    return redirect(
        url_for("networks.index")
        )

@networks_bp.route(
    "/networks/<uuid:network_id>/image",
    methods=["GET"]
)
def network_image(network_id):
    network = get_network(network_id)

    if not network:
        return ("Not found", 404)

    immuneNetwork = NetworkMapper.toImmuneNetwork(network)

    fig, ax = plt.subplots(figsize=(4,8))

    ax.set_title("")
    graphVisualization(immuneNetwork, ax)

    buf = io.BytesIO()

    fig.savefig(
        buf,
        format="png",
        bbox_inches="tight"
    )

    plt.close(fig)

    buf.seek(0)

    return Response(
        buf.read(),
        mimetype="image/png",
        headers={
            "Cache-Control": "no-cache"
        }
    )

@networks_bp.route(
    "/networks/<uuid:network_id>",
    methods=["GET"]
)
def view_network(network_id):

    network = get_network(network_id)

    if not network:
        flash("Network not found.", "error")
        return redirect(url_for("networks.index"))

    networks = get_all_networks()
    repertoires = get_all_repertoires()

    return render_template(
        "networks.html",
        networks=networks,
        repertoires=repertoires,
        selected_network=network
    )

    






