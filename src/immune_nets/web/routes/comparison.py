from flask import Blueprint
from flask import render_template
from flask import request


from immune_nets.services.networks_service import get_all_networks, delete_network, add_network, get_network
from immune_nets.mappers import NetworkMapper
from immune_nets.entities import ImmuneNetwork, GraphStats
from immune_nets.analysis.statDistances import PairwiseDistributionDistance

comparison_bp = Blueprint("comparison", __name__)

@comparison_bp.route("/compare", methods=["GET", "POST"])
def index():

    result = None

    networks = get_all_networks() 
    if request.method == "POST":

        net_a = request.form["network_a"]
        net_b = request.form["network_b"]

        network_a = get_network(net_a)
        network_b = get_network(net_b)

        immuneNetwork_a = NetworkMapper.toImmuneNetwork(network_a)
        immuneNetwork_b = NetworkMapper.toImmuneNetwork(network_b)

        stats_a = GraphStats(immuneNetwork_a)
        stats_b = GraphStats(immuneNetwork_b)

        distrubtionDistance = PairwiseDistributionDistance("degreeDistribution")



        result = distrubtionDistance.stat_dist(stats_a, stats_b)


    return render_template(
        "compare.html",
        networks=networks,
        result=result
    )
