from flask import (
    Blueprint,
    render_template,
    request,
    redirect,
    url_for,
    flash
)

from immune_nets.services.dataset_service import (
    get_dataset,
)

from immune_nets.services.repertoire_service import (
    get_repertoire,
)

clonotypes_bp = Blueprint("clonotypes", __name__)

@clonotypes_bp.route(
    "/datasets/<uuid:dataset_id>/repertoires/<uuid:repertoire_id>/clonotypes"
)
def repertoire_view(dataset_id, repertoire_id):

    dataset = get_dataset(
        dataset_id
    )

    repertoire = get_repertoire(
        repertoire_id
    )

    query = request.args.get("q", "").strip()

    if query:
        repertoire.clonotypes = [
            clone for clone in repertoire.clonotypes
            if query.lower() in str(clone.tcrb_aa).lower()
        ]

    return render_template(
        "clonotypes.html",
        dataset=dataset,
        repertoire=repertoire
    )


