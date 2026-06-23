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
    delete_repertoire,
)

repertoires_bp = Blueprint("repertoires", __name__)

@repertoires_bp.route(
    "/datasets/<uuid:dataset_id>/repertoires"
)
def dataset_view(dataset_id):

    dataset = get_dataset(
        dataset_id
    )

    query = request.args.get("q", "").strip()

    if query:
        dataset.repertoires = [
            rep for rep in dataset.repertoires
            if query.lower() in rep.name.lower()
        ]

    return render_template(
        "dataset_detail.html",
        dataset=dataset
    )

@repertoires_bp.route(
    "/datasets/<uuid:dataset_id>/reprtoires/<uuid:repertoire_id>/delete",
    methods=["POST"]
)
def delete(dataset_id, repertoire_id):

    success = delete_repertoire(repertoire_id)
    
    if not success:

        flash(
            "Dataset not found.",
            "error"
        )

        return redirect(
            url_for("repertories.index")
        )

    flash(
        "Dataset deleted.",
        "success"
    )

    return redirect(
        url_for("repertoires.dataset_view",
                dataset_id=dataset_id)
        )
