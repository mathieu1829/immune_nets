from flask import (
    Blueprint,
    render_template,
    request,
    redirect,
    url_for,
    flash
)

from immune_nets.services.dataset_service import (
    get_all_datasets,
    create_dataset,
    delete_dataset
)

datasets_bp = Blueprint("datasets", __name__)

@datasets_bp.route("/datasets",
          methods=["GET", "POST"])
def index():

    if request.method == "POST":

        try:

            create_dataset(
                request.form["name"],
                request.form["description"]
            )

            flash(
                "Dataset created.",
                "success"
            )

        except ValueError as e:

            flash(
                str(e),
                "error"
            )

        return redirect(
            url_for("datasets.index")
        )

    datasets = get_all_datasets()

    query = request.args.get("q")

    if query:
        datasets = [
            d for d in datasets
            if query.lower() in str(d.name).lower()
        ]

    return render_template(
        "datasets.html",
        datasets=datasets
    )

@datasets_bp.route(
    "/datasets/<uuid:dataset_id>/delete",
    methods=["POST"]
)
def delete(dataset_id):

    success = delete_dataset(dataset_id)
    
    if not success:

        flash(
            "Dataset not found.",
            "error"
        )

        return redirect(
            url_for("datasets.index")
        )

    flash(
        "Dataset deleted.",
        "success"
    )

    return redirect(
        url_for("datasets.index")
        )


