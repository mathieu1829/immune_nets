from flask import (
    Blueprint,
    render_template,
    request,
    redirect,
    url_for,
    flash
)

from pathlib import Path

from immune_nets.services.dataset_service import (
    get_all_datasets,
)

from immune_nets.services.repertoire_service import (
    add_repertoire_to_dataset,
)

import pandas as pd
import io

upload_bp = Blueprint("upload", __name__)

@upload_bp.route("/upload", methods=["GET", "POST"])
def index():

    datasets = get_all_datasets()

    if request.method == "POST":

        dataset_id = request.form.get(
            "dataset_id"
        )

        form_repertoire_name = request.form.get(
            "repertoire_name"
        )

        repertoire_name = form_repertoire_name if form_repertoire_name is not None else ""

        form_repertoire_desc = request.form.get(
            "repertoire_desc"
        )


        repertoire_desc = form_repertoire_desc if form_repertoire_desc is not None else ""


        use_filename = (
            request.form.get(
                "use_filename"
            )
            is not None
        )

        uploaded_file = request.files.get(
            "repertoire_file"
        )

        if use_filename and uploaded_file is not None:
            repertoire_name = Path(uploaded_file.filename).stem

        if not dataset_id:

            flash(
                "Please select a dataset.",
                "error"
            )

            return render_template(
                "upload.html",
                datasets=datasets
            )

        if not uploaded_file:

            flash(
                "No file selected.",
                "error"
            )

            return render_template(
                "upload.html",
                datasets=datasets
            )

        if not uploaded_file.filename.endswith(
            ".csv"
        ):

            flash(
                "Only CSV files are allowed.",
                "error"
            )

            return render_template(
                "upload.html",
                datasets=datasets
            )

        content = uploaded_file.read()
        df = pd.read_csv(io.BytesIO(content))
        if not set(["proportion", "cdr3s_aa"]).issubset(set(df.columns)):
            flash(
                "Wrong csv structure",
                "error"
            )

            return render_template(
                "upload.html",
                datasets=datasets
            )
        uploaded_file.seek(0)

        add_repertoire_to_dataset(
            dataset_id,
            repertoire_name,
            repertoire_desc,
            uploaded_file
        )

        flash(
            "Repertoire uploaded.",
            "success"
        )

        return redirect(
            url_for("repertoires.dataset_view",
                    dataset_id=dataset_id)
        )

    return render_template(
        "upload.html",
        datasets=datasets
    )
