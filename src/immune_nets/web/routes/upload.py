from flask import Blueprint
from flask import render_template
from flask import request
from flask import redirect
from flask import url_for

upload_bp = Blueprint("upload", __name__)

@upload_bp.route("/upload", methods=["GET", "POST"])
def index():

    if request.method == "POST":

        uploaded_file = request.files["dataset"]

        # call your import service
        # dataset_service.import_file(uploaded_file)

        return redirect(
            url_for("datasets.index")
        )

    return render_template(
        "upload.html"
    )
