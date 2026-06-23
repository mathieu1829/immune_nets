from flask import Flask
from .routes import main_bp, upload_bp, datasets_bp, networks_bp, comparison_bp, repertoires_bp, clonotypes_bp

import immune_nets.db

def create_app():
    app = Flask(__name__)

    app.secret_key = "development-key"

    app.register_blueprint(main_bp)
    app.register_blueprint(upload_bp)
    app.register_blueprint(clonotypes_bp)
    app.register_blueprint(repertoires_bp)
    app.register_blueprint(datasets_bp)
    app.register_blueprint(networks_bp)
    app.register_blueprint(comparison_bp)


    return app

if __name__ == "__main__":
    app = create_app()
    app.run(debug=True)
