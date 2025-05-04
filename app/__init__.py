from flask import Flask
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

from app import routes  
# Register the Blueprint
from app.routes import bp
app.register_blueprint(bp)