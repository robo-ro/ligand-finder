from flask import Flask, Response, render_template, request
from chemistry2 import smiles_to_svg

app = Flask(__name__)


@app.route("/", methods=["GET"])
def root():
    all_smiles = []
    with open("smiles.txt", "r") as file:
        content = file.read()
        all_smiles = content.split("\n")
    return render_template("index.html", all_smiles=all_smiles)


@app.route("/r_enumeration_auto", methods=["GET"], endpoint="r_enumeration_auto")
def r_enumeration_auto():
    smiles = request.args.get("smiles", "")
    return 'hello'

@app.route("/compound-image", methods=["GET"], endpoint="get_compound_image")
def get_compound_image():
    smiles = request.args.get("smiles", "")
    return Response(smiles_to_svg(smiles), mimetype="image/svg+xml")