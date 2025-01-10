from flask import Flask, Response, render_template, request
from chemistry import r_enumeration_auto
from chemistry2 import smiles_to_svg

app = Flask(__name__)


@app.route("/", methods=["GET"])
def root():
    all_smiles = []
    with open("smiles.txt", "r") as file:
        content = file.read()
        all_smiles = content.split("\n")
    return render_template("index.html", all_smiles=all_smiles)


@app.route(
    "/r_enumeration_auto_caller", methods=["GET"], endpoint="r_enumeration_auto_caller"
)
def r_enumeration_auto_caller():
    smiles = request.args.get("smiles", "")
    results = r_enumeration_auto(smiles)
    with open("smiles.txt", "r") as file:
        content = file.read()
        all_smiles = content.split("\n")
    return render_template("resultpage.html", results=results)


@app.route("/compound-image", methods=["GET"], endpoint="get_compound_image")
def get_compound_image():
    smiles = request.args.get("smiles", "")
    return Response(smiles_to_svg(smiles), mimetype="image/svg+xml")


@app.route("/scaffoldpage", methods=["GET"])
def x():
    all_smiles = []
    with open("smiles.txt", "r") as file:
        content = file.read()
        all_smiles = content.split("\n")
    return render_template("scaffoldpage.html", all_smiles=all_smiles)


@app.route("/proteinpage", methods=["GET"])
def p():
    all_smiles = []
    with open("smiles.txt", "r") as file:
        content = file.read()
        all_smiles = content.split("\n")
    return render_template("proteinpage.html", all_smiles=all_smiles)


@app.route("/resultpage", methods=["GET"])
def r():
    all_smiles = []
    with open("smiles.txt", "r") as file:
        content = file.read()
        all_smiles = content.split("\n")
    return render_template("resultpage.html", all_smiles=all_smiles)


@app.route("/renumpage", methods=["GET"])
def ren():
    all_smiles = []
    with open("smiles.txt", "r") as file:
        content = file.read()
        all_smiles = content.split("\n")
    return render_template("renumpage.html", all_smiles=all_smiles)
