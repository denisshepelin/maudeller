# pylint: disable=no-value-for-parameter
import cobra
import toml
import click
import re
from pathlib import Path


def get_reactions_dict(model, reactions, unwanted_metabolites, gpr):
    """
    Given COBRA model, reactions and metabolites that should not be included in the Maud model
    prepare python dict containing information about reactions and set of metabolites that are in those reactions
    """
    result = []
    metabolites = set()
    for rid in reactions:
        r = model.reactions.get_by_id(rid)
        metabolites.update([m.id for m in r.metabolites])
        rdict = {}
        rdict["id"] = rid
        rdict["name"] = r.name
        # add stoichiometry
        rdict["stoichiometry"] = {}
        for m, sign in r.metabolites.items():
            if m.id not in unwanted_metabolites:
                rdict["stoichiometry"][m.id] = int(sign)
        # add enzymes information
        if len(r.genes) == 1:
            gene_name = r.gene_name_reaction_rule
        else:
            if rid in gpr:
                gene_name = gpr[rid]
            else:
                gene_name = "FIX"
                print(
                    f"Reaction with ID {rid} has no gene annotation potentially because there are too many genes to set it automatically. Candidates are: {r.gene_name_reaction_rule}"
                )
        rdict["enzymes"] = [
            {
                "id": rid,
                "name": r.name,
                "mechanism": "modular_rate_law",
                "gene": gene_name,
            }
        ]
        # It's possible to specify mechanism for each reaction
        # rdict["enzymes"] = [{"id" : rid, "name": r.name, "mechanism": mechanism_map[rid]}]
        result.append(rdict)

    return (result, metabolites)


def get_metabolites_dict(
    model, metabolites, unwanted_metabolites, unbalanced_metabolites
):
    """
    Given model, metabolite IDs and list of metabolites to be excluded and unbalanced metabolites
    Prepare python dicts representing metabolites information
    """
    result = []
    for mid in metabolites:
        if mid not in unwanted_metabolites:
            m = model.metabolites.get_by_id(mid)
            mdict = {}
            match = re.match(r"(?P<met_id>\S+)_(?P<compartment>\w)$", m.id).groupdict()
            mdict["id"] = match["met_id"]
            mdict["name"] = m.name
            mdict["compartment"] = match["compartment"]
            if m.id in unbalanced_metabolites:
                mdict["balanced"] = False
            else:
                mdict["balanced"] = True
            result.append(mdict)
    return result


def get_priors(
    model_dict,
    data_path,
    kinetic_file="priors_kinetic_parameters.toml",
    thermodynamics_file="thermodynamic_priors.toml",
):
    """
    Given model reactions and set of metabolites 
    select from database relevant parts of kinetic laws and formation energies
    """

    kinetic_toml = toml.load(Path(data_path) / kinetic_file)
    thermodynamic_toml = toml.load(Path(data_path) / thermodynamics_file)

    # select kinetic parameters for reactions. They are stored in
    # ["priors"]["kinetic_parameters"][REACTION_ID]
    rxn_ids = [r["id"] for r in model_dict["reactions"]]
    kinetic_parameters = {
        rid: kinetic_toml["priors"]["kinetic_parameters"][rid] for rid in rxn_ids
    }

    # select formation energies. They are stored in
    # ["priors"]["thermodynamic_parameters"]["formation_energies"] list
    met_ids = [m["id"] for m in model_dict["metabolites"]]
    thermodynamic_parameters = [
        fe
        for fe in thermodynamic_toml["priors"]["thermodynamic_parameters"][
            "formation_energies"
        ]
        if fe["target_id"] in met_ids
    ]
    return (kinetic_parameters, thermodynamic_parameters)


def make_backbone(
    model, reactions, unwanted_metabolites, unbalanced_metabolites, priors_data_path, gpr
):
    """
    model - cobra Model object

    reactions - list of reaction IDs that are in this model

    unwanted_metabolites - list of metabolites that are going to be removed from reactions
    """
    print("Loading model...")
    print(model)

    # We want to make a dict for toml
    toml_dict = {}

    toml_dict["compartments"] = [{"id": "c", "name": "cytosol", "volume": 1}]

    toml_dict["reactions"], metabolites = get_reactions_dict(
        model, reactions, unwanted_metabolites, gpr
    )

    toml_dict["metabolites"] = get_metabolites_dict(
        model, metabolites, unwanted_metabolites, unbalanced_metabolites
    )

    toml_dict["priors"] = {}
    toml_dict["priors"]["thermodynamic_parameters"] = {}
    toml_dict["priors"]["kinetic_parameters"], toml_dict["priors"][
        "thermodynamic_parameters"
    ]["formation_energies"] = get_priors(toml_dict, priors_data_path)

    return toml_dict


@click.command()
@click.option(
    "-R",
    "--reactions",
    nargs=1,
    help="Set of reaction IDs separated by comma looking like this: PGI,GND,TPI",
)  # sadly CLI passes string instead of list because Click doesn't directly support list of options by design
@click.option(
    "-U",
    "--unbalanced-metabolites",
    nargs=1,
    help="Set of metabolite IDs which will be marked as unbalanced ",
)
@click.option(
    "--remove-metabolites",
    nargs=1,
    help="Set of metabolite IDs to be remove, by default it's h2o, H+. Use None as option to bypass all removals",
    show_default=True,
    default="h2o_c,h_c",
)
@click.option(
    "--cobra-model-path",
    help="Path to the COBRA model",
    default="/Users/denshe/Work/DataAnalysis/DataIntegrationProject/models/iML1515.json",
    type=click.Path(exists=True),
)
@click.option(
    "--priors-data-path",
    help="Path to the folder containing prior information",
    default="/Users/denshe/Work/KineticalModeling/maudeller/data/common/",
    type=click.Path(exists=True),
)
@click.option(
    "-C",
    "--toml-config",
    "toml_config",
    help="Load configurtation from toml file",
    type=click.Path(exists=True)
)
@click.option(
    "-o",
    "--output",
    help="Filename where model will be writen in toml format. Default is output.toml",
    default="output.toml",
)
def make_model(
    reactions=None,
    unbalanced_metabolites=None,
    remove_metabolites=["h2o_c", "h_c"],
    cobra_model_path="/Users/denshe/Work/DataAnalysis/DataIntegrationProject/models/iML1515.json",
    priors_data_path="/Users/denshe/Work/KineticalModeling/maudeller/data/common/",
    toml_config=None,
    output="output.toml",
):
    # Check if there is configuration file, if there is then load config from there
    if toml_config is not None:
        config = toml.load(Path(toml_config))
        if "reactions" in config:
            reactions = config["reactions"]
        if "unbalanced_metabolites" in config:
            unbalanced_metabolites = config["unbalanced_metabolites"]

    gpr = config.get("gpr",[])
    cobra_model = cobra.io.load_json_model(cobra_model_path)

    # Check types of inputs and convert them to lists!
    # Decide what reactions are we going to use
    if reactions is not None:
        if type(reactions) is str:
            reaction_list = reactions.split(",")
        if type(reactions) is list:
            reaction_list = reactions
        print(f"Going to work on this {len(reaction_list)} reactions:")
        print(",".join(reaction_list))

    # Are there any metabolites that should be filtered out?
    if type(remove_metabolites) is str:
        if remove_metabolites == "None":
            remove_metabolites = []
        else:
            remove_metabolites = remove_metabolites.split(",")

    # In theory there should be at least some unbalanced metabolites
    if type(unbalanced_metabolites) is str:
        unbalanced_metabolites = unbalanced_metabolites.split(",")
    if unbalanced_metabolites is None:
        print("There are no unbalanced metabolites, that is strange!")

    model_dict = make_backbone(
        cobra_model,
        reaction_list,
        remove_metabolites,
        unbalanced_metabolites,
        priors_data_path,
        gpr,
    )

    with open(output, "w") as out:
        toml.dump(model_dict, out)


if __name__ == "__main__":
    make_model()

