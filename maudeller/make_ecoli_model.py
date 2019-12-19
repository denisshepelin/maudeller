# pylint: disable=no-value-for-parameter
import cobra
import toml
import click
import re

full_reactions = [
    "PGI",
    "FBP",
    "FBA",
    "TPI",
    "GAPD",
    "PGK",
    "PGM",
    "ENO",
    "G6PDH2r",
    "PGL",
    "GND",
    "RPI",
    "RPE",
    "TKT1",
    "TKT2",
    "TALA",
    "EDD",
    "EDA",
]

noloop_reactions = ["PGI", "G6PDH2r", "PGL", "GND", "RPI", "RPE", "TKT1"]

loop_reactions = ["PGI", "G6PDH2r", "PGL", "GND", "RPI", "RPE", "TKT1", "TALA"]

mechanism_map = {
    "PGI": "uniuni",
    "FBP": "uniuni",
    "FBA": "uniuni",
    "TPI": "ordered_unibi",
    "GAPD": "uniuni",
    "PGK": "ordered_bibi",
    "PGM": "uniuni",
    "ENO": "uniuni",
    "G6PDH2r": "ordered_bibi",
    "PGL": "uniuni",
    "GND": "ordered_bibi",  # CO2 is ignored
    "RPI": "uniuni",
    "RPE": "uniuni",
    "TKT1": "ordered_bibi",
    "TKT2": "ordered_bibi",
    "TALA": "ordered_bibi",
    "EDD": "uniuni",
    "EDA": "ordered_unibi",
}


def make_backbone(model, reactions, unwanted_metabolites, unbalanced_metabolites):
    """
    model - cobra Model object

    reactions - list of reaction IDs that are in this model

    unwanted_metabolites - list of metabolites that are going to be removed from reactions
    """
    print("Loading model...")
    print(model)

    metabolites = set()

    # We want to make a dict for toml
    toml_dict = {}

    toml_dict["compartments"] = [{"id": "c", "name": "cytosol", "volume": 1}]
    
    toml_dict["reactions"] = []
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

        toml_dict["reactions"].append(rdict)

    toml_dict["metabolites"] = []
    for mid in metabolites:
        if mid not in unwanted_metabolites:
            m = model.metabolites.get_by_id(mid)
            mdict = {}
            match = re.match(r"(?P<met_id>\S+)_(?P<compartment>\w)$",m.id).groupdict()
            mdict["id"] = match["met_id"]
            mdict["name"] = m.name
            mdict["compartment"] = match["compartment"]
            if m.id in unbalanced_metabolites:
                mdict["balanced"] = False
            else:
                mdict["balanced"] = True            
            toml_dict["metabolites"].append(mdict)

    return toml_dict


@click.command()
@click.option(
    "-P",
    "--preset",
    type=click.Choice(["noloop", "ecoli", "loop"]),
    help="Allows to choose from some prespecified set of reaction lists",
)
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
    help="Set of metabolite IDs to be remove, by default it's h2o, H+, Pi. Use None as option to bypass all removals",
    default="h2o_c,h_c,pi_c",
)
@click.option(
    "--cobra-model-path",
    help="Path to the COBRA model",
    default="/Users/denshe/Work/DataAnalysis/DataIntegrationProject/models/iML1515.json",
    type=click.Path(exists=True),
)
@click.option(
    "-o",
    "--output",
    help="Filename where model will be writen in toml format. Default is output.toml",
    default="output.toml",
)
def make_model(
    preset=None,
    reactions=None,
    unbalanced_metabolites=None,
    remove_metabolites=["h2o_c", "h_c", "pi_c"],
    cobra_model_path="/Users/denshe/Work/DataAnalysis/DataIntegrationProject/models/iML1515.json",
    output="output.toml",
):
    # Check types of inputs and convert them to lists!

    cobra_model = cobra.io.load_json_model(cobra_model_path)

    # Decide what reactions are we going to use
    if reactions is not None:
        if type(reactions) is str:
            reaction_list = reactions.split(",")
        print(f"Going to work on this {len(reaction_list)} reactions:")
        print(",".join(reaction_list))
    else:
        if preset is not None:
            print(f"Using preset {preset}")
            if preset == "noloop":
                reaction_list = noloop_reactions
            elif preset == "ecoli":
                reaction_list = full_reactions
            elif preset == "loop":
                reaction_list = loop_reactions

    # Are there any metabolites that should be filtered out?
    if type(remove_metabolites) is str:
        if remove_metabolites == "None":
            remove_metabolites = []
        else:
            remove_metabolites = remove_metabolites.split(",")

    # In theory there should be at least some unbalanced metabolites
    if type(unbalanced_metabolites) is str:
        unbalanced_metabolites = unbalanced_metabolites.split(",")
    if len(unbalanced_metabolites) == 0:
        print("There are no unbalanced metabolites, that is strange!")

    model_dict = make_backbone(
        cobra_model, reaction_list, remove_metabolites, unbalanced_metabolites
    )

    with open(output, "w") as out:
        toml.dump(model_dict, out)


if __name__ == "__main__":
    make_model()

