# pylint: disable=no-value-for-parameter
from pathlib import Path
import toml
import click

from escher import Builder


@click.group()
def cli():
    pass

@cli.command(name="list")
@click.argument("model_path", type=click.Path(exists=True))
def list_experiments(model_path):
    """
   Lists all the possible experiment IDs that are in the specified <model>.maud.toml file
    """
    model = toml.load(model_path)
    for e in model["experiments"]:
        print(e["id"])

@cli.command()
@click.argument("model_path", type=click.Path(exists=True))
@click.argument("experiment")
@click.option("--map-path",'map_json_path', help="For custom maps specify path to json file")
@click.option("--map-name", "map_name", help="For built-in maps specify name")
@click.option("--big-fonts", is_flag=True, help="Tweak css for escher to make fonts bigger")
@click.option("-O","--output", help="Name of html file to store output. Default is <model-name>_<experiment_id>")
def viz(model_path, experiment, map_json_path, map_name, big_fonts, output):
    """
    Take model.maud.toml file, one experiment and overlay it onto specified escher-compatible json map.
    Currently it's only one experiment per time. 
    """
    model = toml.load(model_path)
    exp_data = next((x for x in model["experiments"] if x["id"]==experiment), None)

    if experiment is not None:
        # Reaction data is in experimental measurements
        reaction_data = {x["target_id"]:x["value"] for x in exp_data["reaction_measurements"]}
        exp_metabolite_data = {x["target_id"]:x["value"] for x in exp_data["metabolite_measurements"]}

        # Enzyme abundance and unbalanced metabolites are now also specified in experimental measurements
        metabolite_data = exp_metabolite_data    

        # Enzymes are trickier
        # Priors are stored right now at reaction level, so we need to recover
        # Besides that Escher currently doesn't allow to disentangle proteomics and fluxomics

        if output is None:
            output = f"{Path(model_path).stem}_{Path(map_json_path).stem}_{experiment}.html"
        

        # potentially change css for svg.escher-svg .node-label and .reaction-label to make biger font-size
        # it's super buggy right now

        big_fonts_css = """
        svg.escher-svg .reaction-label {
            font-size: 30px;
            fill: rgb(32, 32, 120);
            text-rendering: optimizelegibility;
        }
        svg.escher-svg .node-label {
            font-size: 20px;
        }
        """
        if big_fonts:
            css = big_fonts_css
        else:
            css = None

        if map_name:
            Builder(map_name=map_name, reaction_data=reaction_data, metabolite_data=metabolite_data, embedded_css=css).save_html(output)
        if map_json_path:            
            Builder(map_json=map_json_path, reaction_data=reaction_data, metabolite_data=metabolite_data, embedded_css=css).save_html(output)
    else:
        print("Unable to find experiment")



if __name__ == "__main__":
    cli()
