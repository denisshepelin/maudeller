# maudeller
Tooling to support Maud kinetic modeling package

## make_backbone
This tool helps to generate maud-compatible backbones from COBRA model. This tool expects user to have a COBRA model in `json` format.

Run `python make_backbone --help"` to see list of potential options and arguments.

It's recommended to use .toml configs (--toml-config" option) that hold information about what backbone you want to create

## visualize_experiment
This tool takes full Maud model that contains experimental values and overlays experimental values on the Escher map.

Run `python visualize_experiment.py --help` to see list of options and arguments
