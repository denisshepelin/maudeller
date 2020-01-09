# Usage compose_model.sh <model>.toml <path to prior data without trailing slash>
# Uses model name which is filename but without .toml as argument 
echo "Composing the model $1 ..."

model_file=$1
model_name=${model_file/".toml"/}
maud_file="$model_name.maud.toml"
additional_data_path=$2


echo "Will use $model_file"
echo "Creating $maud_file.."

cp $model_file $maud_file

echo "Adding experimental data (concentrations and fluxes) ..."

cat ${model_name}_experiments.toml >> $maud_file

echo "Adding priors ..."


echo "Adding kinetic parameters ..."
cat ${additional_data_path}/priors_kinetic_parameters.toml >> $maud_file


echo "Adding marginal dGs ..."
cat ${additional_data_path}/thermodynamic_priors.toml >> $maud_file


# echo "Adding unbalanced metabolites ..."
# cat ${model_name}_unbalanced_metabolites.toml >> $maud_file

# echo "Adding enzyme priors based on proteomics"
# cat ${model_name}_enzyme_priors.toml >> $maud_file