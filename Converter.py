import SeqToSDF
import json


path_to_config = "config\\1000.json"
path_to_peptides = "C:\\Users\\Zh_Nikolay\\Desktop\\peptides\\"
file_name_prefix = "pept_"
out_dir_name = "out"

with open(path_to_config, 'r') as config:
    properties = config.read()

prop_dict = json.loads(properties)
prop_dict["output"] = out_dir_name


def config_redactor(number):
    peptide_file_name = "/" + file_name_prefix + str(number) + ".csv"
    print("\nCurrent file: " + peptide_file_name)
    input_property = path_to_peptides + peptide_file_name
    prop_dict["input"] = input_property
    prop_dict["output"] = out_dir_name
    new_json = json.dumps(prop_dict)

    with open(path_to_config, 'w') as config:
        config.write(new_json)


def start():
    for i in range(15):
        config_redactor(i+1)
        SeqToSDF.start_fun(path_to_config)
