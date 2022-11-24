import os.path

import SeqToSDF
import json

class SeqToSDF_manager:
    _path_to_config = "config\\1000.json"
    _path_to_peptides = "C:\\Desktop\\peptides\\"
    _file_name_prefix = "pept_"
    _out_dir_name = "out"

    def __init__(self):
        self._prop_dict = None

    def _load(self):
        with open(self._path_to_config, 'r') as config:
            properties = config.read()

        self._prop_dict = json.loads(properties)
        self._prop_dict["output"] = self._out_dir_name

    def set_manager_properties(self, path_to_config, path_to_peptides, peptides_name_prefix, out_dir_name="out"):
        if os.path.exists(path_to_config):
            self._path_to_config = path_to_config
        else:
            raise OSError("Do not exist file with name: " + path_to_config)
        if os.path.exists(path_to_peptides):
            self._path_to_peptides = path_to_peptides
        else:
            raise OSError("Do not exist file with name: " + path_to_peptides)
        if not os.path.exists(out_dir_name) and not out_dir_name == "out":
            os.mkdir(out_dir_name)
        self._out_dir_name = out_dir_name
        self._file_name_prefix = peptides_name_prefix
        self._load()

    def _config_redactor(self, number):
        peptide_file_name = "/" + self._file_name_prefix + str(number) + "/" + self._file_name_prefix + str(number) + ".csv"
        print("\nCurrent file: " + peptide_file_name)
        input_property = self._path_to_peptides + peptide_file_name
        self._prop_dict["input"] = input_property
        self._prop_dict["output"] = self._path_to_peptides + "/" + self._file_name_prefix + str(number)
        self._prop_dict["filename"] = self._file_name_prefix + str(number) + "_sdf"
        new_json = json.dumps(self._prop_dict)

        with open(self._path_to_config, 'w') as config:
            config.write(new_json)

    def start(self):
        for i in range(15):
            self._config_redactor(i+1)
            SeqToSDF.start_fun(self._path_to_config)
