import os

import pandas as pd
import RowDataProcessor as em


"""
Simple data-class, which just encapsulate the entry data. 
"""
class Entry_data:
    __sequence = ''
    __abl = ''
    __transcrypt = ''
    __drug = ''

    @property
    def sequence_path(self):
        return self.__sequence

    @sequence_path.setter
    def sequence_path(self, path):
        if os.path.exists(path):
            self.__sequence = path
        else:
            raise OSError("No such file with path: " + path)

    @property
    def abl_path(self):
        return self.__abl

    @abl_path.setter
    def abl_path(self, path):
        if os.path.exists(path):
            self.__abl = path
        else:
            raise OSError("No such file with path: " + path)

    @property
    def transcrypt_name(self):
        return self.__transcrypt

    @transcrypt_name.setter
    def transcrypt_name(self, name):
        self.__transcrypt = name

    @property
    def explorable_drug_name(self):
        return self.__drug

    @explorable_drug_name.setter
    def explorable_drug_name(self, name):
        self.__drug = name


class Peptides_creator:

    def __init__(self):
        self._pept_tables = list()
        self._flag = False
        for i in range(15):
            self._pept_tables.append(list())

    def convert_to_peptide_tables(self, sequence_source: em.Sequence_entity, drug_and_mutation_dict: dict[str, set]):
        for drug_and_mutate in drug_and_mutation_dict.items():
            drug = drug_and_mutate[0]
            transcrypt = sequence_source.transcrypt_name
            mutations = drug_and_mutate[1]
            sequence_source.add_mutations_to_sequence_mapping(mutations)

            for mutation in mutations:
                left_part = ""
                right_part = ""
                pos = em.Row_data_utils.get_position_by_mutation(mutation)
                sequence_one_str = sequence_source.get_sequence_by_mutation(mutation)
                peptide = sequence_one_str[pos - 1]
                seq_len = len(sequence_one_str)

                if pos <= 15 or (pos + 15) > seq_len:
                    continue

                # counter += 1

                for i in range(1, 16):
                    left_part = sequence_one_str[pos - 1 - i]
                    right_part = sequence_one_str[pos - 1 + i]
                    peptide = left_part + peptide + right_part
                    self._pept_tables[i - 1].append(list([peptide, drug, mutation, transcrypt]))

        self._flag = True

    def peptides_to_csv(self, files_name_prefix, out_dir="out"):

        if self._flag:
            df_list = list()

            for table in self._pept_tables:
                df = pd.DataFrame(table, columns=['PEPTIDE', 'DRUG_NAME', 'AA_MUTATION', 'TRANSCRYPT'])
                df.index.name = 'NUMBERS'
                df_list.append(df)

            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

            for i in range(len(df_list)):
                finish_dir = out_dir + "/" + files_name_prefix + str(i + 1)
                if not os.path.exists(finish_dir):
                    os.mkdir(finish_dir)
                name = out_dir + "/" + files_name_prefix + str(i + 1) + "/" + files_name_prefix + str(i + 1) + ".csv"
                df_list[i].to_csv(name)
                print("---> Create new file: " + str(name))

        else:
            print("Empty peptide tables, please create peptide tables and try again!")

    @property
    def peptide_tables(self):
        return self._pept_tables
