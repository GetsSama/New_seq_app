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


class Peptides:

    def __init__(self, entry_data=None):
        self.__pept_tables = list()
        self.__is_tables_created = False

        for i in range(15):
            self.__pept_tables.append(list())

        if entry_data is not None:
            self.__this_entry_data = entry_data

    def __parse_entry_data(self):
        entry_data = self.__this_entry_data
        self.__sequence = em.Sequence_entity(entry_data.sequence_path, entry_data.transcrypt_name)
        self.__abl_data = em._Mutations_ABL(entry_data.abl_path, entry_data.explorable_drug_name)
        self.__sequence.create_replaced_dict(self.__abl_data.get_drug_and_mutations_dict)

    def create_peptide_tables_by_entry_data(self):
        self.__parse_entry_data()
        if not self.__is_tables_created:
            drug = self.__abl_data.drug_name
            transcrypt = self.__sequence.transcrypt_name

            for mutation in self.__abl_data.get_drug_and_mutations_dict:
                left_part = ""
                right_part = ""
                pos = int(self.__abl_data.get_pos_by_mutation(mutation))
                sequence_one_str = self.__sequence.get_sequence_by_mutation(mutation)
                peptide = sequence_one_str[pos - 1]
                seq_len = len(sequence_one_str)

                if pos <= 15 or (pos + 15) > seq_len:
                    continue

                # counter += 1

                for i in range(1, 16):
                    left_part = sequence_one_str[pos - 1 - i]
                    right_part = sequence_one_str[pos - 1 + i]
                    peptide = left_part + peptide + right_part
                    self.__pept_tables[i - 1].append(list([peptide, drug, mutation, transcrypt]))

                self.__is_tables_created = True

    def create_peptide_tables(self, abl_inst, sequence_inst):
        drug = abl_inst.drug_name
        transcrypt = sequence_inst.transcrypt_name

        for mutation in abl_inst.get_drug_and_mutations_dict:
            left_part = ""
            right_part = ""
            pos = int(abl_inst.get_pos_by_mutation(mutation))
            sequence_one_str = sequence_inst.get_sequence_by_mutation(mutation)
            peptide = sequence_one_str[pos - 1]
            seq_len = len(sequence_one_str)

            if pos <= 15 or (pos + 15) > seq_len:
                continue

            # counter += 1

            for i in range(1, 16):
                left_part = sequence_one_str[pos - 1 - i]
                right_part = sequence_one_str[pos - 1 + i]
                peptide = left_part + peptide + right_part
                self.__pept_tables[i - 1].append(list([peptide, drug, mutation, transcrypt]))

            self.__is_tables_created = True

    def peptides_to_csv(self, files_name, out_dir="out"):

        if self.__is_tables_created:
            df_list = list()

            for table in self.__pept_tables:
                df = pd.DataFrame(table, columns=['PEPTIDE', 'DRUG_NAME', 'AA_MUTATION', 'TRANSCRYPT'])
                df.index.name = 'NUMBERS'
                df_list.append(df)

            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

            for i in range(len(df_list)):
                name = out_dir + "/" + files_name + str(i + 1) + ".csv"
                df_list[i].to_csv(name)
                print("---> Create new file: " + str(name))

        else:
            print("Empty peptide tables, please create peptide tables and try again!")

    @property
    def peptide_tables(self):
        return self.__pept_tables
