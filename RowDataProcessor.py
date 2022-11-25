from abc import abstractmethod, ABCMeta

import pandas as pd

"""
    Utility class for work module API
"""


class _Sequence_tools:

    @staticmethod
    def get_replaced_table_for_ABL(path_to_csv, drug):
        df = pd.read_csv(path_to_csv)
        # Выделяем только смену позиций для иматиниба
        right_data = pd.DataFrame(df[[' DRUG_NAME', ' AA_MUTATION']])
        filter_name = right_data[' DRUG_NAME'] == drug
        data_imat = right_data.loc[filter_name]
        replace_mass = data_imat[' AA_MUTATION']
        # Оставляем только сами АК и позиции
        command_mass = []

        for st in replace_mass:
            command_mass.append(st[2:])

        # print(command_mass)

        # Оставляем только замененные АК
        count = 0
        for st in command_mass:
            for let in st[1:]:
                if not let.isdigit():
                    count += 1
            if count > 1:
                command_mass.remove(st)
            count = 0

        for st in command_mass:
            if st == '?':
                command_mass.remove(st)

        # Создаем таблицу с уникальными значениями
        unique_mass = pd.DataFrame(command_mass)[0].unique()

        unique_comm_table = []
        for st in unique_mass:
            unique_comm_table.append(list([st[1:-1], st[0], st[-1]]))

        df_uni_commands = pd.DataFrame(unique_comm_table, columns=['position', 'replaced_letter',
                                                                   'new_letter'])  # Позиция, заменяемая буква, заменяющая буква
        return df_uni_commands

    @staticmethod
    def get_replaced_table_for_any(path_to_csv, drug_column, mutations_column, drug_name):
        df = pd.read_csv(path_to_csv)
        # Выделяем только смену позиций
        right_data = pd.DataFrame(df[[drug_column, mutations_column]])
        filter_name = right_data[drug_column] == drug_name
        data_imat = right_data.loc[filter_name]
        replace_mass = data_imat[mutations_column]
        # Оставляем только сами АК и позиции
        command_mass = []

        for st in replace_mass:
            command_mass.append(st[2:])

        # print(command_mass)

        # Оставляем только замененные АК
        count = 0
        for st in command_mass:
            for let in st[1:]:
                if not let.isdigit():
                    count += 1
            if count > 1:
                command_mass.remove(st)
            count = 0

        for st in command_mass:
            if st == '?':
                command_mass.remove(st)

        # Создаем таблицу с уникальными значениями
        unique_mass = pd.DataFrame(command_mass)[0].unique()

        unique_comm_table = []
        for st in unique_mass:
            unique_comm_table.append(list([st[1:-1], st[0], st[-1]]))

        df_uni_commands = pd.DataFrame(unique_comm_table, columns=['position', 'replaced_letter',
                                                                   'new_letter'])  # Позиция, заменяемая буква, заменяющая буква
        return df_uni_commands

    @staticmethod
    def get_sequence_one_str(path_to_sequence):
        with open(path_to_sequence, 'r') as source_data:
            lines = source_data.readlines()
        sequence_one_str = ""
        for i in range(1, len(lines)):
            sequence_one_str += lines[i][:-1]

        return sequence_one_str

    @staticmethod
    def replace_letter(position, new_letter, sequence):
        position = int(position)
        new_str = ""
        new_str += sequence[:(position - 1)]
        new_str += new_letter
        new_str += sequence[position:]
        return new_str

    @staticmethod
    def create_mutations(positions, old_letters, new_letters):
        mutations_list = list()
        for i in range(len(positions)):
            mutation = "p." + old_letters[i] + positions[i] + new_letters[i]
            mutations_list.append(mutation)
        return mutations_list

    @staticmethod
    def create_pos_mutation(positions, mutations):
        mutations_dict = dict()
        for i in range(len(positions)):
            mutation = mutations[i]
            mutations_dict[mutation] = positions[i]
        return mutations_dict

    @staticmethod
    def mutations_parser(mutation):
        st = mutation[2:]
        data = list([st[1:-1], st[0], st[-1]])
        pos = data[0]
        old_let = data[1]
        new_let = data[2]

        return pos, old_let, new_let


"""
    Class for work with sequence. 
"""


class Sequence_entity:
    _replaced_dict = dict()

    def __init__(self, path_to_sequence_file, transcrypt_name):
        self.__original_sequence = _Sequence_tools.get_sequence_one_str(path_to_sequence_file)
        self._transcrypt_name = transcrypt_name

    def get_sequence_by_mutation(self, mutation):
        return self._replaced_dict[mutation]

    @property
    def original_sequence(self):
        return self.__original_sequence

    @property
    def transcrypt_name(self):
        return self._transcrypt_name

    def get_mutation_to_sequence_mapping(self):
        return self._replaced_dict

    def add_mutations_to_sequence_mapping(self, mutations):
        for mutation in mutations:
            replaced_position, old_letter, new_letter = _Sequence_tools.mutations_parser(mutation)
            new_seq = _Sequence_tools.replace_letter(replaced_position, new_letter, self.__original_sequence)
            self._replaced_dict[mutation] = new_seq


class Mutations_data_source(metaclass=ABCMeta):
    @abstractmethod
    def get_drug_and_mutations_dict(self) -> dict[str, set]:
        pass


class Mutations_data_factory:
    @staticmethod
    def create_mutations_data_source(path_to_source_table, table_pattern="ABL", mutation_column="default",
                                     drug_column="default"):
        if table_pattern == "ABL":
            return _Mutations_ABL(path_to_source_table)
        elif table_pattern == "AnyTable":
            return _Mutations_any_table(path_to_source_table, mutation_column, drug_column)
        else:
            raise Exception("Unknown table pattern!")

class _Mutations_ABL_one_drug:
    _default_drug_name = "Imatinib"
    _mutation_pos_dict = None

    def __init__(self, path_to_ABL_csv, drug=_default_drug_name):
        if not drug == self._default_drug_name:
            self.__drug_name = drug
        else:
            self.__drug_name = self._default_drug_name

        replaced_table = _Sequence_tools.get_replaced_table_for_ABL(path_to_ABL_csv, self.__drug_name)
        self.__replaced_positions = replaced_table['position'].tolist()
        self.__replaced_letters = replaced_table['replaced_letter'].tolist()
        self.__new_letters = replaced_table['new_letter'].tolist()
        self.__mutations = _Sequence_tools.create_mutations(self.__replaced_positions, self.__replaced_letters,
                                                            self.__new_letters)

    @property
    def replaced_positions(self):
        return self.__replaced_positions

    @property
    def replaced_letters(self):
        return self.__replaced_letters

    @property
    def new_letters(self):
        return self.__new_letters

    @property
    def drug_name(self):
        return self.__drug_name

    @property
    def get_mutations(self):
        return self.__mutations

    def get_pos_by_mutation(self, mutation):
        if self.__mutation_pos_dict is None:
            self.__mutation_pos_dict = _Sequence_tools.create_pos_mutation(self.__replaced_positions, self.__mutations)
        return self.__mutation_pos_dict[mutation]

    def remove_mutations(self, removed_mutations):
        for mutation in removed_mutations:
            pos, old_let, new_let = _Sequence_tools.mutations_parser(mutation)
            try:
                self.__mutations.remove(mutation)
                self.__replaced_positions.remove(pos)
                self.__replaced_letters.remove(old_let)
                self.__new_letters.remove(new_let)
            except ValueError:
                pass
        self.__mutation_pos_dict = _Sequence_tools.create_pos_mutation(self.__replaced_positions, self.__mutations)

class _Mutations_any_table(Mutations_data_source):

    def __init__(self, path_to_any_table_csv, mutations_column, drug_column):
        self.__mutations_column = mutations_column
        self.__drug_column = drug_column
        self.__path_to_source = path_to_any_table_csv
        self._drug_and_mutations_dict = dict()

    def with_drug_parser(self):
        list_uni_drugs = Row_data_utils.get_unique_values_in_column(self.__path_to_source, self.__drug_column)
        for drug in list_uni_drugs:
            replaced_table = _Sequence_tools.get_replaced_table_for_any(self.__path_to_source, self.__drug_column, self.__mutations_column, drug)
            mutations = _Sequence_tools.create_mutations(replaced_table['positions'], replaced_table['replaced_letter'], replaced_table['new_letter'])
            self._drug_and_mutations_dict[drug] = mutations
        Row_data_utils.map_to_dict_of_sets(self._drug_and_mutations_dict)

    def no_drug_parser(self):
        list_mutations = pd.read_csv(self.__path_to_source, sep=";")
        list_mutations = list_mutations[self.__mutations_column]
        self._drug_and_mutations_dict["No drug"] = list_mutations.tolist()
        Row_data_utils.map_to_dict_of_sets(self._drug_and_mutations_dict)

    def get_drug_and_mutations_dict(self) -> dict[str, set]:
        if self.__drug_column == "default":
            self.no_drug_parser()
        else:
            self.with_drug_parser()
        return self._drug_and_mutations_dict

class _Mutations_ABL(Mutations_data_source):
    def __init__(self, path_to_ABL_table_csv):
        self.__path_to_source = path_to_ABL_table_csv
        self.__all_drugs = Row_data_utils.get_unique_values_in_column(self.__path_to_source, " DRUG_NAME")
        self.__res_drug_and_mutations_dict = dict()
        self.__fill_drug_and_mutations_dict()

    def get_drug_and_mutations_dict(self) -> dict[str, set]:
        return self.__res_drug_and_mutations_dict

    def __fill_drug_and_mutations_dict(self):
        for drug in self.__all_drugs:
            for_drug = _Mutations_ABL_one_drug(self.__path_to_source, drug)
            mutations_for_drug = for_drug.get_mutations
            self.__res_drug_and_mutations_dict[drug] = mutations_for_drug
        Row_data_utils.map_to_dict_of_sets(self.__res_drug_and_mutations_dict)

class Drug_oriented_mutations(Mutations_data_source):

    def __init__(self, mutations_data_source: Mutations_data_source):
        self.__data_source = mutations_data_source
        self.__row_drug_and_mutations_dict = self.__data_source.get_drug_and_mutations_dict()
        self.__decorated_dict = dict()
        self.__concrete_drug = None

    def set_drug(self, drug):
        self.__concrete_drug = drug

    def _most_freq_drug(self):
        max_count = 0
        most_freq = ""
        for pair in self.__row_drug_and_mutations_dict.items():
            if len(pair[1]) > max_count:
                max_count = len(pair[1])
                most_freq = pair[0]
        return most_freq

    def _orient_row_data_to_drug(self):
        if self.__concrete_drug is not None:
            concrete_drug_mutations = self.__row_drug_and_mutations_dict[self.__concrete_drug]
        else:
            self.__concrete_drug = self._most_freq_drug()
            concrete_drug_mutations = self.__row_drug_and_mutations_dict[self.__concrete_drug]

        for pair in self.__row_drug_and_mutations_dict.items():
            if pair[0] == self.__concrete_drug:
                self.__decorated_dict[self.__concrete_drug] = concrete_drug_mutations
            else:
                unique_mutations = set()
                for mutate in pair[1]:
                    if mutate not in concrete_drug_mutations:
                        unique_mutations.add(mutate)
                self.__decorated_dict[pair[0]] = unique_mutations
        Row_data_utils.map_to_dict_of_sets(self.__decorated_dict)

    def get_drug_and_mutations_dict(self) -> dict[str, set]:
        if len(self.__decorated_dict) == 0:
            self._orient_row_data_to_drug()
        return self.__decorated_dict


class Row_data_utils:

    @staticmethod
    def get_unique_values_in_column(path_to_table, column_name):
        any_table = pd.read_csv(path_to_table)
        any_column = any_table[column_name]
        unique_list = any_column.unique().tolist()
        return unique_list

    @staticmethod
    def get_position_by_mutation(mutation) -> int:
        pos, old, new = _Sequence_tools.mutations_parser(mutation)
        return int(pos)

    @staticmethod
    def concatenate_drug_and_mutations_dicts(dict1: dict, dict2: dict) -> dict[str, set]:
        Row_data_utils.map_to_dict_of_sets(dict1)
        Row_data_utils.map_to_dict_of_sets(dict2)

        res_dict = dict1.copy()
        for key in dict2.keys():
            if key in res_dict.keys():
                res_dict[key].union(dict2[key])
            else:
                res_dict[key] = dict2[key]

        return res_dict

    @staticmethod
    def map_to_dict_of_sets(_dict_: dict) -> dict[str, set]:
        for _key in _dict_.keys():
            values = _dict_[_key]
            val_set = set(values)
            _dict_[_key] = val_set
        return _dict_


