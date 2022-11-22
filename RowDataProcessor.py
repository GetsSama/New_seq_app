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
        self.__transcrypt_name = transcrypt_name

    def get_sequence_by_mutation(self, mutation):
        return self._replaced_dict[mutation]

    @property
    def original_sequence(self):
        return self.__original_sequence

    @property
    def transcrypt_name(self):
        return self.__transcrypt_name

    @property
    def get_mutation_and_sequence_mapping(self):
        return self._replaced_dict

    def create_replaced_dict(self, mutations):
        for mutation in mutations:
            replaced_position, old_letter, new_letter = _Sequence_tools.mutations_parser(mutation)
            new_seq = _Sequence_tools.replace_letter(replaced_position, new_letter, self.__original_sequence)
            self._replaced_dict[mutation] = new_seq


class Mutations_data_source(metaclass=ABCMeta):
    @abstractmethod
    def get_mutations(self):
        pass

class Mutations_data_factory:
    @staticmethod
    def create_mutations_data_source(table_pattern):


class Mutations_ABL(Mutations_data_source):
    __default_drug_name = "Imatinib"
    __mutation_pos_dict = None

    def __init__(self, path_to_ABL_csv, drug=__default_drug_name):
        if not drug == Mutations_ABL.__default_drug_name:
            self.__drug_name = drug
        else:
            self.__drug_name = Mutations_ABL.__default_drug_name

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


class Table_analyzer:

    @staticmethod
    def get_unique_values_in_column(path_to_table, column_name):
        any_table = pd.read_csv(path_to_table)
        any_column = any_table[column_name]
        unique_list = any_column.unique().tolist()
        return unique_list
