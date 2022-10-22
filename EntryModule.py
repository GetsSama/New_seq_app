import pandas as pd

class Sequence_tools:

    @staticmethod
    def get_replaced_table(path_to_csv, drug):
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

        df_uni_commands = pd.DataFrame(unique_comm_table, columns=['position', 'replaced_letter', 'new_letter']) # Позиция, заменяемая буква, заменяющая буква
        return df_uni_commands

    @staticmethod
    def get_sequence_str(path_to_sequence):
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

        return  pos, old_let, new_let

class Sequence_entity:
    _replaced_dict = dict()

    def __init__(self, path_to_sequence_file, transcrypt_name):
        self.__original_sequence = Sequence_tools.get_sequence_str(path_to_sequence_file)
        self.__transcrypt_name = transcrypt_name

    def get_sequence_by_mutation(self, mutation):
        return self._replaced_dict[mutation]

    @property
    def original_sequence(self):
        return self.__original_sequence

    @property
    def transcrypt_name(self):
        return self.__transcrypt_name

    def create_replaced_dict(self, replaced_positions, new_letters, mutations):
        for i in range(len(replaced_positions)):
            new_seq = Sequence_tools.replace_letter(replaced_positions[i], new_letters[i], self.__original_sequence)
            self._replaced_dict[mutations[i]] = new_seq


class ABL_entity:
    __default_drug_name = "Imatinib"
    __mutation_pos_dict = None

    def __init__(self, path_to_ABL_csv, drug):
        if drug:
            self.__drug_name = drug
        else:
            self.__drug_name = ABL_entity.__default_drug_name

        replaced_table = Sequence_tools.get_replaced_table(path_to_ABL_csv, self.__drug_name)
        self.__replaced_positions = replaced_table['position'].tolist()
        self.__replaced_letters = replaced_table['replaced_letter'].tolist()
        self.__new_letters = replaced_table['new_letter'].tolist()
        self.__mutations = Sequence_tools.create_mutations(self.__replaced_positions, self.__replaced_letters, self.__new_letters)

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
    def mutations(self):
        return self.__mutations

    def get_pos_by_mutation(self, mutation):
        if self.__mutation_pos_dict is None:
            self.__mutation_pos_dict = Sequence_tools.create_pos_mutation(self.__replaced_positions, self.__mutations)
        return self.__mutation_pos_dict[mutation]

    def remove_mutations(self, removed_mutations):
        for mutation in removed_mutations:
            pos, old_let, new_let = Sequence_tools.mutations_parser(mutation)
            try:
                self.__mutations.remove(mutation)
                self.__replaced_positions.remove(pos)
                self.__replaced_letters.remove(old_let)
                self.__new_letters.remove(new_let)
            except ValueError:
                pass

        self.__mutation_pos_dict = Sequence_tools.create_pos_mutation(self.__replaced_positions, self.__mutations)


class ABL_table_analyzer:

    @staticmethod
    def get_unique_values_in_column(path_to_table, column_name):
        abl_table = pd.read_csv(path_to_table)
        drug_column = abl_table[column_name]
        unique_list = drug_column.unique().tolist()
        return unique_list