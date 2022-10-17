import pandas as pd

class Sequence_tools:

    @staticmethod
    def get_replased_table(path_to_csv, drug):
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
    def create_mutation_dict(positions, old_letters, new_letters):
        mutation_dict = dict()
        for i in range(len(positions)):
            mutation = "p." + old_letters[i] + positions[i] + new_letters[i]
            mutation_dict[positions[i]] = mutation
        return mutation_dict


class Sequence_entity:
    __replaced_dict = dict()

    def __init__(self, path_to_sequence_file, transcrypt_name):
        self.__original_sequence = Sequence_tools.get_sequence_str(path_to_sequence_file)
        self.__transcrypt_name = transcrypt_name

    def get_sequence_by_position(self, position):
        return self.__replaced_dict[position]

    @property
    def original_sequence(self):
        return self.__original_sequence

    @property
    def transcrypt_name(self):
        return self.__transcrypt_name

    def create_replaced_dict(self, replaced_positions, new_letters):
        for i in range(len(replaced_positions)):
            new_seq = Sequence_tools.replace_letter(replaced_positions[i], new_letters[i], self.__original_sequence)
            self.__replaced_dict[replaced_positions[i]] = new_seq


class ABL_entity:
    __default_drug_name = "Imatinib"

    def __init__(self, path_to_ABL_csv, drug):
        if drug:
            self.__drug_name = drug
        else:
            self.__drug_name = ABL_entity.__default_drug_name

        replaced_table = Sequence_tools.get_replased_table(path_to_ABL_csv, self.__drug_name)
        self.__replaced_positions = replaced_table['position'].tolist()
        self.__replaced_letters = replaced_table['replaced_letter'].tolist()
        self.__new_letters = replaced_table['new_letter'].tolist()
        self.__mutation_dict = Sequence_tools.create_mutation_dict(self.__replaced_positions, self.__replaced_letters, self.__new_letters)

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

    def get_mutation_by_position(self, position):
        return self.__mutation_dict[position]