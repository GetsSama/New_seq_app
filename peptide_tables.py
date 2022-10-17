import pandas as pd
from numpy import double


def replased_fun(path_to_csv, path_to_consiq):
    df_uni_commands = get_replased_table(path_to_csv)

    # Имя новой последовательности
    path_to_result = path_to_consiq[:-4] + '_new' + path_to_consiq[-4:]

    # Чтение исходного файла последовательности
    # -----------------------------------------------
    lines = list()
    with open(path_to_consiq, 'r') as source_data:
        lines = source_data.readlines()
    # -----------------------------------------------

    simbols_in_line = int(len(lines[1]) - 1)
    # print("Symbols in line: " + str(simbols_in_line))
    edit_list = list()
    count = 0
    repl_count = 0
    letters = df_uni_commands['new_letter'].tolist()
    old_letters = df_uni_commands['replaced_letter'].tolist()
    real_replased_positions = list()

    for pos in df_uni_commands['position']:
        repl_str = ""
        concrete_line = line_by_pos(pos, simbols_in_line)
        concrete_pos = pos_in_line(pos, simbols_in_line)

        for i in range(len(lines[concrete_line])):
            if i != (concrete_pos - 1):
                repl_str += lines[concrete_line][i]
            else:
                if str(lines[concrete_line][i]) == str(old_letters[count]):
                    repl_str += letters[count]
                    # print(letters[count])
                else:
                    repl_str += lines[concrete_line][i]
                # repl_str += letters[count]

        if (lines[concrete_line]) != repl_str:
            repl_count += 1
            real_replased_positions.append(pos)
        edit_list = list()

        for i in range(len(lines)):
            if i == (concrete_line):
                edit_list.append(repl_str)
                # repl_count += 1
                # print('Line ' + str(i) + ' pos ' + str(concrete_pos) + ': ' + repl_str)
            else:
                edit_list.append(lines[i])

        count += 1
        lines = edit_list

    with open(path_to_result, 'w') as exit_data:
        for line in lines:
            exit_data.write(line)

    print("---> Replaced " + str(repl_count) + " letters")
    print("---> Created new file: " + path_to_result)
    return real_replased_positions


def line_by_pos(position, symbols_in_line):
    return int(int(position) / int(symbols_in_line) + 1)


def pos_in_line(position, symbols_in_line):
    return int(int(position) - (int(int(position) / int(symbols_in_line))) * int(symbols_in_line))


def get_replased_table(path_to_csv):
    df = pd.read_csv(path_to_csv)
    # Выделяем только смену позиций для иматиниба
    right_data = pd.DataFrame(df[[' DRUG_NAME', ' AA_MUTATION']])
    filter_name = right_data[' DRUG_NAME'] == 'Imatinib'
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

    df_uni_commands = pd.DataFrame(unique_comm_table, columns=['position', 'replaced_letter', 'new_letter'])
    # print(df_uni_commands)
    # df_uni_commands[['replaced_letter', 'new_letter']] = df_uni_commands[['replaced_letter', 'new_letter']].astype(str)
    return df_uni_commands


"""

# Вторая часть задания

"""


def get_pept_tables(replased_positions, replased_table, path_to_new_seq):
    result_tables = list()
    lines = list()
    transcrypt = path_to_new_seq[:8]
    AA_mutation = get_mutat_dict(replased_positions, replased_table)
    counter = 0

    int_list = list()
    for pos in replased_positions:
        int_list.append(int(pos))

    replased_positions = sorted(int_list)

    for i in range(15):
        result_tables.append(list())

    with open(path_to_new_seq, 'r') as source_data:
        lines = source_data.readlines()

    symbols_in_line = len(lines[1]) - 1

    sequence_one_str = ""
    for i in range(1, len(lines)):
        sequence_one_str += lines[i][:-1]
    seq_len = len(sequence_one_str)

    for pos in replased_positions:
        left_part = ""
        right_part = ""
        peptide = sequence_one_str[pos - 1]

        if pos <= 15 or (pos + 15) > seq_len:
            continue

        flag = True
        for i in range(1, 16):
            if replased_positions.count(pos + i) == 0 and replased_positions.count(pos - i) == 0:
                None
            else:
                flag = False
                break
        if not flag:
            continue

        counter += 1

        for i in range(1, 16):
            left_part = sequence_one_str[pos - 1 - i]
            right_part = sequence_one_str[pos - 1 + i]
            peptide = left_part + peptide + right_part
            result_tables[i - 1].append(list([peptide, "Immatinib", AA_mutation[str(pos)], transcrypt]))
            # print (result_tables[i-1])

    coverage_percentage = double(31 * counter) / double(seq_len) * 100
    print("---> Coverage percentage: " + str(coverage_percentage))

    return result_tables


def get_mutat_dict(replased_positions, replased_table):
    mutat_dict = dict()

    for pos in replased_positions:
        old_letter = replased_table[replased_table.position == pos]['replaced_letter'].tolist()[0]
        new_letter = replased_table[replased_table.position == pos]['new_letter'].tolist()[0]
        mutat = "p." + old_letter + pos + new_letter
        mutat_dict[pos] = mutat

    # print(mutat_dict)
    return mutat_dict


def peptides_to_csv(result_tables):
    df_list = list()

    for table in result_tables:
        df = pd.DataFrame(table, columns=['PEPTIDE', 'DRUG_NAME', 'AA_MUTATION', 'TRANSCRYPT'])
        df.index.name = 'NUMBERS'
        df_list.append(df)

    for i in range(len(df_list)):
        name = "pept_" + str(i + 1) + ".csv"
        df_list[i].to_csv(name)
        print("---> Create new file: " + str(name))


#csv_file = str(input("Enter path to CSV file...\n"))
#seq_file = str(input("Enter path to sequence file...\n"))
#replased_positions = replased_fun(csv_file, seq_file)
#replased_table = get_replased_table(csv_file)
#new_seq = seq_file[:-4] + '_new' + seq_file[-4:]

#result_table = get_pept_tables(replased_positions, replased_table, new_seq)
#peptides_to_csv(result_table)
