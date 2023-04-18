import PeptideAssembly
import RowDataProcessor as rdp
import Converter as conv
import os
import pandas as pd

path_abl = "C:\\Users\\Zh_Nikolay\\Desktop\\Sources\\ABL1_1.csv"
path_seq = "C:\\Users\\Zh_Nikolay\\Desktop\\Sources\\P00519-1.fasta.txt"
out_dir = "out"
path_config = "config/1000.json"
pept_name_prefix = "peptides_"
path_dop = "C:\\Users\\Zh_Nikolay\\Desktop\\Sources\\ABL_Gene_mutations.csv"
transcrypt = "ABL"


def peptides():
    seq_source = rdp.Sequence_entity(path_seq, transcrypt)
    mutate_source = rdp.Drug_oriented_mutations(rdp.Mutations_data_factory.create_mutations_data_source(path_abl))
    drug_and_mut_dict = mutate_source.get_drug_and_mutations_dict()

    pept_creator = PeptideAssembly.Peptides_creator()
    pept_creator.convert_to_peptide_tables(seq_source, drug_and_mut_dict)
    pept_creator.peptides_to_csv(pept_name_prefix)


def convert():
    converter = conv.SeqToSDF_manager()
    converter.set_manager_properties(path_config, out_dir, pept_name_prefix)
    converter.start()


def any_test():
    mutate_source = rdp.Drug_oriented_mutations(
        rdp.Mutations_data_factory.create_mutations_data_source(path_dop, "AnyTable", "AA Mutation"))
    drug_and_mut_dict = mutate_source.get_drug_and_mutations_dict()

    mutate_source2 = rdp.Drug_oriented_mutations(rdp.Mutations_data_factory.create_mutations_data_source(path_abl))
    drug_and_mut_dict2 = mutate_source2.get_drug_and_mutations_dict()

    result_dict = rdp.Row_data_utils.concatenate_drug_and_mutations_dicts(drug_and_mut_dict, drug_and_mut_dict2)
    for pair in result_dict.items():
        print(pair)

    seq_source = rdp.Sequence_entity(path_seq, transcrypt)
    pept_creator = PeptideAssembly.Peptides_creator()
    pept_creator.convert_to_peptide_tables(seq_source, result_dict)
    pept_creator.peptides_to_csv(pept_name_prefix)

def get_differences_on_mutation(path_to_comparable_csv: str, mutations_column_name: str):
    mutate_source = rdp.Drug_oriented_mutations(
        rdp.Mutations_data_factory.create_mutations_data_source(path_dop, "AnyTable", "AA Mutation"))
    drug_and_mut_dict = mutate_source.get_drug_and_mutations_dict()

    mutate_source2 = rdp.Drug_oriented_mutations(rdp.Mutations_data_factory.create_mutations_data_source(path_abl))
    drug_and_mut_dict2 = mutate_source2.get_drug_and_mutations_dict()

    result_dict = rdp.Row_data_utils.concatenate_drug_and_mutations_dicts(drug_and_mut_dict, drug_and_mut_dict2)
    all_mutations = set()

    for key, value in result_dict.items():
        for mut in value:
            all_mutations.add(mut)

    print(all_mutations)
    print(len(all_mutations))

    df = pd.read_csv(path_to_comparable_csv, sep=';')
    print(df.keys())
    other_mutat_table = pd.DataFrame(df[[mutations_column_name]])
    other_mutat_table = other_mutat_table[mutations_column_name]
    other_mutat_table = other_mutat_table.tolist()

    print(other_mutat_table)

    differences_set = set(other_mutat_table).difference(all_mutations)
    print(differences_set)
    drug_and_differences = dict()
    drug_and_differences["No drug"] = differences_set

    seq_source = rdp.Sequence_entity(path_seq, transcrypt)
    pept_creator = PeptideAssembly.Peptides_creator()
    pept_creator.convert_to_peptide_tables(seq_source, drug_and_differences)
    pept_creator.peptides_to_csv(pept_name_prefix)

    df_out = pd.DataFrame(differences_set, columns=['DIFF_MUTATIONS'])
    df_out.index.name = 'NUMBERS'
    df_out.to_csv(out_dir + "\\differences.csv")

def create_olmpass_command_files(gene_name: str):
    task_file_content = ""
    abs_out_path = os.path.abspath(out_dir)
    for i in range(3, 16):
        for j in range(1, 16):
            first_row = "BaseCreate=" + str(i) + ";" + gene_name + "_descript_" + str(i) + "_pept_len_" + str(j) + "\n"
            second_row = ""
            second_row += "BaseAddNewData=" + pept_name_prefix + str(j) + "_sdf.sdf;DRUG_NAME\n"
            other_rows = "BaseSave\nBaseTraining\nBaseValidation\nBaseClose\n"
            task_file_content += first_row + second_row + other_rows + "\n"

        # if not os.path.exists(task_file_path):
        #     os.mkdir(task_file_path)
    task_file_path = abs_out_path + "\\task" + ".txt"
    with open(task_file_path, "w") as task_file:
        task_file.write(task_file_content)
    #print(task_file_content)


if __name__ == '__main__':
    #any_test()
    get_differences_on_mutation("C:\\Users\\Zh_Nikolay\\Desktop\\Sources\\GnomAD_extra.csv", "Ptein")
    #peptides()
    convert()
    create_olmpass_command_files("ABL")
