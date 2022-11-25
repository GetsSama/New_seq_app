import PeptideAssembly
import RowDataProcessor as rdp
import Converter as conv

path_abl = "C:\\Users\\Zh_Nikolay\\Desktop\\Sources\\ABL_2.csv"
path_seq = "C:\\Users\\Zh_Nikolay\\Desktop\\Sources\\P00519-2.fasta.txt"
out_dir = "out"
path_config = "config/1000.json"
pept_name_prefix = "peptides_"
path_dop = "Gene_mutations.csv"
transcrypt = "P00519-2"


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


if __name__ == '__main__':
    any_test()
    convert()
