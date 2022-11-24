import PeptideAssembly
import RowDataProcessor as rdp
import Converter as conv

path_abl = "ABL.csv"
path_seq = "P00519-2.fasta.txt"
out_dir = "out"
path_config = "config/1000.json"
pept_name_prefix = "peptides_"

def peptides():
    seq_source = rdp.Sequence_entity(path_seq, "P00519-2")
    mutate_source = rdp.Drug_oriented_mutations(rdp.Mutations_data_factory.create_mutations_data_source(path_abl))

    pept_creator = PeptideAssembly.Peptides_creator()
    pept_creator.convert_to_peptide_tables(seq_source, mutate_source)
    pept_creator.peptides_to_csv(pept_name_prefix)

def convert():
    converter = conv.SeqToSDF_manager()
    converter.set_manager_properties(path_config, out_dir, pept_name_prefix)
    converter.start()


if __name__ == '__main__':
    peptides()
    convert()