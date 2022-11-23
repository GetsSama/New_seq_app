import PeptideAssembly
import RowDataProcessor as rdp

path_abl = "ABL.csv"
path_seq = "P00519-2.fasta.txt"
out_dir = "out"

seq_source = rdp.Sequence_entity(path_seq, "P00519-2")
mutate_source = rdp.Drug_oriented_mutations(rdp.Mutations_data_factory.create_mutations_data_source(path_abl))

pept_creator = PeptideAssembly.Peptides_creator()
pept_creator.convert_to_peptide_tables(seq_source, mutate_source)
pept_creator.peptides_to_csv("peptides_")