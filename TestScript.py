import RowDataProcessor

path_abl = "ABL.csv"
path_seq = "P00519-2.fasta.txt"

abl = RowDataProcessor.Mutations_ABL(path_abl)
seq = RowDataProcessor.Sequence_entity(path_seq, "P00519-2")

seq.create_replaced_dict(abl.replaced_positions, abl.new_letters, abl.get_mutations)
dict1 = seq.get_mutation_and_sequence_mapping
dict2 = seq.alternative(abl.get_mutations)

print(dict1 == dict2)

