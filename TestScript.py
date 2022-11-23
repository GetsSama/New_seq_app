import RowDataProcessor

path_abl = "ABL.csv"
path_seq = "P00519-2.fasta.txt"

abl = RowDataProcessor.Mutations_data_factory.create_mutations_data_source(path_abl, "ABL")
seq = RowDataProcessor.Sequence_entity(path_seq, "P00519-2")

test_dict = abl.get_drug_and_mutations_dict()

for pair in test_dict.items():
    print(pair)

print("---------------------------------------------------")
oriented_abl = RowDataProcessor.Drug_oriented_mutations(abl)
test_dict2 = oriented_abl.get_drug_and_mutations_dict()

for pair in test_dict2.items():
    print(pair)