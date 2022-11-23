import RowDataProcessor

path_abl = "ABL.csv"
path_seq = "P00519-2.fasta.txt"

abl = RowDataProcessor.Mutations_data_factory.create_mutations_data_source(path_abl, "ABL")
seq = RowDataProcessor.Sequence_entity(path_seq, "P00519-2")

test_dict = abl.get_drug_and_mutations_dict()

for pair in test_dict.items():
    print(pair)

freq_drug = RowDataProcessor.Drug_oriented_mutations(abl).most_freq_drug()
print(freq_drug)