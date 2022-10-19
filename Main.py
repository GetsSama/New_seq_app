import peptide_assembly
import EntryModule as em
import pandas as pd

# print("Please enter the present configuration:")
# seq = str(input("Enter the path to sequence file...\n"))
# abl = str(input("Enter the path to gene file...\n"))
# transcrypt = str(input("Enter the transcrypt name...\n"))
# drug = str(input("Enter the drug name...\n"))
seq = 'P00519-2.fasta.txt'
abl = 'ABL.csv'
transcrypt = 'P00519-2'
drug = 'Imatinib'

peptides = peptide_assembly.Peptides(seq, abl, transcrypt, drug)
peptides.create_peptide_tables()
peptides.peptides_to_csv()

# df = em.Sequence_tools.get_replased_table(abl, drug)
# replace = list()
# # for row in df.iterrows():
# #     row_ser = row[1]
# #     print(row_ser['replaced_letter'] + row_ser['position'] + row_ser['new_letter'])
#
# entity = em.ABL_entity(abl, drug)
# replace = entity.get_mutations()
# seq_entity = em.Sequence_entity(seq, transcrypt)
#
# for pair in replace:
#     print(pair)
#
# seq_entity.create_replaced_dict(entity.replaced_positions, entity.replaced_letters, replace)
# replace = seq_entity._replaced_dict
#
# for pair in replace.items():
#     print(pair)