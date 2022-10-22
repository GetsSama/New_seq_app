import Negative_sample
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

entry_data_Imatinib = peptide_assembly.Entry_data()
entry_data_Imatinib.abl_path = abl
entry_data_Imatinib.sequence_path = seq
entry_data_Imatinib.transcrypt_name = transcrypt
entry_data_Imatinib.explorable_drug_name = drug

# pept = peptide_assembly.Peptides(entry_data_Imatinib)
# pept.create_peptide_tables_by_entry_data()
# pept.peptides_to_csv()

Negative_sample.create_negative_sample_tables(entry_data_Imatinib)