import peptide_assembly

print("Please enter the present configuration:")
seq = str(input("Enter the path to sequence file...\n"))
abl = str(input("Enter the path to gene file...\n"))
transcrypt = str(input("Enter the transcrypt name...\n"))
drug = str(input("Enter the drug name...\n"))

peptides = peptide_assembly.Peptides(seq, abl, transcrypt, drug)
peptides.create_peptide_tables()
peptides.peptides_to_csv()
