import os

import Converter
import peptide_assembly as ps
import Negative_sample as ns


work_dir = "source_data"
seq_main_result = "peptides_source"
seq_negative_result = "negative_sample"
main_pept_result_file_name = "pept_"
negative_sample_file_name = "neg_sam_"

out_sdf_pept = "sdf_peptide_result"
out_sdf_negative_sample = "sdf_negative_result"

if not os.path.exists(work_dir):
    os.mkdir(work_dir)
    os.mkdir(work_dir + "/" + seq_main_result)
    os.mkdir(work_dir + "/" + seq_negative_result)

if not os.path.exists(out_sdf_pept):
    os.mkdir(out_sdf_pept)

if not os.path.exists(out_sdf_negative_sample):
    os.mkdir(out_sdf_negative_sample)


def create_peptides():
    sequence_data = ps.Entry_data()
    while True:
        try:
            abl_path = str(input("Entry path to ABL file...\n"))
            sequence_data.abl_path = abl_path
            break
        except OSError as e:
            print(e)

    while True:
        try:
            seq_path = str(input("Entry path to sequence file...\n"))
            sequence_data.sequence_path = seq_path
            break
        except OSError as e:
            print(e)

    sequence_data.transcrypt_name = str(input("Entry transcrypt name...\n"))
    sequence_data.explorable_drug_name = str(input("Entry drug name...\n"))

    peptides = ps.Peptides(sequence_data)
    peptides.create_peptide_tables_by_entry_data()
    peptides.peptides_to_csv(main_pept_result_file_name, work_dir + "/" + seq_main_result)

    ns.create_negative_sample_tables(sequence_data, negative_sample_file_name, work_dir + "/" + seq_negative_result)

    conv = Converter
    conv.path_to_peptides = work_dir + "/" + seq_main_result
    conv.file_name_prefix = main_pept_result_file_name
    conv.out_dir_name = out_sdf_pept
    conv.start()

    Converter.path_to_peptides = work_dir + "/" + seq_negative_result
    Converter.file_name_prefix = negative_sample_file_name
    Converter.out_dir_name = out_sdf_negative_sample
    Converter.start()


if __name__ == '__main__':
    create_peptides()


