import RowDataProcessor as em
import PeptideAssembly as pept


class __Sequence_and_ABL_data:
    __seq = None
    __abl = None

    def __init__(self, sequence_entity, abl_entity):
        self.__seq = sequence_entity
        self.__abl = abl_entity

    @property
    def sequence_entity(self):
        return self.__seq

    @property
    def abl_entity(self):
        return self.__abl


def create_negative_sample_tables(entry_data, file_name="negative_pept", out_dir="out"):
    all_drugs = em.Row_data_utils.get_unique_values_in_column(entry_data.abl_path, ' DRUG_NAME')
    all_drugs.remove(entry_data.explorable_drug_name)

    drugs_data = dict()

    for drug in all_drugs:
        seq_data = em.Sequence_entity(entry_data.sequence_path, entry_data.transcrypt_name)
        abl_data = em._Mutations_ABL(entry_data.abl_path, drug)
        data_wrap = __Sequence_and_ABL_data(seq_data, abl_data)
        drugs_data[drug] = data_wrap

    table_data_expl_drug = em._Mutations_ABL(entry_data.abl_path, entry_data.explorable_drug_name)
    expl_drug_mutations_list = table_data_expl_drug.get_drug_and_mutations_dict

    for other_drugs_data in drugs_data.values():
        this_abl = other_drugs_data.abl_entity
        this_seq = other_drugs_data.sequence_entity

        this_abl.remove_mutations(expl_drug_mutations_list)
        this_seq.create_replaced_dict(this_abl.get_drug_and_mutations_dict)

    negative_peptides = pept.Peptides()
    for negative_drug in drugs_data.values():
        negative_peptides.convert_to_peptide_tables()
    negative_peptides.peptides_to_csv(file_name, out_dir)