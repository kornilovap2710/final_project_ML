import pandas as pd
import os
import re

print(os.getcwd())

input_file = "../data/UNIFI_for_Polina/UNIFI_subset.tsv"
# output_file = "/Users/korni/Documents/Dissertation/results/UC_without_placebo_vs_HC_metadata.tsv"
expression_datafile = "../data/UNIFI_for_Polina/clean_UNIFI_biopsy_mRNA_microarray.csv"
# geboes_antiTNF_file = "UNIFI_Geboes_AntiTNF.tsv"
# out_expression_datafile = "/Users/korni/Documents/Dissertation/results/UC_without_placebo_vs_HC_expression_datafile.tsv"
output_dir = "../data/UNIFI_for_Polina/filtered_data"
# Set option to display all columns
pd.set_option('display.max_columns', None)


# all_metadata = {}
# samples_in_interest = ["Gene.Symbol"]
# final_samples = ["Gene.Symbol"]
# patient_alt_dict = {}
# samples_order_dict = {}
# dosages = []

doses = ["Ustekinumab 130 mg IV (I-0)", "Ustekinumab 6 mg/kg (390 mg) (I-0)", "Ustekinumab 6 mg/kg (520 mg) (I-0)","Ustekinumab 6 mg/kg (260 mg) (I-0)"]


def filter_metadata(input_file, dose, include_healthy=True):
    """Filters the metadata based on the given dose and optionally includes healthy volunteers."""
    samples_order_dict = {}
    final_samples = ["Gene.Symbol"]
    sample_weeks = {}

    with open(input_file, 'r') as i:
        i.readline()
        i.readline()

        for line in i:
            line = line.strip().split('\t')
            patient_ID = line[2]
            treatment = line[25]
            period = line[26]

            if "Placebo" in treatment:
                continue
            if dose in treatment:
                condition = dose
            elif include_healthy==True and "Normal" in treatment:
                condition = "Healthy"
            else:
                continue
            samples_order_dict[patient_ID] = (condition,period)
            final_samples.append(patient_ID)
            sample_weeks[patient_ID] = period

    return samples_order_dict, final_samples, sample_weeks

def filter_expression_data(expression_datafile, final_samples, output_expression_file):
    """Filters the expression data to include only relevant samples."""
    df = pd.read_csv(expression_datafile, sep=",")
    df = df[final_samples]
    df.to_csv(output_expression_file, sep='\t', index=False)
    print(f"Expression data saved to {output_expression_file}")

def save_metadata(samples_order_dict, output_metadata_file):
    """Saves the metadata file with sample conditions."""
    with open(output_metadata_file, 'w') as out:
        out.write("SampleID\tCondition\tWeek\n")
        for patient, (condition,period) in samples_order_dict.items():
            out.write(f"{patient}\t{condition}\t{period}\n")
    print(f"Metadata saved to {output_metadata_file}")


def process_dose(input_file, expression_datafile, dose, output_dir, include_healthy=True):
    """Processes the given dose, generating metadata and expression files, optionally including healthy volunteers."""
    suffix = "_with_HC" if include_healthy else ""
    dose_cleaned = re.sub(r"[ /()-]","_",dose)

    output_metadata_file = f"{output_dir}/metadata_{dose_cleaned}{suffix}.tsv"
    output_expression_file = f"{output_dir}/expression_{dose_cleaned}{suffix}.tsv"

    samples_order_dict, final_samples, _ = filter_metadata(input_file, dose, include_healthy)
    filter_expression_data(expression_datafile, final_samples, output_expression_file)
    save_metadata(samples_order_dict, output_metadata_file)


def generate_all_dose_files(input_file, expression_datafile, doses, output_dir):
    """Generates metadata and expression data files for each dose with and without healthy volunteers."""
    for dose in doses:
        process_dose(input_file, expression_datafile, dose, output_dir, include_healthy=True)
        process_dose(input_file, expression_datafile, dose, output_dir, include_healthy=False)

generate_all_dose_files(input_file, expression_datafile, doses, output_dir)

# with open(input_file, 'r') as i:
#
#     i.readline()
#     i.readline()
#
#     for line in i:
#         line = line.strip().split('\t')
#         patient_alt_ID = line[1]
#         patient_ID = line[2]
#         treatment = line[25]
#         period = line[26]
#         mayo_score = line[31]
#
#         if "Placebo" in treatment:
#             continue
#         if "Normal" in treatment:
#             condition = "Healthy"
#         else:
#             condition=treatment
#
#         samples_order_dict[patient_ID] = {"condition":condition,"week":period}
#         final_samples.append(patient_ID)
#
#
#
#         samples_in_interest.append(patient_ID)
#
# print(len(samples_in_interest))

#         if patient_alt_ID not in patient_alt_dict:
#             patient_alt_dict[patient_alt_ID] = []

#         patient_alt_dict[patient_alt_ID].append(patient_ID)

# for pat in patient_alt_dict:

#     if len(patient_alt_dict[pat]) == 2:

#         for sampleid in patient_alt_dict[pat]:
#             samples_in_interest.append(sampleid)


# all_metadata_samples = ["Gene.Symbol"]

# with open(geboes_antiTNF_file, 'r') as geboes:

#     geboes.readline()

#     for line in geboes:
#         line = line.strip().split('\t')
#         sampleID = line[1]
#         antiTNF = line[2]
#         geboes_total = line[3]

#         if sampleID not in all_metadata:
#             all_metadata[sampleID] = [antiTNF, geboes_total]

#         if sampleID in samples_in_interest:
#             all_metadata_samples.append(sampleID)
#
# final_samples = ["Gene.Symbol"]
#
# with open(input_file, 'r') as i:
#
#     i.readline()
#     i.readline()
#
#     for line in i:
#         line = line.strip().split('\t')
#         patient_ID = line[2]
#         treatment = line[25]
#         period = line[26]
#         mayo_score = line[31]
#
#         if patient_ID in samples_in_interest:
#             # actual_antiTNF = all_metadata[patient_ID][0]
#             # geboes_score = all_metadata[patient_ID][1]
#
#             # if actual_antiTNF != "Normal" and geboes_score == "NA":
#             #     continue
#
#             # if geboes_score != "NA":
#             #     geboes_score = int(geboes_score)
#
#             #     if geboes_score > 0 and geboes_score <= 5:
#             #         actual_geboes_score = "Low"
#
#             #     if geboes_score >= 6 and geboes_score <= 10:
#             #         actual_geboes_score = "Medium"
#
#             #     if geboes_score >= 11 and geboes_score <= 15:
#             #         actual_geboes_score = "High"
#
#             #     if geboes_score >= 16 and geboes_score <= 20:
#             #         actual_geboes_score = "Ultra_high"
#
#             # if geboes_score == "NA":
#             #     actual_geboes_score = "NA"
#
#             # if mayo_score != "NA":
#             #     mayo_score = int(mayo_score)
#
#             #     # if mayo_score <= 1:
#             #     #     actual_mayo_score = "0-1"
#
#             #     # if mayo_score >= 2 and mayo_score <= 3:
#             #     #     actual_mayo_score = "2-3"
#
#             #     if mayo_score == 5:
#             #         actual_mayo_score = "5"
#
#             #     if mayo_score >= 6 and mayo_score <= 7:
#             #         actual_mayo_score = "6-7"
#
#             #     if mayo_score >= 8 and mayo_score <= 9:
#             #         actual_mayo_score = "8-9"
#
#             #     if mayo_score == 10:
#             #         actual_mayo_score = "10"
#
#             #     if mayo_score == 11:
#             #         actual_mayo_score = "11"
#
#             #     if mayo_score == 12:
#             #         actual_mayo_score = "12"
#
#             if "Ustekinumab" in treatment:
#                 condition = "UC"
#
#             # if "Ustekinumab 6 mg/kg (390 mg)" in treatment:
#             #     condition = "UC_6_390"
#
#             # if "Ustekinumab 6 mg/kg (260 mg)" in treatment:
#             #     condition = "UC_6_260"
#
#             # if "Ustekinumab 6 mg/kg (520 mg)" in treatment:
#             #     condition = "UC_6_520"
#
#             if "Normal" in treatment:
#                 condition = "HC"
#
#             # new_info = all_metadata[patient_ID]
#             # new_info.append(condition)
#             samples_order_dict[patient_ID] = condition
#             final_samples.append(patient_ID)
#
# print(len(final_samples))

# df = pd.read_csv(expression_datafile, sep = ",")
# df = df[final_samples]
# print(df)
# df.to_csv(out_expression_datafile, sep = '\t', index = False)
#
# with open(output_file, 'w') as out:
#
#     out.write(f"SampleID\tCondition\tWeek\n")
#
#     for patient in samples_order_dict:
#         out.write(patient + '\t' + samples_order_dict[patient]["condition"] + '\t' + samples_order_dict[patient]["week"] + '\n')
