import pandas as pd
import os



input_file = "UNIFI_subset.tsv"
output_file = "/Users/korni/Documents/Dissertation/results/UC_without_placebo_vs_HC_metadata.tsv"
expression_datafile = "clean_UNIFI_biopsy_mRNA_microarray.csv"
geboes_antiTNF_file = "UNIFI_Geboes_AntiTNF.tsv"
out_expression_datafile = "/Users/korni/Documents/Dissertation/results/UC_without_placebo_vs_HC_expression_datafile.tsv"

# Set option to display all columns
pd.set_option('display.max_columns', None)

d1 = pd.read_csv(input_file, sep="\t")
# values_trt = d1["TRT01A"].unique()
print(d1.shape)
# print(values_trt)

# d2 = pd.read_csv(geboes_antiTNF_file, sep="\t")
# print(d2.shape)
#
# print(d1['X'].duplicated().sum())
# print(d2['ids'].duplicated().sum())
#
# # Find the duplicate IDs in d2
# duplicate_ids = d2[d2['ids'].duplicated(keep=False)]  # keep=False to show all duplicates, not just the first occurrence
# print(duplicate_ids)
# merged_df = pd.merge(d1,d2,left_on="X",right_on="ids",how="left")
# print(merged_df.shape)



all_metadata = {}
samples_in_interest = ["Gene.Symbol"]
final_samples = ["Gene.Symbol"]
patient_alt_dict = {}
samples_order_dict = {}
dosages = []

dose = "Ustekinumab 130 mg IV (I-0)"

with open(input_file, 'r') as i:

    i.readline()
    i.readline()

    for line in i:
        line = line.strip().split('\t')
        patient_alt_ID = line[1]
        patient_ID = line[2]
        treatment = line[25]
        period = line[26]
        mayo_score = line[31]

        if "Placebo" in treatment:
            continue
        if "Normal" in treatment:
            condition = "healthy"
        else:
            condition=treatment

        samples_order_dict[patient_ID] = {"condition":condition,"week":period}
        final_samples.append(patient_ID)



        samples_in_interest.append(patient_ID)

print(len(samples_in_interest))

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

df = pd.read_csv(expression_datafile, sep = ",")
df = df[final_samples]
print(df)
df.to_csv(out_expression_datafile, sep = '\t', index = False)

with open(output_file, 'w') as out:

    out.write(f"SampleID\tCondition\tWeek\n")

    for patient in samples_order_dict:
        out.write(patient + '\t' + samples_order_dict[patient]["condition"] + '\t' + samples_order_dict[patient]["week"] + '\n')
