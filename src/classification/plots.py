import pandas as pd
import collections
import matplotlib.pyplot as plt

helices = pd.read_csv("training_data_new.csv")
helices = helices.drop_duplicates()
helices = helices.dropna()


aas = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'H', 'I', 'G', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def sort_dict_by_key(mydict):
    sorted_dict = {}
    for key in sorted(mydict):
           sorted_dict[key] = mydict[key]
    return sorted_dict

nontm_protein_helices = pd.read_csv("helices_in_non_tm_proteins.csv")
nontm_protein_helices = nontm_protein_helices.replace(0, -1.0)
nontm_protein_helices = nontm_protein_helices.drop_duplicates()
nontm_protein_helices = nontm_protein_helices.dropna()
tm_helices = helices.loc[helices['Is Transmembrane'] == 1]
nontm_helices = helices.loc[helices['Is Transmembrane'] == -1]
partly_tm_helices = helices.loc[helices['Is Transmembrane'] == 0]
print("TM Helices:",tm_helices.shape[0])
print("Non-TM helices",nontm_helices.shape[0])
print("Partly TM helices:",partly_tm_helices.shape[0])
print("Non TM protein helices:",nontm_protein_helices.shape[0])


tm_helix_sequences = list(tm_helices["Helix Sequence"])
tm_helix_aas = dict(collections.Counter("".join(tm_helix_sequences)))
tm_helix_aas_count = sum(tm_helix_aas.values())
for key in tm_helix_aas.keys():
    tm_helix_aas[key] = tm_helix_aas.get(key) / tm_helix_aas_count * 100
tm_helix_aas = sort_dict_by_key(tm_helix_aas)

nontm_helix_sequences = list(nontm_helices["Helix Sequence"])
nontm_helix_aas = dict(collections.Counter("".join(nontm_helix_sequences)))
nontm_helix_aas_count = sum(nontm_helix_aas.values())
for key in nontm_helix_aas.keys():
    nontm_helix_aas[key] = nontm_helix_aas.get(key) / nontm_helix_aas_count * 100
nontm_helix_aas = sort_dict_by_key(nontm_helix_aas)

nontm_protein_helix_sequences = list(nontm_protein_helices["Helix Sequence"])
nontm_protein_helix_aas = dict(collections.Counter("".join(nontm_protein_helix_sequences)))
nontm_protein_helix_aas_count = sum(nontm_protein_helix_aas.values())
for key in nontm_protein_helix_aas.keys():
    nontm_protein_helix_aas[key] = nontm_protein_helix_aas.get(key) / nontm_protein_helix_aas_count * 100
nontm_protein_helix_aas = sort_dict_by_key(nontm_protein_helix_aas)

print("There are", len(tm_helix_sequences), "Transmembrane Helices vs", len(nontm_helix_sequences), "Non-TM Helices")
f1 = plt.figure()
plt.bar(range(len(nontm_helix_aas)), list(nontm_helix_aas.values()), align='center', color='blue', alpha=0.5)
plt.bar(range(len(tm_helix_aas)), list(tm_helix_aas.values()), align='center', color='red', alpha=0.5)
plt.xticks(range(len(tm_helix_aas)), list(tm_helix_aas.keys()))
plt.title("Comparison of Amino Acid Frequency")
plt.legend(["Helices", "Transmembrane Helices"])
plt.xlabel("Amino Acids")
plt.ylabel("Frequency (%)")
plt.ylim([0, 20])
plt.savefig('freq.png')
plt.show

diff_dict = {}
for key in tm_helix_aas:
    diff_dict[key] = tm_helix_aas[key] - nontm_helix_aas[key]

f2 = plt.figure()
plt.bar(range(len(nontm_protein_helix_aas)), list(diff_dict.values()), align='center', color='blue')
plt.xticks(range(len(tm_helix_aas)), list(tm_helix_aas.keys()))
plt.title("Comparison of Amino Acid Frequency")
plt.xlabel("Amino Acids")
plt.ylabel("Frequency Diff (%)")
plt.ylim([-4, 4])
plt.savefig('freqdiff.png')
plt.show

print("There are", len(tm_helix_sequences), "Transmembrane Helices vs", len(nontm_protein_helix_sequences), "Non-TM Protein Helices")
f3 = plt.figure()
plt.bar(range(len(nontm_protein_helix_aas)), list(nontm_protein_helix_aas.values()), align='center', color='blue', alpha=0.5)
plt.bar(range(len(tm_helix_aas)), list(tm_helix_aas.values()), align='center', color='red', alpha=0.5)
plt.xticks(range(len(tm_helix_aas)), list(tm_helix_aas.keys()))
plt.title("Comparison of Amino Acid Frequency")
plt.legend(["Helices", "Transmembrane Helices"])
plt.xlabel("Amino Acids")
plt.ylabel("Frequency (%)")
plt.ylim([0, 20])
plt.savefig('freq_full.png')
plt.show

diff2_dict = {}
for key in tm_helix_aas:
    diff2_dict[key] = tm_helix_aas[key] - nontm_protein_helix_aas[key]

f4 = plt.figure()
plt.bar(range(len(nontm_protein_helix_aas)), list(diff2_dict.values()), align='center', color='blue')
plt.xticks(range(len(tm_helix_aas)), list(tm_helix_aas.keys()))
plt.title("Comparison of Amino Acid Frequency")
plt.xlabel("Amino Acids")
plt.ylabel("Frequency Diff (%)")
plt.ylim([-10, 10])
plt.savefig('freqdiff_full.png')
plt.show