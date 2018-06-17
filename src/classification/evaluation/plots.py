import pandas as pd
import numpy as np
import collections
import matplotlib.pyplot as plt
import sklearn

helices = pd.read_csv("../data/training_data_new.csv")
helices = helices.drop_duplicates()
helices = helices.dropna()

aas = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'H', 'I', 'G', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def sort_dict_by_key(mydict):
    sorted_dict = {}
    for key in sorted(mydict):
           sorted_dict[key] = mydict[key]
    return sorted_dict

helices = pd.read_csv("../data/training_data_new.csv")
helices = helices.drop_duplicates()
helices = helices.dropna()
nontm_helices = helices.loc[helices["Is Transmembrane"] == -1.0]
tm_helices = helices[helices["Is Transmembrane"] == 1.0]
nontm_protein_helices = pd.read_csv("../data/helices_in_non_tm_proteins.csv")
lengths = []
for index, row in nontm_protein_helices.iterrows():
    lengths.append(row["Helix Length"] + 1)

lengths = np.array(lengths)
nontm_protein_helices['Helix Length'] = lengths
nontm_protein_helices = nontm_protein_helices[nontm_protein_helices['Helix Length'] > 4]
nontm_protein_helices = nontm_protein_helices.replace(0, -1.0)
nontm_protein_helices = nontm_protein_helices.drop_duplicates()
nontm_protein_helices = nontm_protein_helices.dropna()

frames = [nontm_helices, nontm_protein_helices]
negative_samples = pd.concat(frames)

print("TM Helices:",tm_helices.shape[0])
print("Non-TM helices",nontm_helices.shape[0])
print("Non TM protein helices:",nontm_protein_helices.shape[0])

# Plot positive vs negative samples
f1 = plt.figure()
x = np.arange(2)
tm, nontm = plt.bar(x, [tm_helices.shape[0], nontm_helices.shape[0]+ nontm_protein_helices.shape[0]])
tm.set_color('r')
plt.xticks(x, ('Transmembrane', 'Non-Transmembrane'))
plt.ylim(0,60000)
plt.savefig("../figures/samples.pdf")
plt.savefig("../figures/samples.png")
plt.show()

# Plot amino acid distribution of transmembrane helices
tm_helix_sequences = list(tm_helices["Helix Sequence"])
tm_helix_aas = dict(collections.Counter("".join(tm_helix_sequences)))
tm_helix_aas_count = sum(tm_helix_aas.values())
for key in tm_helix_aas.keys():
    tm_helix_aas[key] = tm_helix_aas.get(key) / tm_helix_aas_count * 100
tm_helix_aas = sort_dict_by_key(tm_helix_aas)

f2 = plt.figure()
plt.bar(range(len(tm_helix_aas)), list(tm_helix_aas.values()), align='center', color='r')
plt.xticks(range(len(tm_helix_aas)), list(tm_helix_aas.keys()))
plt.title("TM helices Amino Acid Frequency")
plt.xlabel("Amino Acids")
plt.ylabel("Frequency (%)")
plt.ylim([0, 20])
plt.savefig('../figures/tm_aa_distribution.png')
plt.savefig('../figures/tm_aa_distribution.pdf')
plt.show()

# Plot amino acid distribution of non-transmembrane helices
negative_sequences = list(negative_samples["Helix Sequence"])
negative_aas = dict(collections.Counter("".join(negative_sequences)))
negative_aas_count = sum(negative_aas.values())
for key in negative_aas.keys():
    negative_aas[key] = negative_aas.get(key) / negative_aas_count * 100
negative_aas = sort_dict_by_key(negative_aas)

f3 = plt.figure()
plt.bar(range(len(negative_aas)), list(negative_aas.values()), align='center')
plt.xticks(range(len(negative_aas)), list(negative_aas.keys()))
plt.title("Non-TM helices Amino Acid Frequency")
plt.xlabel("Amino Acids")
plt.ylabel("Frequency (%)")
plt.ylim([0, 20])
plt.savefig('../figures/nontm_aa_distribution.png')
plt.savefig('../figures/nontm_aa_distribution.pdf')
plt.show()

# Plot comparison of amino acid distribution of helices
diff_dict = {}
for key in tm_helix_aas:
    diff_dict[key] = tm_helix_aas[key] - negative_aas[key]

dic = {'R': '#ad0020', 'L': '#204cfd', 'N': '#fd3116', 'D': '#fd3116', 'Q': '#fd3116', 'E': '#fd3116', 'H': '#ff4a26', 'P': '#ffa07b', 'Y': '#ffac87', 'W': '#ffb895', 'S': '#fabc9a', 'T': '#f8be9f', 'G': '#e9c9b5', 'A': '#89b0ff', 'M': '#85adff', 'C': '#6d96ff', 'F': '#5d8bff', 'V': '#1639ed', 'I': '#1733dd'}
dic = {"D":'#ffa07b',
"E":'#d5cfcd',
"N":'#ff885f',
"Q":'#ff9f79',
"K":'#ac0020',
"R":'#ff784f',
"H":'#f31d0e',
"G":'#dbcec7',
"P":'#2855ff',
"S":'#e2ccbe',
"T":'#deccc2',
"C":'#97bcff',
"M":'#5e8cff',
"A":'#86aeff',
"V":'#5b89ff',
"I":'#2150ff',
"L":'#1a42f3',
"F":'#1732d9',
"W":'#739eff',
"Y":'#4571ff'}
dic = sort_dict_by_key(dic)

aas = []
for key in dic.keys():
    aas.append(key)

print(aas)
print(dic)

f4 = plt.figure()
barlist = plt.bar(range(len(negative_aas)), list(diff_dict.values()), align='center')
for i in range(0, len(barlist)-1):
    barlist[i].set_color(dic.get(aas[i]))
plt.xticks(range(len(tm_helix_aas)), list(tm_helix_aas.keys()))
plt.title("Amino Acid Differences of Transmembrane Helices")
plt.xlabel("Amino Acids")
plt.ylabel("Frequency Diff (%)")
plt.ylim([-7.5, 7.5])
plt.savefig('figures/comp_aa_distribution.png')
plt.savefig('figures/comp_aa_distribution.pdf')
plt.show()

# Plot Helix Length
tm_helix_length = list(tm_helices["Helix Length"])
neg_helix_length = list(negative_samples["Helix Length"])
f5 = plt.figure()
plt.hist(tm_helix_length, bins=range(5,50), color='red', density=True)
plt.hist(neg_helix_length, bins=range(5,50), alpha=0.5, density=True)
plt.savefig('../figures/helix_length.png')
plt.savefig('../figures/helix_length.pdf')
plt.show()

f6 = plt.figure()
plt.bar(range(len(negative_aas)), list(negative_aas.values()), align='center')
plt.bar(range(len(tm_helix_aas)), list(tm_helix_aas.values()), align='center', color='r', alpha=0.5)

plt.xticks(range(len(negative_aas)), list(negative_aas.keys()))
plt.title("Non-TM helices Amino Acid Frequency")
plt.xlabel("Amino Acids")
plt.ylabel("Frequency (%)")
plt.ylim([0, 20])
plt.savefig('../figures/nontm_aa_distribution.png')
plt.savefig('../figures/nontm_aa_distribution.pdf')
plt.show()