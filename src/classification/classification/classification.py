from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn.metrics import accuracy_score, matthews_corrcoef, f1_score, classification_report, confusion_matrix

import pandas as pd
import numpy as np
from sklearn.externals import joblib

# Load training data

helices = pd.read_csv("../data/training_data_new.csv")
helices = helices.drop_duplicates()
helices = helices.dropna()
nontm_helices = helices.loc[helices["Is Transmembrane"] == -1.0]
print(nontm_helices.shape)
helices = helices[helices["Is Transmembrane"] != -1.0]
print(helices.shape)
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
print(nontm_protein_helices.shape)

print("Positive Samples:", str(helices.shape[0]))
print("Negative Samples:", str(nontm_helices.shape[0]+nontm_protein_helices.shape[0]))
print("Overall Samples:", str(helices.shape[0]+nontm_protein_helices.shape[0]+nontm_helices.shape[0]))



#full_data.to_csv("full_data.csv", sep=',', encoding='utf-8', index=False)

#TODO Add non-tm protein helices to training data

aas = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'H', 'I', 'G', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def count_aas(helix):
    try:
        counts = []
        for aa in aas:
            counts.append(helix.count(aa))
    except:
        print(helix)
    return counts


def create_aa_frequency_features(df):
    X = []
    Y = list(df["Is Transmembrane"])
    helices = list(df["Helix Sequence"])
    for helix in helices:
        X.append(count_aas(helix))

    return np.array(X), np.array(Y)


X_tm, y_tm = create_aa_frequency_features(helices)
X_nontm_protein_helices, y_nontm_protein_helices = create_aa_frequency_features(nontm_protein_helices)
X_nontm, y_nontm = create_aa_frequency_features(nontm_helices)

X = np.concatenate((X_tm, X_nontm_protein_helices, X_nontm), axis=0)
y = np.concatenate((y_tm, y_nontm_protein_helices, y_nontm), axis=0)

print("Split data into train and test")
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3)

print("Linear SVM")
svc = svm.SVC(kernel="linear", class_weight="balanced")
svc.fit(X_train, y_train)

print("Predict test data")
y_pred = svc.predict(X_test)
tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()

print("Results:")
print("Accuracy:", accuracy_score(y_test, y_pred))
print("TN:",tn,"TP:",tp,"FN:",fn,"FP:",fp)
print("MCC:", matthews_corrcoef(y_test, y_pred))
print("F1-Score:", f1_score(y_test, y_pred))
print(classification_report(y_test, y_pred))

print(svc.coef_[0])

print("RBF SVM")
svc = svm.SVC(kernel="rbf", class_weight="balanced", C=100, gamma=0.1)
svc.fit(X_train, y_train)
#joblib.dump(svc, '../tm_helix_predictor/svm_c100_g01.pkl')
print("Predict test data")
y_pred = svc.predict(X_test)
tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()

print("Results:")
print("Accuracy:", accuracy_score(y_test, y_pred))
print("TN:",tn,"TP:",tp,"FN:",fn,"FP:",fp)
print("MCC:", matthews_corrcoef(y_test, y_pred))
print("F1-Score:", f1_score(y_test, y_pred))
print(classification_report(y_test, y_pred))

#print("RBF optimized SVM")
#svc = svm.SVC(kernel="rbf", C=100, gamma=0.1,  class_weight="balanced")
#svc.fit(X_train, y_train)

#print("Predict test data")
#y_pred = svc.predict(X_test)
#tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()

#print("Results:")
#print("Accuracy:", accuracy_score(y_test, y_pred))
#print("TN:",tn,"TP:",tp,"FN:",fn,"FP:",fp)
#print("MCC:", matthews_corrcoef(y_test, y_pred))
#print("F1-Score:", f1_score(y_test, y_pred))
#print(classification_report(y_test, y_pred))