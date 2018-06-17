import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.svm import SVC

# Load training data

helices = pd.read_csv("../data/training_data_new.csv")
helices = helices.drop_duplicates()
helices = helices.dropna()
nontm_helices = helices.loc[helices["Is Transmembrane"] == -1.0]
print(nontm_helices.shape)
helices = helices[helices["Is Transmembrane"] != -1.0]
print(helices.shape)
nontm_protein_helices = pd.read_csv("../data/helices_in_non_tm_proteins.csv")
nontm_protein_helices = nontm_protein_helices.replace(0, -1.0)
nontm_protein_helices = nontm_protein_helices.drop_duplicates()
nontm_protein_helices = nontm_protein_helices.dropna()
print(nontm_protein_helices.shape)



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

# Split the dataset in two equal parts
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.5, random_state=0)

# Set the parameters by cross-validation
tuned_parameters = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4],'C': [1, 10, 100, 1000], 'class_weight':['balanced', None]},
                    {'kernel': ['linear'], 'C': [1, 10, 100, 1000], 'class_weight':['balanced', None]}]

scores = ['f1', 'accuracy']

for score in scores:
    print("# Tuning hyper-parameters for %s" % score)
    print()

    clf = GridSearchCV(SVC(), tuned_parameters, cv=5,
                       scoring=score)
    clf.fit(X_train, y_train)

    print("Best parameters set found on development set:")
    print()
    print(clf.best_params_)
    print()
    print("Grid scores on development set:")
    print()
    means = clf.cv_results_['mean_test_score']
    stds = clf.cv_results_['std_test_score']
    for mean, std, params in zip(means, stds, clf.cv_results_['params']):
        print("%0.3f (+/-%0.03f) for %r"
              % (mean, std * 2, params))
    print()

    print("Detailed classification report:")
    print()
    print("The model is trained on the full development set.")
    print("The scores are computed on the full evaluation set.")
    print()
    y_true, y_pred = y_test, clf.predict(X_test)
    print(classification_report(y_true, y_pred))
    print()