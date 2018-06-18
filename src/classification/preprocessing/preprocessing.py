import xml.etree.ElementTree as ET
import os
import pandas as pd

class PDBTM:
    def __init__(self, id=str, chains=[]):
        self.id = id
        self.chains = chains


class Chain:
    def __init__(self, id=str, num_tm=int, typ=str, seq=str, regions=[]):
        self.id = id
        self.num_tm = num_tm
        self.typ = typ
        self.seq = seq
        self.regions = regions


class Region:
    def __init__(self, seq_beg=int, pdb_beg=int, seq_end=int, pdb_end=int, typ=str):
        self.seq_beg = seq_beg
        self.pdb_beg = pdb_beg
        self.seq_end = seq_end
        self.pdb_end = pdb_end
        self.typ = typ


def parse_pdbtm_xmls(paths):
    '''
    Parses PDBTM XML files to a list od PDBTM Objects.
    A PDBTM Object contains the PDB ID as well as a list of Chain objects.
    A Chain Object contains the chain ID, the tm number, the type, the sequence and a list of Region objects.
    A Region object contains the begin and end indices of the sequence and the pdb file, as well as the type.
    '''
    ns = {'pdbtm': 'http://pdbtm.enzim.hu'}
    pdbtms = []
    for path in paths:
        pdbtm_xml = ET.parse(path)
        pdbtm_root = pdbtm_xml.getroot()
        temp = pdbtm_root.attrib.get('ID')
        pdbtm_id = temp.upper()

        if pdbtm_root.attrib.get('TMP') == 'yes':
            chains = []
            for chain_xml in pdbtm_root.findall('pdbtm:CHAIN', ns):
                chain_id = chain_xml.attrib.get('CHAINID')
                num_tm = chain_xml.attrib.get('NUM_TM')
                typ = chain_xml.attrib.get('TYPE')
                seq = chain_xml.find('pdbtm:SEQ', ns)

                regions = []
                for region_xml in chain_xml.findall('pdbtm:REGION', ns):
                    seq_beg = region_xml.attrib.get('seq_beg')
                    pdb_beg = region_xml.attrib.get('pdb_beg')
                    seq_end = region_xml.attrib.get('seq_end')
                    pdb_end = region_xml.attrib.get('pdb_end')
                    typ_region = region_xml.attrib.get('type')
                    region = Region(seq_beg, pdb_beg, seq_end, pdb_end, typ_region)
                    regions.append(region)
                chain = Chain(chain_id, num_tm, typ, seq.text.replace(" ", "").replace("\n", ""), regions)
                chains.append(chain)

            pdbtm = PDBTM(pdbtm_id, chains)
            pdbtms.append(pdbtm)

        else:
            print(pdbtm_id, "is no TMP and was ignored.")

    return pdbtms


paths = []
for file in os.listdir("../pdbtm_xmls"):
    if file.endswith(".xml"):
        path = os.path.join("pdbtm_xmls", file)
        paths.append(path.strip())
pdbtms = parse_pdbtm_xmls(paths)
print(len(pdbtms))

dic = {}
for pdbtm in pdbtms:
    dic[pdbtm.id] = {}
    for chain in pdbtm.chains:
        if chain.typ == "alpha":
            dic[pdbtm.id][chain.id] = []
            for region in chain.regions:
                if region.typ == "H":
                    dic[pdbtm.id][chain.id].append([int(region.pdb_beg), int(region.pdb_end)])

pdb_data = pd.read_csv("../data/test_helices.csv")
print(pdb_data.shape)


def overlap_lst(l1, l2):
    mi = max(l1[0], l2[0])
    ma = min(l1[1], l2[1])
    overls = max(ma - mi + 1, 0)
    return overls

for index, row in pdb_data.iterrows():
    #print(row["PDB ID"], row["Chain"])
    try:
        pdb_data.set_value(index, "Helix Length", row["Helix Length"] + 1)
        SElst = dic[row["PDB ID"]][row["Chain"]]
        SE = [row["Helix Start"], row["Helix End"]]
        overlap = max(list(map(lambda x: overlap_lst(SE, x), SElst)))
        if row["Helix Length"] > 0:
            if overlap / row["Helix Length"] > 0.4 and row["Helix Length"] > 3:  # 70% overlap
                #print(row["PDB ID"],row["Chain"], "overlap!")
                # the +1 should be removed once the original csv is fixed
                pdb_data.set_value(index, "Is Transmembrane", 1)  # depreceated
            elif overlap / row["Helix Length"] == 0 and row["Helix Length"] > 3:
                pdb_data.set_value(index, "Is Transmembrane", -1.0)
    except:
        pass

pdb_data.sort_values(["PDB ID", "Chain", "Is Transmembrane"], ascending=[True, True, True], inplace=True)
pdb_data.to_csv("../data/training_data_new.csv", sep=',', encoding='utf-8', index=False)

