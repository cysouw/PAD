# author      : Jelena Prokić <prokic@uni-marburg.de>
# created     : 2013-11-02
# last update : 2013-12-15

"""
converts PAD original data format to qlc lingpy format
usage: pad2qlc.py ../sources/originaldata_corrected.csv ../doculects.csv
"""

import sys
import re
import unicodedata
from time import gmtime, localtime, strftime


pad = open(sys.argv[1], "r").readlines() # pad_data_original_corrected.csv
pad = pad[1:] #remove header
names = open(sys.argv[2], 'r').readlines() #doculects.csv
qlc = open("data.qlc", "w")

id_dict = {}

for l in names:
    if not (l.startswith("#") or l.startswith("ID")):
        loc_id = l.split("\t")[0]
        name = l.split("\t")[4].rstrip()
        id_dict[loc_id] = name


qlc.write("@source: " + sys.argv[1] +"\n")
qlc.write("@date: "+strftime("%Y-%m-%d %H:%M:%S", localtime())+"\n")
qlc.write("ID"+"\t"+"CONCEPT"+"\t"+"COUNTERPART"+"\t"+"CONCEPTID"+"\t"+"DOCULECT"+"\n")

id = 1
for line in pad:

    line = line.strip()
    line = unicodedata.normalize("NFD", line)
    tokens = line.split(";")
    concept_id = tokens[0].strip()
    concept = tokens[1].strip().strip("\"")
    doculect_id = tokens[2].strip()
    ipa_string = tokens[3].strip().strip("\"")
    sampa_string = tokens[4].strip()


    #cases where there are multiple forms in the same row separated by "/" or " "
    if "/" in ipa_string.strip() or " " in ipa_string.strip():
        forms = re.split( '/|\s',ipa_string)
        forms = [x for x in forms if x != ""]
        for form in forms:
            form = form.strip()
            qlc.write(str(id)+"\t"+concept+"\t"+form+"\t"+concept_id+"\t"+id_dict[doculect_id]+"\n")
            id += 1


    #all other cases, except a missing entry line
    elif ("|" not in ipa_string) and (len(ipa_string) > 0):
        qlc.write(str(id)+"\t"+concept+"\t"+ipa_string+"\t"+concept_id+"\t"+id_dict[doculect_id]+"\n")
        id += 1
