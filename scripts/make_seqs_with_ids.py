# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2014-01-15 15:27
# modified : 2014-01-15 15:27
"""
get cognate ids from processed files
"""

__author__="Johann-Mattis List"
__date__="2014-01-15"

from lingpyd import *
from lingpyd.read.qlc import read_msa
from glob import glob

files = glob('../alignments/corrected/*.msa')

wl = Wordlist('data_parsed.qlc')

cids = dict(
        zip(files,range(1001,1000+len(files)+1))
        )

id2cid = dict([(k,0) for k in wl])
cid2con = {}

wrongs = []
for f in files:
    print(f)
    cid = cids[f]
    msa = MSA(read_msa(f, header=True, ids=True))
    con = msa.seq_id
    cid2con[cid] = con
    for i,num in enumerate(msa.ID):
        id2cid[int(num)] = cid

        taxA = msa.taxa[i]
        taxB = wl[int(num),'taxa']
        if taxA.lower() != taxB.lower():
            wrongs += [(str(num),taxA,taxB)]

wl.add_entries("cogid", id2cid, lambda x:x)
wl.output('csv', filename='words', ignore=['taxa', 'json'])

missing = []
for k in wl:
    if wl[k,'cogid'] == 0:
        missing += [(
            str(k),
            wl[k,'concept'],
            wl[k,'doculect'],
            wl[k,'counterpart']
            )]
with open('missing_cogs.txt','w') as f:
    for line in missing:
        f.write('\t'.join(line)+'\n')

with open('wrong_ids.txt','w') as f:
    for line in wrongs:
        f.write('\t'.join(line)+'\n')

