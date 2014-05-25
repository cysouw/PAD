# author  : Jelena Prokić & Johann-Mattis List
# created : 2013-11-27

 
"""
parses PAD data using PAD.prf orthography profile
and aligns it with column annotation on top
"""

from lingpy import *

wl = Wordlist("data.qlc")
wl.tokenize('PAD.prf', column="Simplified")

wl.output('qlc', filename="data_parsed", ignore=['taxa', 'json'])

msa = Alignments('data_parsed.qlc', ref="conceptid")
msa.align(method='library', model=rc('asjp'), gop=-5, iteration=True)
msa.get_consensus(gaps=True)
msa.output('msa')



