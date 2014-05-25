PAD data
========


PROBLEMATIC FILES: 45, 42, 37, 38, 27, 3 53, 89, 104


DETAILED TO DO LIST:

- first draft of abstract (JP: 24 Dec 2013)
- final abstract (JP, JML, MC: 10 Jan 2014)
- DONE: go through files manually corrected by MC (JP & JML: 22 Dec 2013) 
- DONE: add Wenker sentences (JP & JML: 22 Dec 2013)
- DONE: add translations to concepts.csv (JML: 24 Dec 2013)
- DONE: change village names (JP: 15 Dec 2013)
- DONE: add json file (JP: 22 Dec 2013)
- DONE: change annotation of rows (JP: 22 Dec 2013)
- DONE: change file names (JP: 22 Dec 2013)


ANNOTATION CONVENTIONS:

- 0      CROSSED...	. . + - - + . .
- 0      COMPLEX...	. . < - - > . .
- 0      IGNORE....	. . x . . x . .


TODO:

- add Wenker sentence at the beginning of each file (e.g. # wenker sentence : "Im Winter fliegen die trockenen Blätter in der Luft herum.")
- DONE: change village names
- DONE: rename manually corrected files into concept\_id.msa (e.g. wort\_19.msa)
- POSTPONED: find a better solution for the annotation rows (0	CROSS --> :CROSS	...)


PAD data set comes from ...


In this package, the data can be found in the following files and formats:

- Sources/originaldata_.csv
	- original file downloaded form the Deutscher Sprachatlas website (?) without any changes from our side
	
- Sources/originaldata\_corrected.csv 
	- the same format as the originaldata.csv file
	- ʋ changed to ɒ where symbol 'Q' was used in the sampa transcription (considered typo)
	- lɤ changed to ɤ (considered typo)

- Strings/data.qlc
	- qlc file format of the originaldata\_corrected.csv 
	- multiple pronunciations split into separate rows

- Strings/data_parsed.qlc
	- format: qlc wordlist format
	- columns: ID, CONCEPT,COUNTERPART, CONCEPTID, DOCULECT, TOKENS
	- column TOKENS contains parsed strings where the following tokens are removed:
		- symbol for syllable break (".")  (not consistently annotated)
		- symbol for continuous speech ("‿")  (not consistently annotated)
		- symbol for primary stress 
		- symbol for secondary stress 


- Alignments
	- format: directory with the msa files for each concept  separately 

The data is parsed using PAD.prf file. This file lists all unique phonemes that occur in the data.

Scripts

- pad2qlc.py: transforms original pad data to qlc format
	- input: pad\_original\_data\_corrected.csv
	- output: pad\_data.csv
- pad\_align.py: parses the data and multi-aligns it
	- input: data.qlc
	- output:
		- data\_parsed.qlc
		- PAD\_parsed-msa/* 


