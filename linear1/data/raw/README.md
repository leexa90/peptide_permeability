Mammal_signal_peptide.csv 

Data collection : (actual what I done)
There are two manually curated databases, one which contains a positive dataset of 1000 permeable peptides. (link) We chose mainly signal peptides as the negative dataset as they exit instead of enter cells and also protein data-bank sequences which have experimentally determined 3D structures.
The negative experimentally validated dataset was obtained from three public databases.

1. 5679 protein data-bank sequences of length 4 to 30 
2. 1474 signal peptides from mostly non mammalian sources (link)
3. 1761 signal peptides from mammals (link ,chose those which were 'confirmed')
I removed non-natural amino acids, rare structural bonds and low complexity sequences since the former two lack data and the latter is biologically uninteresting. There were ~500 positive and ~9200 negative peptide sequences left after the removal. The large drop in positive data was because I removed polycationic peptides (> 33.4% R or K), since these peptides are toxic to cells. 
Since there was a ratio of 1 : 18 positive to negative data-points, I duplicated the positive dataset 18 times for the training set to balance the train set. This improved the AUC metric notwithstanding that AUC is already a balanced metric.

sites (web scrape scripts are included with the processing)
http://crdd.osdd.net/raghava/cppsite/

http://www.peptides.be/index.php?p=search&accession_number=&name=&organism_group=&organism_species=&length_from=4&length_to=30&mass_from=&mass_to=&family_group=&family=&uniprot=&aminoacid=&submitbutton=Submit

http://www.signalpeptide.de/index.php
