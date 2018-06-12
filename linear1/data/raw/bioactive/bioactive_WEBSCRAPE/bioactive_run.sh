python web_scrape_bioactive.py
sed -e 's/<[^>]*>//g' fasta.txt |sed '/^\s*$/d' > bioactive_peptide.csv
