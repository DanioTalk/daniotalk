if ! [ -x "$(command -v curl)" ]; then
    echo 'Error: curl is not installed.' >&2
    exit 1
fi

if ! [ -x "$(command -v wget)" ]; then
    echo 'Error: wget is not installed.' >&2
    exit 1
fi

if ! [ -x "$(command -v gzip)" ]; then
    echo 'Error: unzip is not installed.' >&2
    exit 1
fi

if ! [ -x "$(command -v tr)" ]; then
    echo 'Error: tr is not installed.' >&2
    exit 1
fi


echo "Downloading Zebrafish genes and associtation files..."
curl https://zfin.org/downloads/all_rna_accessions.txt > zfin_genes.txt
curl https://zfin.org/downloads/gene.txt > ncbi_assoc.txt
curl https://zfin.org/downloads/uniprot.txt > uni_assoc.txt
curl https://zfin.org/downloads/ensembl_1_to_1.txt > ensembl_assoc.txt
rm aliases.txt
curl https://zfin.org/downloads/aliases.txt | grep -i "zdb-gene" | cut -f3,4 >> aliases.txt
curl https://download.alliancegenome.org/5.2.1/GENE-DESCRIPTION-TXT/ZFIN/GENE-DESCRIPTION-TXT_ZFIN_0.txt.gz | gunzip -c | grep -i "zfin:zdb-gene" | cut -f2 >> aliases.txt

cat aliases.txt | sort | uniq > temp123456789
rm aliases.txt && mv temp123456789 aliases.txt
echo "Done.\n"

echo "Downloading orthology data..."
wget https://zfin.org/downloads/human_orthos.txt
echo "Done."

echo "Downloading receptor-ligand pairs data from IID..."
wget http://iid.ophid.utoronto.ca/static/download/human_annotated_PPIs.txt.gz
gzip -d human_annotated_PPIs.txt.gz
cut -f1,2,3,4,8 < human_annotated_PPIs.txt > human_annotated_PPIs_cut.txt
rm human_annotated_PPIs.txt
echo "Done."

echo "Downloading receptor-ligands pairs data from STRING"
wget -O string-protein-links.txt.gz https://stringdb-static.org/download/protein.physical.links.v11.5/7955.protein.physical.links.v11.5.txt.gz
gzip -d string-protein-links.txt.gz
tr " " "\t" < string-protein-links.txt > string-protein-links.tsv
rm string-protein-links.txt

wget -O string-protein-info.txt.gz https://stringdb-static.org/download/protein.info.v11.5/7955.protein.info.v11.5.txt.gz
gzip -d string-protein-info.txt.gz
tr " " "\t" < string-protein-info.txt | sed '1d' > string-protein-info.tsv
rm string-protein-info.txt
echo "Done"

echo "Downloading data from DrugCentral"
wget "https://unmtid-shinyapps.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz" 
gzip -d drug.target.interaction.tsv.gz

echo "Done"
echo "To build database run create_pairs.py script!"