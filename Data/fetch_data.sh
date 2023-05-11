# This checks for `curl` presence
if ! [ -x "$(command -v curl)" ]; then
    echo 'Error: curl is not installed.' >&2
    exit 1
fi
# This checks for `wget` presence
if ! [ -x "$(command -v wget)" ]; then
    echo 'Error: wget is not installed.' >&2
    exit 1
fi
# This checks for `gzip` presence
if ! [ -x "$(command -v gzip)" ]; then
    echo 'Error: unzip is not installed.' >&2
    exit 1
fi
# This checks for `tr` presence
if ! [ -x "$(command -v tr)" ]; then
    echo 'Error: tr is not installed.' >&2
    exit 1
fi

# Part below downloads all the association files (metadata for each zebrafish gene)
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

# This downloads file containing conservation scores
wget -O conservation.tsv.gz https://fms.alliancegenome.org/download/ORTHOLOGY-ALLIANCE_COMBINED.tsv.gz
gunzip conservation.tsv.gz
echo "Done.\n"

# This downloads human orthologs for zebrafish genes
echo "Downloading orthology data..."
wget https://zfin.org/downloads/human_orthos.txt
echo "Done."

# This downloads IID pairs and cuts the columns of interest
echo "Downloading receptor-ligand pairs data from IID..."
wget http://iid.ophid.utoronto.ca/static/download/human_annotated_PPIs.txt.gz
gzip -d human_annotated_PPIs.txt.gz
cut -f3,4,8 < human_annotated_PPIs.txt > human_annotated_PPIs_cut.txt
rm human_annotated_PPIs.txt
echo "Done."

# This downloads receptor-ligand pairs from STRING database
echo "Downloading receptor-ligands pairs data from STRING"
wget -O string-protein-links.txt.gz https://stringdb-static.org/download/protein.physical.links.v11.5/7955.protein.physical.links.v11.5.txt.gz
gzip -d string-protein-links.txt.gz
tr " " "\t" < string-protein-links.txt > string-protein-links.tsv
rm string-protein-links.txt
# Second part of STRING download
wget -O string-protein-info.txt.gz https://stringdb-static.org/download/protein.info.v11.5/7955.protein.info.v11.5.txt.gz
gzip -d string-protein-info.txt.gz
tr " " "\t" < string-protein-info.txt | sed '1d' > string-protein-info.tsv
rm string-protein-info.txt
echo "Done"

# This downloads data from DrugCentral
echo "Downloading data from DrugCentral"
wget "https://unmtid-shinyapps.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz" 
gzip -d drug.target.interaction.tsv.gz

# Downloads data from CellTalkDB
echo "Downloading data from CellTalkDB"
wget "https://github.com/ZJUFanLab/CellTalkDB/releases/download/v1.0/scsrctdb-1.0.tar.gz"
tar xvzf "scsrctdb-1.0.tar.gz"
cp scsrctdb-1.0/data/LRdb.rda ./LRdb.rda
rm -rf scsrctdb-1.0 
rm scsrctdb-1.0.tar.gz

echo "Downloading data from CellCellInteractions"
wget "https://zenodo.org/record/7589953/files/receptor_ligand_interactions_mitab_v1.0_April2017.txt.gz"
gunzip -c "receptor_ligand_interactions_mitab_v1.0_April2017.txt.gz" | cut -f1,2 >> cellcellinteractions.txt

echo "Downloading Ramilowski dataset"
curl https://fantom.gsc.riken.jp/5/suppl/Ramilowski_et_al_2015/data/PairsLigRec.txt | cut -f2,4 > ramilowski.txt

echo "Downloading data from CellPhoneDB"
wget "https://github.com/ventolab/cellphonedb-data/archive/refs/tags/v4.1.0.tar.gz"
tar xvzf "v4.1.0.tar.gz"
cat cellphonedb-data-4.1.0/data/interaction_input.csv | sed "s/_HUMAN//g"  > cellphonedb.csv
rm -rf "cellphonedb-data-4.1.0"
rm v4.1.0.tar.gz

echo "Downloading ligand-type annotations for Zebrafish and Human"
wget "http://current.geneontology.org/annotations/zfin.gaf.gz"
wget "http://geneontology.org/gene-associations/goa_human.gaf.gz"
gunzip -c zfin.gaf.gz | cut -f3,5 > zfin.gaf
gunzip -c goa_human.gaf.gz | cut -f3,5 > goa_human.gaf

echo "Downloading CellCellChat annotation"
wget "https://github.com/sqjin/CellChat/archive/refs/tags/v1.5.0.tar.gz"
tar xvzf "v1.5.0.tar.gz"
cp CellChat-1.5.0/data/CellChatDB.human.rda .
rm -rf CellChat-1.5.0
rm v1.5.0.tar.gz

echo "Downloading Matrisome annotation"
wget https://github.com/Matrisome/MatrisomeAnalyzeR/raw/main/data/matrisome.list.rda

echo "Done"
echo "To build database run create_pairs.py script!"