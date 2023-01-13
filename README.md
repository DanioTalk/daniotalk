# DanioTalk <img src="./misc/DaniotalkR logo.jpg" width="139" align="right"/>

# How to cite 
Chodkowski M. et al., A ligand-receptor interactome atlas of the zebrafish. bioRxiv 2022, doi: https://doi.org/10.1101/2022.12.15.520415

# Download
To download it you can go to the Github website and click on the `Download` button.

However if you have a git tool install then you can clone this using a simple command.
```bash
git clone https://github.com/777moneymaker/daniotalk
```

# Installation
There are two ways to install DanioTalk:
- using docker (Linux, Windows, MacOs)
- without docker (Linux only)

## Installation with docker
If you want to use docker go to the docker website and install it 
>`https://docs.docker.com/desktop/install/`

For the purpose of this tutorial we will use the CLI version of the Docker (we will not use the graphical tool, only the command-line interface).

If you have a poor internet connection or you want to build the up-to-date docker image with up-to-date assets files, then you can go into the downloaded/cloned repository:
```bash
cd <path/to/repository/on/your/machine>
docker build -t mlchodkowski/daniotalk:1.0 .
```
Now the docker image creation will begin. 

If you have good internet connection then you can pull the image from dockerhub:
```bash
# image location: https://hub.docker.com/repository/docker/mlchodkowski/daniotalk
docker pull mlchodkowski/daniotalk:1.0
```

## Installation without docker
Enter the directory `daniotalk`. From this directory execute the following commands:
```bash
cd Data/
./fetch_data.sh
cd ..
```
Now you have all the files downloaded and you can run
```bash
python create_pairs.py
```

### Requirements
- python >= 3.9
- pandas
- sqlalchemy
- openpyxl

```bash
python -m pip install pandas sqlalchemy openpyxl
```

# Usage
## Usage with docker
If you want to use docker do the following:
```bash
# If you already have a running container and want to remove it just exec `docker rm -f daniotalk`
docker run -itd --name daniotalk --entrypoint=/bin/bash mlchodkowski/daniotalk:1.0
docker exec -it daniotalk bash
# You should see something like this at the beginning of the line (it means that container was created successfully)
# docker@158db1d9cfa1:/daniotalk$

python create_pairs.py
```
## Usage without docker
If you have the software installed locally (not using docker) then you have to just `cd` into the cloned repository.
```bash
cd <path/to/repository/on/your/machine>
cd Data/
bash fetch_data.sh # If you already have the files downloaded then you can skip this step
cd ..
python create_pairs.py
```

### Output files
If you used the software locally (not installed using docker) then the `Database.xlsx` file will be saved in the `Database/` directory.

However if you used docker then you have to copy the results from the container.
First exit the container
```bash
# docker@158db1d9cfa1:/daniotalk$
exit
```
Then from the powershell or any other terminal execute this command:
```bash
docker cp daniotalk:/daniotalk/Database/Database.csv <path/to/save/location>
# Example
# docker cp daniotalk:/daniotalk/Database/Database.csv C:/Users/user1/Desktop/
```

To remove the running container:
```bash
docker stop daniotalk
docker rm daniotalk
```

**Congratulations!** Now you have successfully executed the scripts using **Docker**.

# Scripts
DanioTalk repository has a directory `Scripts` containing 3 scripts written in `R`. You can open this scripts and edit their content in comment-highlighted sections (you can edit input filenames and other parameters). These scripts need a dependency files which you should put next to the script itself (in the same directory).

1. **Ligand-Receptor finder for DE genes** (`Scripts/script_v12_DanioTalk LR finder for DE genes.R`)
2. **Ligand-Receptor finder for all expressed genes** (`Scripts/script_v13_DanioTalk LR finder for all expressed genes.R`)
2. **Group visualizer** (`Scripts/script-circ_v2_2 group visualizer.R`)
3. **Multiple group visualizer** (`Scripts/script-circ_v5_LR multiple group visualizer.R`)


## Scripts dependencies
1.  
    - Generated database file (`Database.csv`) from `Database/` directory
    - For script_v12: Excel file with your singleCell data with the following columns (`Cell type`, `Gene`, `FC`, `P-value`) - default name: `Data Sheet.xlsx`. Example:
        ```
        Cell type	    Gene		FC	    P-value
        Celltype1       ndr2		7.55	0.00E+00
        Celltype1	    spaw		7.16	0.00E+00
        Celltype1	    ndr1		7.16	0.00E+00
        Celltype2	    ackr3a		7.16	0.00E+00
        Celltype2	    ackr3b		7.16	0.00E+00
        Celltype2	    ackr4a		7.16	0.00E+00
        ```
     - For script_v13: Excel file with your singleCell data with the following columns (`Cell type`, `Gene`, `Expression`) - default name: `Data Sheet.xlsx`. Example:
        ```
        Cell type	    Gene		Expression	    
        Celltype1	    ndr2		3.00	
        Celltype1	    spaw		0.09	
        Celltype1	    ndr1		0.04	
        Celltype2	    ackr3a		0.45	
        Celltype2	    ackr3b		0.31	
        Celltype2	    ackr4a		0.30	
        ```
    - `aliases.txt`, `human_orthos.txt`, `drug.target.interaction.tsv` files from the `Data/` directory
    - `Plasma ligands_expt.xlsx` and `Plasma ligands_predicted.xlsx` from `Assets/` directory
    - `GO ID_db.xlsx` from `Assets/` directory
    - ***Note!*** If you're using docker you can get these files using `docker cp` command just like when you were copying the `Database.csv` file.
---
2. 
    - `Resulting-ligand-receptor-pairs.xlsx` generated by first script (**Ligand-Receptor finder**)
---
3. 
    - `Resulting-ligand-receptor-pairs.xlsx` generated by first script (**Ligand-Receptor finder**)


## Running
You can just copy all the scripts into new directory and add dependency files. From this directory you can either `cd` into it and run the scripts using command
```bash
Rscript <script-name>.R
```
Or you can open `RStudio`, then open script, use `setwd` to set working directory to the path containing scripts & assets files and click on `Source` to run whole script or click `Run` to run single line at the time.
