# DanioTalk <img src="./misc/DaniotalkR logo.jpg" width="139" align="right"/>

# How to cite

Chodkowski M. et al., A ligand-receptor interactome atlas of the zebrafish. iScience 2023, doi: 10.1016/j.isci.2023.107309
URL: https://pubmed.ncbi.nlm.nih.gov/37539027/

# Download

To download it you can go to the Github website and click on the `Download` button.

However if you have a git tool install then you can clone this using a simple command.

```bash
git clone https://github.com/DanioTalk/daniotalk
```

# Installation

There are two ways to install DanioTalk:

- using docker (Linux, Windows, MacOs)
- without docker (Linux only)

## Installation with docker

If you want to use docker go to the docker website and install it 

> `https://docs.docker.com/desktop/install/`

For the purpose of this tutorial we will use the CLI version of Docker.

If you have a poor internet connection or you want to build the up-to-date docker image with up-to-date assets files, then you can go into the downloaded/cloned repository:

```bash
cd <path/to/repository/on/your/machine>
docker build -t mlchodkowski/daniotalk:latest .
```

Now the docker image creation will begin. 

If you have good internet connection then you can pull the image from dockerhub:

```bash
# image location: https://hub.docker.com/repository/docker/mlchodkowski/daniotalk
docker pull mlchodkowski/daniotalk:latest
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
Requirements libraries are provided in the `requirements.txt` file.

```bash
python -m pip install -r requirements.txt
```

# Usage

## Usage with docker

If you want to use docker do the following:

```bash
# If you already have a running container and want to remove it just exec `docker rm -f daniotalk`
docker run -itd --name daniotalk --entrypoint=/bin/bash mlchodkowski/daniotalk:latest
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

##### Example output of 'create_pairs.py'
```bash
docker@ae0aac9f01b5:/daniotalk$ python3 create_pairs.py

______            _     _____     _ _
|  _  \          (_)   |_   _|   | | |
| | | |__ _ _ __  _  ___ | | __ _| | | __
| | | / _` | '_ \| |/ _ \| |/ _` | | |/ /
| |/ | (_| | | | | | (_) | | (_| | |   <
|___/ \__,_|_| |_|_|\___/\_/\__,_|_|_|\_
04:27:58 [INFO 🐠]: Fetching ZF and HUMAN orthology data...
04:28:01 [INFO 🐠]: Loading receptors and ligands...
04:28:01 [INFO 🐠]: Joining association files...
04:28:02 [INFO 🐠]: Reducing redundant records; aggregating frames...
04:28:03 [INFO 🐠]: Joining human orthology...
04:28:03 [INFO 🐠]: Loading Matrisome annotation
04:28:05 [INFO 🐠]: Adding Matrisome annotation
04:28:05 [INFO 🐠]: Adding conservation scores
04:28:08 [INFO 🐠]: Loading data from IID...
04:28:10 [INFO 🐠]: Creating pairs from IID data...
04:28:17 [INFO 🐠]: Creating every possible ligand-receptor pair...
04:28:19 [INFO 🐠]: Creating list of pairs from STRING...
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

DanioTalk repository has a directory `Scripts` containing 4 scripts written in `R`. You can open this scripts and edit their content in comment-highlighted sections (you can edit input filenames and other parameters). These scripts need a dependency files which you should put next to the script itself (in the same directory).

1. **Ligand-Receptor finder for DE genes** (`Scripts/script_v14_DanioTalk LR finder for DE genes.R`)
2. **Ligand-Receptor finder for all expressed genes** (`Scripts/script_v15_DanioTalk LR finder for all expressed genes.R`)
3. **Group visualizer** (`Scripts/script-circ_v3_2 group visualizer.R`)
4. **Multiple group visualizer** (`Scripts/script-circ_v6_LR multiple group visualizer.R`)

## Scripts dependencies

1. For script_v14: Excel file with your singleCell or bulk RNA-seq data with the following columns (`Cell type`, `Gene`, `FC`, `P-value`) - default name: `Data Sheet.xlsx`. Example:
     
     ```
     Cell type        Gene        FC        P-value
     Celltype1        ndr2        7.55    0.00E+00
     Celltype1        spaw        7.16    0.00E+00
     Celltype1        ndr1        7.16    0.00E+00
     Celltype2        ackr3a        7.16    0.00E+00
     Celltype2        ackr3b        7.16    0.00E+00
     Celltype2        ackr4a        7.16    0.00E+00
     ```
   - Generated database file (`Database.csv`) from `Database/` directory
   - `aliases.txt`, `human_orthos.txt` files from the `Data/` directory
   - `Plasma ligands_expt.xlsx` and `Plasma ligands_predicted.xlsx` from `Assets/` directory
   - ***Note!*** If you're using docker you can get these files using `docker cp` command just like when you were copying the `Database.csv` file.

---

2. For script_v15: Excel file with your singleCell or bulk RNA-seq data with the following columns (`Cell type`, `Gene`, `Expression`) - default name: `Data Sheet.xlsx`. Example:
       
       ```
       Cell type        Gene        Expression        
       Celltype1        ndr2        3.00    
       Celltype1        spaw        0.09    
       Celltype1        ndr1        0.04    
       Celltype2        ackr3a        0.45    
       Celltype2        ackr3b        0.31    
       Celltype2        ackr4a        0.30    
       ```
     - Generated database file (`Database.csv`) from `Database/` directory
     - `aliases.txt`, `human_orthos.txt` files from the `Data/` directory
     - `Plasma ligands_expt.xlsx` and `Plasma ligands_predicted.xlsx` from `Assets/` directory
     - 
     - ***Note!*** If you're using docker you can get these files using `docker cp` command just like when you were copying the `Database.csv` file.

---

3. - `Resulting-ligand-receptor-pairs.xlsx` generated by first or second script (**Ligand-Receptor finder**)

---

4. - `Resulting-ligand-receptor-pairs.xlsx` generated by first or second script (**Ligand-Receptor finder**)

## Running

You can just copy all the scripts into new directory and add dependency files. From this directory you can either `cd` into it and run the scripts using command

```bash
Rscript <script-name>.R
```

Or you can open `RStudio`, then open script, use `setwd` to set working directory to the path containing scripts & assets files and click on `Source` to run whole script or click `Run` to run single line at the time.
