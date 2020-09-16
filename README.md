# KDDANet
A computational framework for uncovering molecular roadmap mediating Known Drug-Disease Associations (KDDAs)

# Version
1.0.0

# Author
Hua Yu and Lu Lu

# Description
This software identifies hidden gene interactions and modules mediating KDDAs through implementing minimum cost optimization, combined with depth-first searching and graph clustering on a unified flow network model

# Installation
KDDANet can be installed as follows:<br>
* I) Add KDDANet to PATH variable
```Bash
git clone https://github.com/huayu1111/KDDANet.git
export PATH=path_to_KDDANet:$PATH
cd path_to_KDDANet
chmod +x lpsolver
```
* II) Download MCL software (Version 'mcl-14-137') from https://micans.org/mcl/ and install:
```Bash
tar xzf mcl-14-137.tar.gz
cd mcl-14-137
./configure --prefix=$HOME/local
make install
export PATH=$HOME/local/bin:$PATH
```
# Running
* I) Step 1: Obtaining the solution of minimum cost flow optimization problem for each query drug (disease) and its related diseases (drugs)
```Bash
perl getSolution.pl -dr <DRUG_FILE> -di <DISEASE_FILE> -n <NET_FILE> -dt <DRUG_TARGET_FILE> -dg <DISEASE_GENE_FILE> --dd <DRUG_DISEASE_ASSOCIATION_FILE> -o <OUTPUT_DIR> -gmin <GAMMA_MIN> -gmax <GAMMA_MAX> -gstep GAMA_INCREMENT -opt <CONTEXT> -h
-dr  <DRUG_FILE>, A single column text file (default DrugBank drugs)
-di  <DISEASE_FILE>, A single column text file (default OMIM diseases)
-n  <NET_FILE>, A tab-delimited text file (default HumanNet)
-dt <DRUG_TARGET_FILE> A tab-delimited text file (default DrugBank gene targets)
-dg <DISEASE_GENE_FILE> A tab-delimited text file (default disease-related genes extracted from literature [1])
-dd <DRUG_DISEASE_ASSOCIATION_FILE> A tab-delimited text file (default drug-disease association contained in CTD database)
-o  OUTPUT_DIR for output files (default './result/')
-gmin  GAMMA_MIN	(default 4)
-gmax  GAMMA_MAX	(default 12)
-gstep  GAMMA_INCREMENT (default 1)
-opt <CONTEXT> (SDrTDi or SDiTDr, default SDrDTi)
-h  help
```
Input files format:
1. <DRUG_FILE> / <DISEASE_FILE><br>
A single column text file provides drug or disease information<br>
Format: DRUG_ID/DISEASE_ID
2. <NET_FILE><br>
A tab-delimited text file provides interactome network information<br>
Format: G1_ID G2_ID	Weight

3. <DRUG_TARGET_INFO> / <DISEASE_GENE_INFO><br>
A tab delimited text file provides drug's target information / disease-related gene information<br>
Format: DRUG_ID  TARGET_GENE_ID / DISEASE_ID RELATED_GENE_ID 

4. <DRUG_DISEASE_ASSOCIATION_FILE><br>
A tab delimited text file provides drug-disease association information<br>
Format: DRUG_ID  DISEASE_ID

An example for running this step:
```Bash
perl getSolution.pl -dr ../inputdata/DrugBank.Drug.Info.txt -di ../inputdata/OMIM.Disease.Info.txt -n ../inputdata/HumanNet.txt -dt ../inputdata/Used_Drug_Target_Data.txt -dg ../inputdata/Used_Disease_Gene_Data.txt --dd ../inputdata/KDDAs_Total.txt -o ../result/ -gmin 4 -gmax 4 -gstep 1 -opt SDrTDi
```
This command will obtain the solution for each query drug/disease and its related diseases/drugs in the directory of ../result/$gamma for each gamma value
* II) Step2: Extracting roadmap of genes mediating individual KDDA from the solution using depth-first searching<br>
```Bash
perl getNodeFlowAndEdgeFlow.pl --indir <IN_DIR> --outdir <OUT_DIR>
--indir <IN_DIR> the directory used to place result files of solutions
--outdir <OUT_DIR> the directory used to place node flow file and edge flow file
```
Or:
```Bash
 perl getGeneInteractionSubNetForEachKDDAUsingDFS.pl --infile <IN_FILE> --outDir <OUT_DIR>
--infile <IN_FILE> the edge flow file
--outDir <OUT_DIR> the directory used to place gene interaction subnetwork for each known drug-disease association
```
An example for running this step:
```Bash
perl getNodeFlowAndEdgeFlowFromSolutions.pl --indir ../results/4 --outdir  ../results/4
perl getGeneInteractionSubNetForEachKDDAUsingDFS.pl.pl --infile ../results/4/opti_edge_flow.txt --outDir  ../results/4/subnetworks
```
* III) Identifying gene modules from subnetwork using MCL algorithm
```Bash
Command:ObtainModuleInfoByMCL.pl --indir <IN_DIR> --outdir <OUT_DIR>
--indir <IN_DIR> the directory used to place result files of gene interaction subnetworks
--outdir <OUT_DIR> the directory used to place gene module files
```
An example for running this step:
```Bash
perl ObtainModuleInfoByMCL.pl --indir ../results/4/subnetworks --outDir  ../results/4/modules
```
# Input data
The directory of [./inputdata/](https://github.com/huayu1111/KDDANet/tree/master/inputdata) provides the datasets used in our paper

# License
KDDANet is licensed under the GPL version 3 or any later version

# For any questions, please contact:
Hua Yu (yuhua200886@163.com) or Lu Lu (tkrwy@126.com)
