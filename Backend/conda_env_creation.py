#!/usr/bin/env python3
import os

# source conda thing
# conda init
# source ~/.bashrc

# for each environment, create an empty environment and then install packages 
# exported from Digital Ocean conda environments
# currently hardcoded with paths to yml files
os.system("conda create -n genome_assembly -y")
os.system("conda env update -n genome_assembly -f /projects/groupb/Team2-WebServer/Backend/conda/genome_assembly.yml")

os.system("conda create -n quast -y")
os.system("conda env update -n quast -f /projects/groupb/Team2-WebServer/Backend/conda/quast.yml")

os.system("conda create -n gene_prediction -y")
os.system("conda env update -n gene_prediction -f /projects/groupb/Team2-WebServer/Backend/conda/gene_prediction.yml")

os.system("conda create -n microbeannotator -y")
os.system("conda env update -n microbeannotator -f /projects/groupb/Team2-WebServer/Backend/conda/microbeannotator.yml")

os.system("conda create -n comparative_genomics -y")
os.system("conda env update -n comparative_genomics -f /projects/groupb/Team2-WebServer/Backend/conda/comparative_genomics.yml")
