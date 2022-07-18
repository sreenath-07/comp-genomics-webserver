#!/usr/bin/env python3

"""
The above script is the one which talks to django
This script will take the path of the input files , user name and the output files location pass it to the Genome assembly
Access the output files and throw it to django
It should throw the results and also errors if any
Generate the log file and maintain a record
Check the input file format also check if the files are zipped
 """

"""
Generate a temp folder
Final output to one folder
Conda Environment
Requirements.txt
"""

"""
Input: Django call
1.files are in zipped fastq format (.fq.gz) or not
2.File format
2.Username and email id
3.Files location
"""

"""
Ouput: Django call
1.files are in zipped fastq format (.fq.gz)
2.Username and email id
3.Files location
4.Send it to django back
5.Throw any errors
"""

import sys
import uuid
import argparse
import subprocess
import os
from sys import exit
from Bio import SeqIO
from subprocess import PIPE, run
import argparse
from os import walk
import shutil

def create_directories(output_dir):
    ## Removing Output Directory if it already exists
    if os.path.exists(output_dir) and os.path.isdir(output_dir):
        shutil.rmtree(output_dir)
    ## Creating a new output directory
    os.mkdir(output_dir)

def create_folder_structure(user_output_folder_path):
    ## ToDo: Need to delete folder if it already exists during creation
    try:
        # os.mkdir(user_output_folder_path) #Already created and given to backend from views.py
        create_directories(user_output_folder_path+'/genome_assembly')
        quast_output_path = user_output_folder_path + '/genome_assembly/quast/'
        create_directories(quast_output_path)
        create_directories(user_output_folder_path+'/assembled_contigs')

        create_directories(user_output_folder_path+'/gene_prediction')

        create_directories(user_output_folder_path+'/functional_annotation')

        create_directories(user_output_folder_path+'/comparative_genomics')
        create_directories(user_output_folder_path + '/comparative_genomics/fastani/')
        create_directories(user_output_folder_path + '/comparative_genomics/fastmlst/')
        create_directories(user_output_folder_path + '/comparative_genomics/amr_finder/')

        create_directories(user_output_folder_path+'/final_results')
    except:
        print("Error with creating sub-directories")

def run_genome_assembly(input_path,user_output_folder_path):

    print("Running Genome Assembly Pipeline")
    contigs_path = user_output_folder_path + '/assembled_contigs'
    genome_assembly_output = user_output_folder_path + '/genome_assembly'

    #Running Shovill
    run_shovill(input_path, user_output_folder_path, contigs_path)

    #Running Quast
    quast_output_path = user_output_folder_path + '/genome_assembly' + '/quast/'
    run_quast(contigs_path, quast_output_path)

    tar_command = f'tar -C {genome_assembly_output} -czvf {user_output_folder_path}/genome_assembly.tar.gz .'
    print(tar_command)
    os.system(tar_command)

    copy_command = f'cp {user_output_folder_path}/genome_assembly.tar.gz {user_output_folder_path}/final_results/genome_assembly.tar.gz'
    print(copy_command)
    os.system(copy_command)

def run_shovill(input_path,user_output_folder_path, contigs_path):
    """
    Assumption: Separate Foldes for each isolate
    """
    try:
        output_path = user_output_folder_path + '/genome_assembly'
        conda_env = "conda run -n webserver"
        for isolate_id in os.listdir(input_path):
            print(isolate_id)
            shovill_command = f"{conda_env} shovill --gsize 5M --cpus 4 --trim --outdir {output_path}/{isolate_id} --R1 {input_path}/{isolate_id}/*_1.fq.gz --R2 {input_path}/{isolate_id}/*_2.fq.gz --force"
            copy_command = f"cp {output_path}/{isolate_id}/contigs.fa {contigs_path}/{isolate_id}_contigs.fa"
            print(shovill_command)
            print(copy_command)

            os.system(shovill_command)
            os.system(copy_command)

    except:
        print("Error running Shovill")
        # sys.exit()

def run_quast(input_contigs_dir, output_dir):
    """
    :param output_dir:
    :return:
    """
    try:
        conda_env = "conda run -n quast"
        quast_command = f"{conda_env} quast {input_contigs_dir}/* -o {output_dir}"
        print(f"Running Quast for all assemblies with command: {quast_command} \n \n")
        #subprocess.call(shlex.split(quast_command))
        os.system(quast_command)
    except:
        print("Error: Quast command didnt work")
        # sys.exit()

def run_prokka(user_output_folder_path):
    try:
        print("Running Prokka")
        prokka_output = user_output_folder_path + '/gene_prediction'
        contigs_path = user_output_folder_path + '/assembled_contigs'
        functional_annotation_path  = user_output_folder_path + '/functional_annotation'

        for x in os.listdir(contigs_path):
            print(x)
            prefix=(x.split(".")[0])
            command = f'conda run -n gene_prediction prokka --outdir {prokka_output}/{x} {contigs_path}/{x} --kingdom Bacteria --rfam --prefix {prefix} --force --cpus 4 --compliant --centre --usegenu --genus Escherichia --centre UNCC --quiet'
            print(f"Executing: {command}")
            os.system(command)
            # subprocess.call(shlex.split(command))

        # Copy all the .faa files to functional annotation folder
        copy_faafiles_fa = f'cp -r {prokka_output}/*/*.faa {functional_annotation_path}'
        print(copy_faafiles_fa)
        os.system(copy_faafiles_fa)
        
        tar_command = f'tar -C {prokka_output} -czvf {user_output_folder_path}/gene_prediction_prokka.tar.gz .'
        print(tar_command)
        os.system(tar_command)
        
        copy_command = f'cp {user_output_folder_path}/gene_prediction_prokka.tar.gz {user_output_folder_path}/final_results/gene_prediction_prokka.tar.gz'
        print(copy_command)
        os.system(copy_command)
        
    except:
        print("Error: Prokka failed")
        # with open("/home/Team5/Team2-WebServer/sample_inputs/errors.log","a+") as log:
        #     log.write("Prokka failed \n")

def run_microbeannotator(user_output_folder_path):
    try:
        fa_path = user_output_folder_path+"/functional_annotation"
        os.chdir(fa_path)
        d_path = "/projects/groupb/MicrobeAnnotator_DB"
        command = f'conda run -n microbeannotator microbeannotator -i  $(ls *.faa)  -d {d_path} -o {fa_path}  -m blast -p 1 -t 4  --light --refine'
        print(command)
        os.system(command)

        tar_command = f'tar -C {fa_path} -czvf {user_output_folder_path}/functional_annotation.tar.gz .'
        print(tar_command)
        os.system(tar_command)

        copy_command = f'cp {user_output_folder_path}/functional_annotation.tar.gz {user_output_folder_path}/final_results/functional_annotation.tar.gz'
        print(copy_command)
        os.system(copy_command)

    except:
        print("Error: Microbeannotator failed")
        # with open("/home/Team5/Team2-WebServer/sample_inputs/errors.log","a+") as log:
        #     log.write("microbeannotator failed \n")

def run_comparative_genomics(user_output_folder_path):

    comparative_genomics_output = user_output_folder_path + '/comparative_genomics'

    try:
        run_fastani(user_output_folder_path)
    except:
        print("Error: Fast_ani failed")

    try:
        # os.system("source /projects/groupb/miniconda3/etc/profile.d/conda.sh")
        # os.system("conda init")
        # os.system("conda activate comparative_genomics")
        run_fast_ani_vis()
        # os.system("conda deactivate")
    except:
        print("Error: Fast_ani_vis failed")

    try:
        run_amr_finder(user_output_folder_path)
    except:
        print("Error: AMR_Finder failed")

    try:
        run_fastmlst(user_output_folder_path, "temp_output")
    except:
        print("Error: fastmlst failed")


    tar_command = f'tar -C {comparative_genomics_output} -czvf {user_output_folder_path}/comparative_genomics.tar.gz .'
    # print(tar_command)
    os.system(tar_command)

    copy_command = f'cp {user_output_folder_path}/comparative_genomics.tar.gz {user_output_folder_path}/final_results/comparative_genomics.tar.gz'
    # print(copy_command)
    os.system(copy_command)

def run_fastani(user_output_folder_path):
    """
    input_dir: Input directory of genomes
    """

    contigs_path = user_output_folder_path + '/assembled_contigs'
    global fastani_output
    fastani_output = f'{user_output_folder_path}/comparative_genomics/fastani'
    dirs = os.listdir(contigs_path)

    contigs_path_list=f'find {contigs_path} -type f -name "*.fa" > {fastani_output}/files.txt'
    os.system(contigs_path_list)
    conda_env = "conda run -n comparative_genomics"
    fastani_command = f"{conda_env} fastANI --ql {fastani_output}/files.txt --rl {fastani_output}/files.txt --matrix -o {fastani_output}/output.out"
    os.system(fastani_command)

    global input_matrix_path
    input_matrix_path = f'{fastani_output}/output.out.matrix'
    

def run_fast_ani_vis():
    """
    Author:  Kenji Nishiura
    Description: Generates ANI tree from fastANI output matrix.
    input_matrix_path = input fastani matrix output
    output_prefix = prefix for tree and heatmap
    """
    import pandas as pd
    import csv
    from pathlib import Path
    import scipy.cluster.hierarchy as hc
    from scipy.cluster.hierarchy import ClusterNode
    from Bio import Phylo
    import matplotlib.pyplot as plt

    def convert_scipy_to_newick(
            node=ClusterNode, parent=float, leaf=list([str]), newick=str("")
    ):
        if node.is_leaf():
            return f"{leaf[node.id]}:{(parent - node.dist):.2f}{newick}"
        else:
            if len(newick) > 0:
                newick = f"):{(parent - node.dist):.2f}{newick}"
            else:
                newick = ");"
            newick = convert_scipy_to_newick(node.left, node.dist, leaf, newick)
            newick = convert_scipy_to_newick(node.right, node.dist, leaf, f",{newick}")
            newick = f"({newick}"
            return newick

    SampleID = []
    ani_row = []
    with open(f"{input_matrix_path}", "r") as f:
        reader = csv.reader(f, delimiter="\t")
        genome_num = int(next(reader)[0].rstrip("\n"))

        for row in reader:
            SampleID.append(Path(row[0]).with_suffix("").name.rstrip(".fa"))
            ani_values = list(map(lambda d: float(d), row[1:]))
            ani_values.extend([0] * (genome_num - len(ani_values)))
            ani_row.append(ani_values)

    df = pd.DataFrame(data=ani_row, columns=SampleID, index=SampleID)

    # fill in diagonal and mirror values across diagonal
    for selfmatch in range(genome_num):
        df.iat[selfmatch, selfmatch] = 100
    for i, id in enumerate(SampleID):
        for j, ani in enumerate(df[id][i:]):
            df.iat[i, i + j] = ani

    # Hierarchical clustering
    linkage = hc.linkage(df, method="average", metric="seuclidean", optimal_ordering=True)

    # Output newick tree and save to image
    tree = hc.to_tree(linkage)
    if isinstance(tree, ClusterNode):
        output_file = f"{fastani_output}/fastani.nwk"
        with open(output_file, "w") as f:
            newick_tree = convert_scipy_to_newick(tree, tree.dist, list(df.columns))
            f.write(newick_tree)
        with open(output_file, "r") as f:
            image_tree = Phylo.read(output_file, "newick")
            image_tree.root_at_midpoint()
            fig = plt.figure(figsize=(10, 20), dpi=300)
            axes = fig.add_subplot(1, 1, 1)
            Phylo.draw(image_tree, axes=axes)
            plt.show()
            plt.savefig(f"{fastani_output}/fastani.png")

def run_amr_finder(user_output_folder_path):
    """
    input_contigs_dir: Path to directory contianing contigs
    """
    organism = "Escherichia"
    input_contigs_dir = user_output_folder_path + '/assembled_contigs'
    amr_finder_output = f'{user_output_folder_path}/comparative_genomics/amr_finder'

    create_amr_db = f'conda run -n comparative_genomics amrfinder --force_update'
    os.system(create_amr_db)

    conda_env = "conda run -n comparative_genomics"
    run_amr_finder = f'for assembly in {input_contigs_dir}/*; do {conda_env} amrfinder -n "$assembly" --organism {organism} --output {amr_finder_output}/"$(basename $assembly)".out; done'
    os.system(run_amr_finder)

    try:
        ##ToDO: Kenji please fix this
        # os.system("source /projects/groupb/miniconda3/etc/profile.d/conda.sh")
        # os.system("conda activate comparative_genomics")

        print("Trying to create amr_finder image")
        import dataframe_image as dfi
        import pandas as pd
        import matplotlib

        # load output tsv, store in dataframe, save as image
        for amrfile in os.listdir(amr_finder_output):
           df = pd.read_csv(f"{amr_finder_output}/{amrfile}", sep="\t", header=0, index_col=0)
           dfi.export(df, f"{amr_finder_output}/{amrfile}.png", table_conversion='matplotlib')

        # os.system("conda deactivate")
    except:
        print("Error: AMR Finder image creation failed")

def run_fastmlst(user_output_folder_path, output_file_name):
    """
    Author:  Kenji Nishiura
    Description:  identifies E. Coli MLST profile, and then infers approximately-maximum-likelihood phylogenetic tree using GTR substitution model
    """
    import pandas as pd
    import dataframe_image as dfi

    cpu_threads = 4
    input_contigs_dir = user_output_folder_path + '/assembled_contigs'
    fastmlst_output = f'{user_output_folder_path}/comparative_genomics/fastmlst'
    output_fasta = f'{fastmlst_output}/{output_file_name}_fasta'
    output_tab = f'{fastmlst_output}/{output_file_name}_tab'
    output_tree = f'{fastmlst_output}/{output_file_name}_tree'

    conda_env = 'conda run -n comparative_genomics'

    # update  MLST database
    updateDB = f"{conda_env} fastmlst --update-mlst"
    os.system(updateDB)

    # run fastmlst
    mlst = f"{conda_env} fastmlst -t {cpu_threads} -v 2 -sch ecoli#1 -fo {output_fasta} -to {output_tab} -n {fastmlst_output}/novel.fasta {input_contigs_dir}/*"
    os.system(mlst)

    # load fastmlst tabular output, store in dataframe, save as image
    df = pd.read_csv(f"{output_tab}", header=0, index_col=0)
    dfi.export(df, f"{output_tab}.png", table_conversion='matplotlib')

    # align concatenated gene fasta file
    align = f"{conda_env} mafft --auto --thread {cpu_threads} {output_fasta} > {fastmlst_output}/tmp.aln"
    os.system(align)

    # trim alignment
    trim = f"{conda_env} trimal -in {fastmlst_output}/tmp.aln -out {fastmlst_output}/tmp_trimmed.aln -automated1"
    os.system(trim)

    try:
        ##ToDO: Kenji please fix this
        # build tree
        os.system("conda activate comparative_genomics")
        tree = f'FastTree -nt -gtr < {fastmlst_output}/tmp_trimmed.aln > {output_tree}'
        os.system(tree)
        os.system("conda deactivate")
    except:
        print("Error: fastmlst build tree failed (mostly due to conda init issue")

def tar_output(user_output_folder_path):
    tar_command = f'tar -C {user_output_folder_path}/final_results/ -czvf {user_output_folder_path}/final_results/final_results.tar.gz .'
    os.system(tar_command)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-input_path', help='Enter the Input path of the file')
    args = parser.parse_args()
    user_output_folder_path = args.input_path
    user_output_folder_path = user_output_folder_path.rstrip("/")

    input_path = f'{user_output_folder_path}/input'

    # print(input_path,user_output_folder_path)

    ####### Pipeline #######
    create_folder_structure(user_output_folder_path) #Creating the folders
    run_genome_assembly(input_path,user_output_folder_path) #Running Genome Assembly Pipeline with Shovill and Quast
    run_prokka(user_output_folder_path)
    run_microbeannotator(user_output_folder_path)
    run_comparative_genomics(user_output_folder_path)
    tar_output(user_output_folder_path)

main()


