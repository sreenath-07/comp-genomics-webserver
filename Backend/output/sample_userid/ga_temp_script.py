#!/usr/bin/env python3

import os
for isolate_id in ["CGT1225"]:
    shovill_command = f"shovill --gsize 5M --cpus 4 --trim --outdir ./ga/{isolate_id} --R1 ../../../sample_inputs/{isolate_id}/{isolate_id}_1.fq.gz --R2 ../../../sample_inputs/{isolate_id}/{isolate_id}_2.fq.gz --force"
    copy_command = f"cp ./ga/{isolate_id}/contigs.fa ./ga/{isolate_id}/{isolate_id}_contigs.fa"
    move_command = f"mv ./ga/{isolate_id}/{isolate_id}_contigs.fa ../final_outputs/{isolate_id}_contigs.fa"
    print(shovill_command)
    print(copy_command)
    print(move_command)

    os.system(shovill_command)
    os.system(copy_command)
    os.system(move_command)
