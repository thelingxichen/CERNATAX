# CERNATAX

The workflow to detect ceRNA network mRNA-miRNA-lncRNA triplet axis from RNA expression profiles.

We have manually curated a ceRNA network database by integrating TargetSCAN 8.0, miRTarBase 9.0, miRDB 6.0, NPInter 4.0, ENCORI/starBase 2.0, miRWalk V3, and RNAInter in 2020. 

To run CERNATAX with full reference ceRNA network, please download the full db file [`ceRNA_database.csv`](https://doi.org/10.5281/zenodo.15357964), and place it into `cernatax/db'.

The tutorial and case study of the SCZ cohort is placed at `run_scz.ipynb`.
