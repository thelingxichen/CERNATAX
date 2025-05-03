# CERNATAX

The workflow to detect ceRNA network mRNA-miRNA-lncRNA triplet axis from RNA expression profiles.

We have manually curated a ceRNA network database by integrating TargetSCAN 8.0, miRTarBase 9.0, miRDB 6.0, NPInter 4.0, ENCORI/starBase 2.0, miRWalk V3, and RNAInter in 2020. 

To run CERNATAX with full reference ceRNA network, please download the full db file [`ceRNA_database.csv`](https://www.dropbox.com/scl/fi/cepr33ia5ltm54ru79j35/ceRNA_database.csv?rlkey=l5c0dc8rw3e4vr6ekylxlu7dl&dl=0), and place it into `cernatax/db'.

The tutorial and case study of the SCZ cohort is placed at `run_scz.ipynb`.
