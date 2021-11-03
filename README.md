# BIP
Pipeline for improving BRAKER2 gene predictions with MS/MS data.

# Description



# Instructions

1.  Make a directory X. Make a directory ```scripts``` and ```inputs``` inside the directory X. The directory “scripts” should contain the following scripts: ```find_highly_supp_prot.py, find_and_apply_score_filter.py, find_gene_spec_pept.py, make_tx_scores_tsv.py, unite_gtf.py, select_supp_prot.py, run_bip.sh```. The directory ```inputs``` should contain the following files without headers:   
(a) Two ```.gtf``` files with the default and relaxed BRAKER2 predictions. The files must be named ```default_pred.gtf``` and ```relaxed_pred.gtf```.  
(b) Two ```.tsv``` files with transcript BRAKER2 scores. There should be following columns in the ```.tsv``` files: 1) protein id (=transcriptid from ```.gtf``` files); 2) BRAKER2 transcript score. The files should be named ```tx_scores_default.tsv``` and ```tx_scores_relaxed.tsv```. If a ```.gtf``` file contains transcript scores, this file can be produced by running: 
```    
    python3 /path/to/directory_X/scripts/make_tx_scores_tsv.py\/path/to/directory_X/inputs/pred_file.gtf\/path/to/directoryX/inputs/output_file.tsv
```  
(c) Two directories containing ```.tsv``` files with peptide mapping data. Each ```.tsv``` file stands for one tissue and should be named accordingly. There should be following columns in the ```.tsv``` files: 1) peptide sequence, 2) protein id (=transcriptid from ```.gtf``` files), 3) ```+``` if peptide is unique (found in one protein), ```-``` if not. These two directories should be named ```mapped_default``` and ```mapped_relaxed```. 

2. Go to the directory X and run the pipeline.
```
    bash /path/to/directory_X/scripts/run_bip.sh  
```
3. The final file with the improved prediction is named ```bip.gtf```.

# Contributors
Author: Kateryna Neishsalo
Supervisors: Prof. Dmitrij Frishman, Prof. Mathias Wilhelm, Dr. Nils Rugen.
Additional support from: Prof. Bernhard Küster, Prof. Mark Borodovsky, Dr. Tomáš Brůna.
