# BIP
Pipeline for improving BRAKER2 gene predictions with MS/MS data.


Instructions:

1.  Make a directory X. Make a directory “scripts” and “inputs” inside the directory X. The directory “scripts” should contain the following scripts: findhighlysuppprot.py, findandapplyscorefilter.py, findgenespecpept.py, maketxscorestsv.py, unitegtf.py, selectsuppprot.py, runbip.sh. The directory “inputs” should contain the following files without headers:
(a) Two .gtf files with the default and relaxed BRAKER2 predictions. The files must be named defaultpred.gtf and relaxedpred.gtf.
(b) Two .tsv files with transcript BRAKER2 scores. There should be following columns in the .tsv files: 1) protein id (=transcriptid from .gtf files); 2) BRAKER2 transcript score. The files should be named txscoresdefault.tsv and txscoresrelaxed.tsv. If a .gtf file contains transcript scores, this file can be produced by running: 
    python3 /path/to/directoryX/scripts/maketxscorestsv.py\/path/to/directoryX/inputs/predfile.gtf\/path/to/directoryX/inputs/outputfile.tsv
(c) Two directories containing .tsv files with peptide mapping data. Each .tsv file stands for one tissue and should be named accordingly. There should be following columns in the .tsv files: 1) peptide sequence, 2) protein id (=transcriptid from .gtf files), 3) “+” if peptide is unique (found in one protein), “-” if not. These two directories should be named mappeddefault and mappedrelaxed. 

2. Go to the directory X and run the pipeline.
    bash /path/to/directoryX/scripts/runbip.sh  
    
3. The final file with the improved prediction is named bip.gtf.
