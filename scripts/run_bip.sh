dir=$(pwd)
mkdir outputs

#Default prediction
#Select highly supported proteins (>= 2 protein-specific peptides from one tissue
python3 "$dir"/scripts/find_highly_supp_prot.py \
        "$dir"/inputs/default_pred.gtf \
        "$dir"/inputs/mapped_default \
        "$dir"/outputs/hsp_default.gtf
#Find high-scoring proteins with a BRAKER transcript score cutoff
python3 "$dir"/scripts/find_and_apply_score_filter.py \
        "$dir"/inputs/default_pred.gtf \
        "$dir"/outputs/hsp_default.gtf \
        "$dir"/inputs/tx_scores_default.tsv \
        "$dir"/outputs/hsp_cut_default.gtf

#Relaxed prediction
#Find gene-specific peptides from genes with at least two gene-specific peptides
python3 "$dir"/scripts/find_gene_spec_pept.py \
        "$dir"/inputs/relaxed_pred.gtf \
        "$dir"/inputs/mapped_relaxed \
        "$dir"/outputs/gene_spec_pept_relaxed.tsv
#Select supported proteins
python3 "$dir"/scripts/select_supp_prot.py \
        "$dir"/inputs/relaxed_pred.gtf \
        "$dir"/inputs/tx_scores_relaxed.tsv \
        "$dir"/outputs/gene_spec_pept_relaxed.tsv \
        "$dir"/outputs/sp_relaxed.gtf
#Find high-scoring proteins with a BRAKER transcript score cutoff
python3 "$dir"/scripts/find_and_apply_score_filter.py \
        "$dir"/inputs/relaxed_pred.gtf \
        "$dir"/outputs/sp_relaxed.gtf \
        "$dir"/inputs/tx_scores_relaxed.tsv \
        "$dir"/outputs/sp_cut_relaxed.gtf

#Combine supported and high-scoring proteins from both predictions
python3 "$dir"/scripts/unite_gtf.py \
        "$dir"/outputs/sp_cut_relaxed.gtf \
        "$dir"/outputs/hsp_cut_default.gtf \
        relaxed \
        default \
        20 \
        "$dir"/outputs/bip.gtf
