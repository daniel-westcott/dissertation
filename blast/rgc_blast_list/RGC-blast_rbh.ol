#ReGreenCut - a phylogenenomic cut to identify genes related to thylakoid biogenesis

#Load sequences
#All zofingiensis proteins
Cz_aa = load_faa "/home/lang/ortholang/lab-data/Czofingiensis/Czofingiensis_461_v5.2.3.2.protein.fa"
#Arabidopsis proteome
athal_aa = load_faa "/home/lang/ortholang/lab-data/Athaliana/Athaliana_447_Araport11.protein.fa"
#Only genes involved in re-greening after glucose has been removed
bg_degs = load_faa "/home/lang/ortholang/my-data/genelists/20201005_bg_degs.faa"
#Arabidopsis genes that are coexpressed with CURT1, VIPP1, MGD1, and cpSAR1
at_germ = load_faa "/home/lang/ortholang/my-data/genelists/CURT_VIPP_MGD_CPSAR_coexp.faa"
#Auxenochlorella protothecoides. Can also metabolize glucose. Unfortunately I don't have a transcriptome of glucose grown Aux.
aux_aa = load_faa "/home/lang/ortholang/lab-data/Auxeprot1/Auxeprot1_GeneCatalog_proteins_20170909.aa.fasta"
#Chlamy to compare for photosynthetic genes
cre_aa = load_faa "/home/lang/ortholang/lab-data/Creinhardtii/Creinhardtii_281_v5.6.protein.fa"
#E. coli, S.cer, and melainabacteria to remove non-photosynthetic genes
e_coli = load_faa "/home/lang/ortholang/my-data/genelists/E.coli_0157-H7.UP000000558_83334.fasta"
S_cer = load_faa "/home/lang/ortholang/my-data/genelists/S.cerevisiae_ATCC204508_S288c_UP000002311_559292.fasta"
#Melainabacteria have a few strains that I'm going to concatenate
mela_aa = load_faa "/home/lang/ortholang/lab-data/melainabacteria/soo_origins_2017/MEL_A1.faa"
melb_aa = load_faa "/home/lang/ortholang/lab-data/melainabacteria/soo_origins_2017/MEL_B1.faa"
melc_aa = load_faa "/home/lang/ortholang/lab-data/melainabacteria/soo_origins_2017/MEL_C2.faa"
mel_aa = concat_faa [mela_aa,melb_aa,melc_aa]
#C.re acetate requiring collection. 371 genes or so. 
arc_seq_aa = load_faa "/home/lang/ortholang/my-data/genelists/ARC_Seq.fasta"

#reciprocal blast searches F = A to B and R = B to A
e_col-degsF = blastp 1.0e-5 e_coli bg_degs
e_col-degsR = blastp 1.0e-5 bg_degs e_coli
s_cer-degsF = blastp 1.0e-5 S_cer bg_degs
s_cer-degsR = blastp 1.0e-5 bg_degs S_cer
mel-degsF = blastp 1.0e-5 mel_aa bg_degs
mel-degsR = blastp 1.0e-5 bg_degs mel_aa
atgerm-degs-blastF = blastp 1.0e-5 at_germ bg_degs
atgerm-degs-blastR = blastp 1.0e-5 bg_degs at_germ
aux-degs-blastF = blastp 1.0e-5 aux_aa bg_degs
aux-degs-blastR = blastp 1.0e-5 bg_degs aux_aa
cre-degs-blastF = blastp 1.0e-5 cre_aa bg_degs
cre-degs-blastR = blastp 1.0e-5 bg_degs cre_aa
arc_seq-degs-blastF = blastp 1.0e-5 arc_seq_aa bg_degs
arc_seq-degs-blastR = blastp 1.0e-5 bg_degs arc_seq_aa
ath-degs-blastF = blastp 1.0e-5 bg_degs athal_aa
ath-degs-blastR = blastp 1.0e-5 athal_aa bg_degs

#use reciprocal_best to finish the RBH
e_col-degs = reciprocal_best e_col-degsF e_col-degsR
s_cer-degs = reciprocal_best s_cer-degsF s_cer-degsR
mel-degs = reciprocal_best mel-degsF mel-degsR
atgerm-degs-blast = reciprocal_best atgerm-degs-blastF atgerm-degs-blastR
aux-degs-blast = reciprocal_best aux-degs-blastF aux-degs-blastR
cre-degs-blast = reciprocal_best cre-degs-blastF cre-degs-blastR
arc-degs-blast = reciprocal_best arc_seq-degs-blastF arc_seq-degs-blastR
ath-degs-blast = reciprocal_best ath-degs-blastF ath-degs-blastR

#Create BHTs of organism sets
non-ps-bht = concat_bht [e_col-degs, mel-degs, s_cer-degs]
ps-bht = concat_bht [aux-degs-blast, cre-degs-blast,ath-degs-blast]
biogenic-bht = concat_bht [atgerm-degs-blast,arc-degs-blast]

#Extract IDs for comparisons
rbh-bg-at-olap = extract_targets atgerm-degs-blast
rbh-bg-aux-olap = extract_targets aux-degs-blast
rbh-bg-cre-olap = extract_targets cre-degs-blast
rbh-bg-scer-olap = extract_targets s_cer-degs
rbh-bg-ecol-olap = extract_targets e_col-degs
rbh-bg-mel-olap = extract_targets mel-degs
rbh-bg-arc-olap = extract_targets arc-degs-blast
rbh-bg-ath-olap = extract_queries ath-degs-blast

rbh-bg-ids = extract_ids bg_degs
non-ps-ids = extract_targets non-ps-bht
ps-bht-ids = extract_targets ps-bht
biogenic-bht-ids = extract_targets biogenic-bht

#Set operations to create cuts
#Non-photynthetic are genes that have high homology to those in yeast, e.coli, and melainaibacteria
non-photosynthetic = rbh-bg-scer-olap | rbh-bg-ecol-olap | rbh-bg-mel-olap
#This is the set of all genes that have high homology to the Cz rbh-bg-ids set
all_lap = rbh-bg-at-olap | rbh-bg-aux-olap | rbh-bg-cre-olap | non-photosynthetic | rbh-bg-arc-olap
#These are the genes that don't match anything
nolap = rbh-bg-ids ~ all_lap
#biogenic set includes arc-seq and at_germ
rbh-biogenic = rbh-bg-at-olap | rbh-bg-arc-olap
#regreencut is the overlapping genes with the non-photosynthetic subracted. 
regreencut_blast_rbh = rbh-biogenic ~ non-photosynthetic
#A venn of all photosynthetic genes overlapping with the bg set you need to add a full At blast set here.
venn_ps_blast_rbh = venndiagram [rbh-bg-cre-olap,rbh-bg-aux-olap,rbh-bg-ids]
#A venn of non-photosynthetic genes overlapping with the rbh-bg-set
venn_non_ps_blast_rbh = venndiagram [rbh-bg-scer-olap, rbh-bg-ecol-olap, rbh-bg-mel-olap, rbh-bg-ids]
#A broad level venn diagram showing the original set, non-ps set, photosynthetic, and biogenic
venn_all_blast_rbh = venndiagram [non-ps-ids, ps-bht-ids,rbh-bg-at-olap,rbh-bg-arc-olap, rbh-bg-ids]

tabl_blast_rbh = sets_table [rbh-bg-ids,
                   rbh-bg-at-olap,rbh-bg-arc-olap,
                   rbh-bg-aux-olap,
                   rbh-bg-cre-olap,
                   rbh-bg-scer-olap,
                   rbh-bg-ecol-olap,
                   rbh-bg-mel-olap,rbh-bg-ath-olap]
result = regreencut_blast_rbh
