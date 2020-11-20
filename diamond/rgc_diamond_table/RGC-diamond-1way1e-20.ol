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

#Permissive 1 way diamond searches.
e_col-degs = diamond_blastp 1.0e-5 bg_degs e_coli
s_cer-degs = diamond_blastp 1.0e-5 bg_degs S_cer
mel-degs = diamond_blastp 1.0e-5 bg_degs mel_aa
atgerm-degs = diamond_blastp 1.0e-5 bg_degs at_germ
aux-degs = diamond_blastp 1.0e-5 bg_degs aux_aa
cre-degs = diamond_blastp 1.0e-5 bg_degs cre_aa
arc_seq-degs = diamond_blastp 1.0e-5 bg_degs arc_seq_aa
arc-cre-degs = diamond_blastp 1.0e-5 arc_seq_aa cre_aa
ath-degs = diamond_blastp 1.0e-5 bg_degs athal_aa

#Increase stringency 
atgerm-degs-strict = filter_evalue 1.0e-20 atgerm-degs
aux-degs-strict = filter_evalue 1.0e-20 aux-degs
e_col-degs-strict = filter_evalue 1.0e-20 e_col-degs
s_cer-degs-strict = filter_evalue 1.0e-20 s_cer-degs
mel-degs-strict = filter_evalue 1.0e-20 mel-degs
arc-degs-strict = filter_evalue 1.0e-20 arc_seq-degs
cre-degs-strict = filter_evalue 1.0e-20 cre-degs
ath-degs-strict = filter_evalue 1.0e-20 ath-degs

#Create BHTs of organism sets
non-ps-bht = concat_bht [e_col-degs-strict, mel-degs-strict, s_cer-degs-strict]
ps-bht = concat_bht [aux-degs-strict, cre-degs-strict, ath-degs-strict]
biogenic-bht = concat_bht [atgerm-degs-strict,arc-degs-strict]

#Extract IDs for comparisons
bg-at-olap-e-20 = extract_queries atgerm-degs-strict
bg-aux-olap-e-20 = extract_queries aux-degs-strict
bg-cre-olap-e-20 = extract_queries cre-degs-strict
bg-scer-olap-e-20 = extract_queries s_cer-degs-strict
bg-ecol-olap-e-20 = extract_queries e_col-degs-strict
bg-mel-olap-e-20 = extract_queries mel-degs-strict
bg-arc-olap-e-20 = extract_queries arc-degs-strict
bg-ath-olap-e-20 = extract_queries ath-degs-strict
bg-ids = extract_ids bg_degs
non-ps-ids = extract_queries non-ps-bht
ps-bht-ids = extract_queries ps-bht
biogenic-bht-ids = extract_queries biogenic-bht

#Set operations to create cuts
#Non-photynthetic are genes that have high homology to those in yeast, e.coli, and melainaibacteria
non-photosynthetic = bg-scer-olap-e-20 | bg-ecol-olap-e-20 | bg-mel-olap-e-20
#This is the set of all genes that have high homology to the Cz bg-ids set
all_lap = bg-at-olap-e-20 | bg-aux-olap-e-20 | bg-cre-olap-e-20 | non-photosynthetic | bg-arc-olap-e-20
#These are the genes that don't match anything
nolap-e-20 = bg-ids ~ all_lap
# biogenic is arc-seq and at_germ together
biogenic-e-20 = bg-at-olap-e-20 | bg-arc-olap-e-20
#regreencut is the overlapping genes with the non-photosynthetic subracted. 
regreencut_diamond_blastp_eVal-20 = biogenic-e-20 ~ non-photosynthetic
#A venn of all photosynthetic genes overlapping with the bg set you need to add a full At diamond set here.
venn_ps_diamond_blastp_eVal-20 = venndiagram [bg-cre-olap-e-20,bg-aux-olap-e-20,bg-ids]
#A venn of non-photosynthetic genes overlapping with the bg-set
venn_non_ps_diamond_blastp_eVal-20 = venndiagram [bg-scer-olap-e-20, bg-ecol-olap-e-20, bg-mel-olap-e-20, bg-ids]
#A broad level venn diagram showing the original set, non-ps set, photosynthetic, and biogenic
venn_all_1way_diamond_blastp_eVal-20 = venndiagram [non-ps-ids, ps-bht-ids,bg-at-olap-e-20,bg-arc-olap-e-20,bg-ids]

tabl_1way_diamond_blastp_e-20 = sets_table [bg-ids,
                   bg-at-olap-e-20,bg-arc-olap-e-20,
                   bg-aux-olap-e-20,
                   bg-cre-olap-e-20,
                   bg-scer-olap-e-20,
                   bg-ecol-olap-e-20,
                   bg-mel-olap-e-20,bg-ath-olap-e-20]
result = tabl_1way_diamond_blastp_e-20
