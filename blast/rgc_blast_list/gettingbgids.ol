#ReGreenCut - a phylogenenomic cut to identify genes related to thylakoid biogenesis
#Adusting the previous eValues because they resulted in an empty set. 
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

#Permissive 1 way blast searches.
e_col-degs = blastp 1.0e-5 bg_degs e_coli
s_cer-degs = blastp 1.0e-5 bg_degs S_cer
mel-degs = blastp 1.0e-5 bg_degs mel_aa
atgerm-degs = blastp 1.0e-5 bg_degs at_germ
aux-degs = blastp 1.0e-5 bg_degs aux_aa
cre-degs = blastp 1.0e-5 bg_degs cre_aa
arc_seq-degs = blastp 1.0e-5 bg_degs arc_seq_aa
arc-cre-degs = blastp 1.0e-5 arc_seq_aa cre_aa
ath-degs = blastp 1.0e-5 bg_degs athal_aa

#Evolutionarily distant organisms will be less strict. Closely related will be strict.  
atgerm-degs-strict = filter_evalue 1.0e-30 atgerm-degs
aux-degs-strict = filter_evalue 1.0e-75 aux-degs
cre-degs-strict = filter_evalue 1.0e-75 cre-degs
e_col-degs-strict = filter_evalue 1.0e-10 e_col-degs
s_cer-degs-strict = filter_evalue 1.0e-10 s_cer-degs
mel-degs-strict = filter_evalue 1.0e-50 mel-degs
arc-degs-strict = filter_evalue 1.0e-75 arc_seq-degs
ath-degs-strict = filter_evalue 1.0e-30 ath-degs

#Create BHTs of organism sets
non-ps-bht = concat_bht [e_col-degs-strict, mel-degs-strict, s_cer-degs-strict]
ps-bht = concat_bht [aux-degs-strict, cre-degs-strict,ath-degs-strict]
biogenic-bht = concat_bht [atgerm-degs-strict,arc-degs-strict]

#Extract IDs for comparisons
bg-at-olap_FS = extract_queries atgerm-degs-strict
bg-aux-olap_FS = extract_queries aux-degs-strict
bg-cre-olap_FS = extract_queries cre-degs-strict
bg-scer-olap_FS = extract_queries s_cer-degs-strict
bg-ecol-olap_FS = extract_queries e_col-degs-strict
bg-mel-olap_FS = extract_queries mel-degs
bg-arc-olap_FS = extract_queries arc-degs-strict
bg-ath-olap_FS =extract_queries ath-degs-strict
bg-ids_FS = extract_ids bg_degs
non-ps-ids = extract_queries non-ps-bht
ps-bht-ids = extract_queries ps-bht
biogenic-bht-ids = extract_queries biogenic-bht

#Set operations to create cuts
#Non-photynthetic are genes that have high homology to those in yeast, e.coli, and melainaibacteria
non-photosynthetic = bg-scer-olap_FS | bg-ecol-olap_FS | bg-mel-olap_FS
#This is the set of all genes that have high homology to the Cz bg-ids_FS set
all_lap = bg-at-olap_FS | bg-aux-olap_FS | bg-cre-olap_FS | non-photosynthetic | bg-arc-olap_FS
#These are the genes that don't match anything
nolap = bg-ids_FS ~ all_lap
#biogenic is arc-sec and at_germ combined
biogenic_FS = bg-at-olap_FS | bg-arc-olap_FS
#regreencut is the overlapping genes with the non-photosynthetic subracted. 
regreencut_blast_FS = bg-at-olap_FS ~ non-photosynthetic
#A venn of all photosynthetic genes overlapping with the bg set you need to add a full At blast set here.
venn_ps_blast_FS = venndiagram [bg-cre-olap_FS,bg-aux-olap_FS,bg-ids_FS]
#A venn of non-photosynthetic genes overlapping with the bg-set
venn_non_ps_blast_FS = venndiagram [bg-scer-olap_FS, bg-ecol-olap_FS, bg-mel-olap_FS, bg-ids_FS]
#A broad level venn diagram showing the original set, non-ps set, photosynthetic, and biogenic
venn_all_1way_blast_FS_evalFS = venndiagram [non-ps-ids, ps-bht-ids,bg-at-olap_FS,bg-arc-olap_FS, bg-ids_FS]

tabl_blast_FS = sets_table [bg-ids_FS,
                   bg-at-olap_FS,bg-arc-olap_FS,
                   bg-aux-olap_FS,
                   bg-cre-olap_FS,
                   bg-scer-olap_FS,
                   bg-ecol-olap_FS,
                   bg-mel-olap_FS,bg-ath-olap_FS]
result =bg-ids_FS
