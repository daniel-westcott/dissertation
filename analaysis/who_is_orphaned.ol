#I need to know how many genes of CzCORE are not in other regreencuts
#Load sequences
# this is an ortholang script with the sole purpose of comparing results
# Load all result lists:
#All zofingiensis proteins
Cz_aa = load_faa "/home/lang/ortholang/lab-data/Czofingiensis/Czofingiensis_461_v5.2.3.2.protein.fa"
#mmseqs
mFS = load_list "mmseqs/rgc_mmseqs_list/mmseqs_rgc_FS.txt"
mRBH = load_list "mmseqs/rgc_mmseqs_list/mmseqs_rgc_rbh.txt"
mKC = load_list "mmseqs/rgc_mmseqs_list/mmseqs_rgc_krisCut.txt"
m-10 = load_list "mmseqs/rgc_mmseqs_list/mmseqs_rgc_1waye10.txt"
m-20 = load_list "mmseqs/rgc_mmseqs_list/mmseqs_rgc_1waye20.txt"
m-50 = load_list "mmseqs/rgc_mmseqs_list/mmseqs_rgc_1waye50.txt"
m-100 = load_list "mmseqs/rgc_mmseqs_list/mmseqs_rgc_1waye100.txt"
m-200 = load_list "mmseqs/rgc_mmseqs_list/mmseqs_rgc_1waye200.txt"
# blast
bFS = load_list "blast/rgc_blast_list/blast_rgc_FS.txt"
bRBH = load_list "blast/rgc_blast_list/blast_rgc_rbh.txt"
bKC = load_list "blast/rgc_blast_list/blast_rgc_krisCut.txt"
b-10 = load_list "blast/rgc_blast_list/blast_rgc_1waye-10.txt"
b-20 = load_list "blast/rgc_blast_list/blast_rgc_1waye-20.txt"
b-50 = load_list "blast/rgc_blast_list/blast_rgc_1waye-50.txt"
b-100 = load_list "blast/rgc_blast_list/blast_rgc_1waye-100.txt"
b-200 = load_list "blast/rgc_blast_list/blast_rgc_1waye-200.txt"
# diamond_blast
dFS = load_list "diamond/rgc_diamond_list/diamond_rgc_FS.txt"
dRBH = load_list "diamond/rgc_diamond_list/diamond_rgc_rbh.txt"
dKC = load_list "diamond/rgc_diamond_list/diamond_rgc_krisCut.txt"
d-10 = load_list "diamond/rgc_diamond_list/diamond_rgc_1waye-10.txt"
d-20 = load_list "diamond/rgc_diamond_list/diamond_rgc_1waye-20.txt"
d-50 = load_list "diamond/rgc_diamond_list/diamond_rgc_1waye-50.txt"
d-100 = load_list "diamond/rgc_diamond_list/diamond_rgc_1waye-100.txt"
d-200 = load_list "diamond/rgc_diamond_list/diamond_rgc_1waye-200.txt"
#kris-cut non-PS IDS
mKC-nonPS = load_list "/home/lang/ortholang/my-data/dw_ortholang_scripts/october/rgc/rgc_output_only/mmseqs_KC_non-ps-ids.txt"
bKC-nonPS = load_list "/home/lang/ortholang/my-data/dw_ortholang_scripts/october/rgc/rgc_output_only/blast_KC_non-ps-ids.txt"
dKC-nonPS = load_list "/home/lang/ortholang/my-data/dw_ortholang_scripts/october/rgc/rgc_output_only/diamond_KC_non-ps-ids.txt"
bg-ids = load_list "/home/lang/ortholang/my-data/dw_ortholang_scripts/october/rgc/rgc_output_only/bg-ids.txt"


# all sets (&)
all-200 = b-200 & m-200 & d-200
fs_all = bFS & dFS & mFS
rbh_all = mRBH & bRBH & dRBH
kc_all = mKC & dKC & bKC
all-10 = b-10 & m-10 & d-10
all-50 = b-50 & m-50 & d-50
all-100 = b-100 & m-100 & d-100
all-200 = b-200 & d-200 & m-200

# any sets (|) 
any-200 = b-200 | m-200 | d-200
fs_any = bFS | dFS | mFS
rbh_any = mRBH | bRBH | dRBH
kc_any = mKC | dKC | bKC
any-10 = b-10 | m-10 | d-10
any-50 = b-50 | m-50 | d-50
any-100 = b-100 | m-100 | d-100
any-200 = b-200 | d-200 | m-200
non-ps-core = bg-ids ~ kc_any

#any RGC combined

big-rgc = rbh_any | fs_any | any-10 | any-50 | any-100 |any-200
orphans =  kc_any ~ big-rgc

#venndiagrams
v-200 = venndiagram [b-200, m-200, d-200]
v-100 = venndiagram [b-100, m-100, d-100]
v-50 = venndiagram [b-50, m-50, d-50]
v-20 = venndiagram [b-20, m-20, d-20]
v-10 = venndiagram [b-10, m-10, d-10]
v-rbh = venndiagram [mRBH, bRBH, dRBH]
v-kc = venndiagram [mKC, dKC, bKC]
v-fs = venndiagram [mFS, dFS, bFS]
all_venn = venndiagram [all-10, kc_all, fs_all, rbh_all, all-200]
any_venn = venndiagram [any-10, kc_any, fs_any, rbh_any, any-200]
all_small_venn = venndiagram [all-10,kc_all,fs_all, rbh_all]
any_small_venn = venndiagram [any-10,kc_any,fs_any,rbh_any]
venn-non-ps = venndiagram [bg-ids, mKC-nonPS, bKC-nonPS, dKC-nonPS]
venn-rbh-fs-core-bgids = venndiagram [bg-ids, fs_any, rbh_any, kc_any]
venn-rbh-fs-core-non-ps = venndiagram [bg-ids, fs_any, rbh_any, non-ps-core]
venn-fs-core-bgids-non-ps-core = venndiagram [bg-ids,fs_any, non-ps-core, kc_any]
venn-rbh-core-bgids-non-ps-core = venndiagram [bg-ids,rbh_any,non-ps-core,kc_any]

big-venn = venndiagram [ bg-ids, kc_any, non-ps-core, big-rgc]
orphan-venn = venndiagram [kc_any, non-ps-core, big-rgc, orphans]
#New_faa
kc_any_faa = extract_seqs Cz_aa kc_any
non_ps_faa = extract_seqs Cz_aa non-ps-core

#_______________________________________#

#bg_degs
bg_degs = load_faa "/home/lang/ortholang/my-data/genelists/20201005_bg_degs.faa"
#Arabidopsis genes that are coexpressed with CURT1, VIPP1, MGD1, and cpSAR1
at_germ = load_faa "/home/lang/ortholang/my-data/genelists/CURT_VIPP_MGD_CPSAR_coexp.faa"
#C.re acetate requiring collection. 371 genes. 
arc = load_faa "/home/lang/ortholang/my-data/genelists/ARC_Seq.fasta"

#search Cz -> At
b_bg_at = blastp 1e-5 bg_degs at_germ
d_bg_at = diamond_blastp 1e-5 bg_degs at_germ
m_bg_at = mmseqs_search 1e-5 bg_degs at_germ

#search Cz -> arc
b_bg_cr = blastp 1e-5 bg_degs arc
d_bg_cr = diamond_blastp 1e-5 bg_degs arc
m_bg_cr = mmseqs_search 1e-5 bg_degs arc

#Ids for venns

b_bg_at_id = extract_queries b_bg_at
d_bg_at_id = extract_queries d_bg_at
m_bg_at_id = extract_queries m_bg_at
b_bg_cr_id = extract_queries b_bg_cr
d_bg_cr_id = extract_queries d_bg_cr
m_bg_cr_id = extract_queries m_bg_cr

at_ids = b_bg_at_id | (d_bg_at_id | m_bg_at_id)
cr_ids = b_bg_cr_id | (d_bg_cr_id | m_bg_cr_id)

#venndiagram
core_venn = venndiagram [bg-ids, kc_any, at_ids,cr_ids]

#make 1 BHT for each At and arc
bg_at = concat_bht [b_bg_at, d_bg_at, m_bg_at]
bg_cr = concat_bht [b_bg_cr, d_bg_cr, m_bg_cr]


result = big-rgc
