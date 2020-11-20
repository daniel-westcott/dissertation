
#The goal of this is to find the top hits so that I can find out information about the C.re and A.thal genes that are flagged.

#Load sequences
#bg_degs
bg_degs = load_faa "/home/lang/ortholang/my-data/genelists/20201005_bg_degs.faa"
#Arabidopsis genes that are coexpressed with CURT1, VIPP1, MGD1, and cpSAR1
at_germ = load_faa "/home/lang/ortholang/my-data/genelists/CURT_VIPP_MGD_CPSAR_coexp.faa"
#C.re acetate requiring collection. 371 genes. 
arc = load_faa "/home/lang/ortholang/my-data/genelists/ARC_Seq.fasta"
Cz_core = load_list "/home/lang/ortholang/my-data/dw_ortholang_scripts/october/rgc/CzCORE_any.txt"

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

at_ids = b_bg_at_id & d_bg_at_id & m_bg_at_id
cr_ids = b_bg_cr_id & d_bg_cr_id & m_bg_cr_id
deg_ids = extract_ids bg_degs
#venndiagram
core_venn = venndiagram [deg_ids, Cz_core, at_ids,cr_ids]

#make 1 BHT for each At and arc
bg_at = concat_bht [b_bg_at, d_bg_at, m_bg_at]
bg_cr = concat_bht [b_bg_cr, d_bg_cr, m_bg_cr]

venn_test = venndiagram [(extract_queries b_bg_at), (extract_queries d_bg_at)]

result = bg_cr
