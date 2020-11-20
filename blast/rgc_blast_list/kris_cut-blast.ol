#Kris Cut is a simplified version of the ReGreenCut that is just the BG-IDs with the non-photosynthetic set removed
#The initial searches will be permissive

#Load sequences

#Only genes involved in re-greening after glucose has been removed
bg_degs = load_faa "/home/lang/ortholang/my-data/genelists/20201005_bg_degs.faa"

#E. coli, S.cer, and melainabacteria to remove non-photosynthetic genes
e_coli = load_faa "/home/lang/ortholang/my-data/genelists/E.coli_0157-H7.UP000000558_83334.fasta"
S_cer = load_faa "/home/lang/ortholang/my-data/genelists/S.cerevisiae_ATCC204508_S288c_UP000002311_559292.fasta"
mela_aa = load_faa "/home/lang/ortholang/lab-data/melainabacteria/soo_origins_2017/MEL_A1.faa"
melb_aa = load_faa "/home/lang/ortholang/lab-data/melainabacteria/soo_origins_2017/MEL_B1.faa"
melc_aa = load_faa "/home/lang/ortholang/lab-data/melainabacteria/soo_origins_2017/MEL_C2.faa"
mel_aa = concat_faa [mela_aa,melb_aa,melc_aa]

#Permissive 1 way blast searches.
e_col-degs = blastp 1.0e-5 bg_degs e_coli
s_cer-degs = blastp 1.0e-5 bg_degs S_cer
mel-degs = blastp 1.0e-5 bg_degs mel_aa

#Create non-photosynthetic set
non-ps-bht = concat_bht [e_col-degs, mel-degs, s_cer-degs]

#Extract IDs for comparisons
#Individual sets for venns
e_col-degs-ids = extract_queries e_col-degs
s_cer-degs-ids = extract_queries s_cer-degs
mel-degs-ids = extract_queries mel-degs

#Cz Biogenic ids
bg-ids = extract_ids bg_degs
#all the concatentated non-photosynthetic ids
non-ps-ids = e_col-degs-ids | s_cer-degs-ids | mel-degs-ids

#Set operations to create cuts

#kris kut
kris_kut_blast = bg-ids ~ non-ps-ids
#genes retained in melainabacteria but not in any other set
mel_cut = mel-degs-ids ~ (e_col-degs-ids| s_cer-degs-ids)
#A broad level venn diagram showing the original set, non-ps set, photosynthetic, and biogenic
vennKris_kut_blast = venndiagram [bg-ids,e_col-degs-ids,s_cer-degs-ids,mel-degs-ids]

tabl_1way_blast_kris_cut = sets_table [bg-ids,
                   e_col-degs-ids,
                   s_cer-degs-ids,
                   mel-degs-ids]
result = non-ps-ids
