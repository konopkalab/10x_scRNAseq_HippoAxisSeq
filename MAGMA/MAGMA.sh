
# Rcode to get the gene anno
#gtf <- rtracklayer::import('gencode.v19.annotation.gtf')
#gtf_df <- as.data.frame(gtf)
#gtf_df <- gtf_df[gtf_df$type == "gene" & gtf_df$transcript_type == "protein_coding",]
#gtf_df$chr = gsub("chr", "", gtf_df$seqnames)
#gtf_df <- gtf_df[!(gtf_df$chr %in% c("X","Y","M"),]
#out = data.frame(ENSG=gtf_df$transcript_name, CHR=gtf_df$chr, START=gtf_df$start, STOP=gtf_df$end)
#out <- out[!(duplicated(out$ENSG)),]
#write.table(out, file="gencode.v19.genes.out", quote=F, row.names = F, col.names = F)


# Run MAGMA
# Annotate with magma 
#magma --annotate --gene-loc /U3/stefano/src/magma_v1.07/GENES/gencode.v19.genes.out --snp-loc /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur.bim --out /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation

# ADHD
# PMID: 29325848
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/adhd_jul2017.bed N=55374 \
--out ADHD_2017

# ALZHEIMER
# PMID: https://www.nature.com/articles/s41588-018-0311-9
# awk -v OFS="\t" '{print $6,$8}' ALZ_sumstats_Jansenetal.txt > ALZ_JENSEN_ForMagma.bed 
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/ALZ_JENSEN_ForMagma.bed N=455253 \
--out ALZ_2019

# AUTISM  
# https://www.biorxiv.org/content/early/2017/11/27/224774
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/iPSYCH-PGC_ASD_Nov2017.bed N=46350 \
--out ASD_2017

# BIPOLAR DISORDER
# PMID: 29906448
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/BDvsCONT.sumstats.bed N=74194 \
--out BIP_2018

# EPILEPSY
# PMID: 30531953
#awk -v OFS="\t" '{print $3,$10}' all_epilepsy_METAL.txt | tail -n +2  > EPILEPSY_ForMagma.bed 
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/EPILEPSY_ForMagma.bed N=44889 \
--out EPILEPSY_2018

# MAJOR DEPRESSIVE DISORDER
#PMID: 29700475
# awk -v OFS="\t" '{print $2,$11}' MDD2018_ex23andMe.bed | tail -n +2 > MDD2018_ex23andMe_ForMagma.bed 
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/MDD2018_ex23andMe_ForMagma.bed N=480359 \
--out MDD_2018

# SCHIZOPHRENIA
# PMID: 29906448
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/CLOZUK_PGC2.bed N=105318 \
--out SZ_2018

# EDUCATION ATTAINMENT
# PMID: 21694764
# awk -v OFS="\t" '{print $1,$9}' GWAS_EA_excl23andMe.txt | tail -n +2  > GWAS_EA_excl23andMe_ForMagma.bed 
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/GWAS_EA_excl23andMe_ForMagma.bed N=1100000  \
--out tEduATT

# General COGNITIVE
# PMID: 29844566
# awk -v OFS="\t" '{print $3,$6}' Davies2018_OPEN_DATASET_summary_results.txt | tail -n +2  > Davies2018_OPEN_DATASET_summary_ForMagma.bed
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/Davies2018_OPEN_DATASET_summary_ForMagma.bed N=300486  \
--out tCognFun

# INTELLIGENCE 
# PMID: 29942086
# awk -v OFS="\t" '{print $1,$11}' SavageJansen_2018_intelligence_metaanalysis.txt | tail -n +2 > SavageJansen_2018_intelligence_ForMagma.bed  
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/SavageJansen_2018_intelligence_ForMagma.bed N=269867  \
--out tIntelligence

# CORONARY HEART DISEASE
# PMID: 21378990
# awk -v OFS="\t" '{print $1,$6}' CARDIoGRAM_GWAS_RESULTS.txt | tail -n +2  > CARDIOGRAM_ForMagma.bed 
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/CARDIOGRAM_ForMagma.bed N=84509  \
--out zCHD

# DIABETES
# PMID: 22885922
# awk -v OFS="\t" '{print $1,$6}' DIAGRAMv3.2012DEC17.txt | tail -n +2  > DIAGRAMv3_ForMagma.bed 
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/DIAGRAMv3_ForMagma.bed N=149821  \
--out zDIAB

# OSTEOPOROSIS
# PMID: 22504420
# awk -v OFS="\t" '{print $1,$4}' GEFOS2_FNBMD_POOLED_GC.txt  | tail -n +2  > GEFOS2_ForMagma.bed 
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/GEFOS2_ForMagma.bed N=153377  \
--out zOSTEO

# HEIGHT
# https://www.biorxiv.org/content/early/2018/07/09/355057
# awk -v OFS="\t" '{print $3,$9}' Meta-analysis_Wood_et_al+UKBiobank_2018 | tail -n +2  > HEIGHT_ForMagma.bed 
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/HEIGHT_ForMagma.bed N=589762  \
--out zHGT

# BMI
# PMID: 30108127
# awk -v OFS="\t" '{print $1,$10}' BMI-GERA+GIANT-2018.tsv | tail -n +2 > BMI_ForMagma.bed 
magma --bfile /U3/stefano/src/magma_v1.07/BKG/g1000_eur/g1000_eur \
--gene-annot /U3/stefano/src/magma_v1.07/GENES/hg19_gencodeV19_Window10kb_magma_annotation.genes.annot \
--set-annot *.txt col=1,2 \
--pval /U5/Stefano/GWAS_DATABASE/BMI_ForMagma.bed N=527927  \
--out zBMI

R CMD BATCH --vanilla GetDataTable.R
























