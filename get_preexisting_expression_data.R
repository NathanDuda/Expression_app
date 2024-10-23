


# import all expression files 
p1_exp <- read.csv("C:/Users/17735/Downloads/Azenta_Normalization_Results/P1_Averaged_RPKM_exp.tsv", sep="")
p2_exp <- read.csv("C:/Users/17735/Downloads/Azenta_Normalization_Results/P2_Averaged_RPKM_exp.tsv", sep="")
p3_exp <- read.csv("C:/Users/17735/Downloads/Azenta_Normalization_Results/P3_Averaged_RPKM_exp.tsv", sep="")
p4_exp <- read.csv("C:/Users/17735/Downloads/Azenta_Normalization_Results/P4_Averaged_RPKM_exp.tsv", sep="")

# merge all expression files into one
exp <- merge(p1_exp, p2_exp, by = 'ID')
exp <- merge(exp, p3_exp, by = 'ID')
exp <- merge(exp, p4_exp, by = 'ID')

# format column names of expression file
cnames <- colnames(exp)
cnames <- gsub("\\.", "_", cnames)
cnames <- gsub("_+$", "", cnames)
cnames <- gsub("Counts_", "", cnames)
cnames <- gsub("_0", "", cnames)
cnames <- gsub("2003_P9", "AB_2003_P9", cnames)
cnames <- gsub("1810_P9", "AB_1810_P9", cnames)
colnames(exp) <- cnames

# use gene symbols instead of ENSG ID
ENSG_symbol <- read.csv("C:/Users/17735/Downloads/Azenta_Analyses/PreExisting_Data/protein_coding_ENSG_ids.tsv", sep="")

exp <- exp %>%
  select(ENSG = ID, everything()) %>%
  merge(ENSG_symbol, by = 'ENSG') %>%
  select(ID = gn_symbol, everything(), -ENSG)


write.table(exp, file = 'C:/Users/17735/Downloads/Expression_app/preexisting_data/Expression_RPKM.tsv')




