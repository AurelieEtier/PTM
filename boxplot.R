library(ggpubr)
library(dplyr)

args<-commandArgs(TRUE)

file = args[1]
outDir = args[2]

df = read.table(file, sep="\t", header=TRUE)

yPos = max(df$Length) + 2
p = ggboxplot(df, x = "Condition", y = "Length", color = "Condition", add = "jitter", shape = "Condition", ylab = "Peptide length (amino acids)") + stat_compare_means(method = "wilcox.test", label.y = yPos)
print(paste(outDir, "/peptide_size.png", sep=""))
ggsave(paste(outDir, "/peptide_size.png", sep=""))