library(RColorBrewer)
library(vcfR) #read vcf
library(pheatmap)
library(gplots)
library(tidyverse)


rm(list=ls())
#####################################
####Get genotype and make heatmap####
#####################################
#load vcf file 
vcf<-read.vcfR("./walleye.crh8.inv.recode.vcf")

#extract genotypes in numeric format
genotypes_numeric<-extract.gt(vcf, element = "GT", mask = FALSE, as.numeric = FALSE,
                                   return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE,
                                   convertNA = TRUE)


#function to convert genotypes to 0,1,2 format used for heatmap
convertGenos<-function(genotype){
  alleles<-as.numeric(str_split_fixed(genotype,"/",2))
  genoCode<-alleles[1]+alleles[2]
  return(genoCode)
}

#convert to 0,1,2 format
genosCode<-apply(genotypes_numeric,1:2,convertGenos)
genosCode<-t(genosCode)

#format for heatmap2
heatmap2<-genosCode
heatmap2[is.na(heatmap2)] <- -1 #replace msising info with -1
heatmap2<-apply(heatmap2,1:2,function(x) as.numeric(x))

###fancy heatmap
coord_class <- read.table("./chr8_inv_pca_cluster.txt",header=TRUE)
class <- coord_class %>% 
  select(sample=name,class=cluster)

class <- class %>% 
  arrange(class)
heatmap2 <- heatmap2 %>% 
  as.data.frame() %>% 
  rownames_to_column("sample")

heatmap2 <- left_join(class, heatmap2) %>% 
  select(-class) %>% 
  column_to_rownames("sample") %>% 
  as.matrix()
heatmap2<-apply(heatmap2,1:2,function(x) as.numeric(x))




pop_colors <- list(pop=brewer.pal(6,"Set2"))
pop_colors
names(pop_colors$pop)=unique(popmap$pop)

class=data.frame(class,row.names=1)
class
class$class <- factor(class$class, levels=c(0,1,2))
class_colors <- list(class=c("red", "purple", "blue"))
names(class_colors$class)=c(0,1,2)



pdf("./chr8_inv_genotype.pheatmap.pdf", width = 12, height = 9)
pheatmap(heatmap2,
         annotation_row = class,
         annotation_colors  = class_colors,
         cluster_rows = F, 
         cluster_cols = F,
         show_colnames = F,
         show_rownames = FALSE,
         legend_breaks = c(-1,0,1,2),
         legend_labels =c("NA","homo2","het","homo1"))
dev.off()
