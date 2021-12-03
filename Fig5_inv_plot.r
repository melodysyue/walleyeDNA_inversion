library(tidyverse)
library(gridExtra)
library(ggsignif)
library(lostruct)
library(grid)
library(ggsci)
library(patchwork)

rm(list=ls())
loci <- read.table("./walleye.loci.list", header=F)
colnames(loci) <- c("chromosome", "position")
popmap <- read.table("./combined_popmap2_edited.txt", header=FALSE)
colnames(popmap) <- c("sample", "pop")
popmap$pop <- factor(popmap$pop, levels=c("Kansas", "Swan_Creek", "Dauphin_River", "Matheson", "Red_River"))
unique(popmap$pop)

#pop colors
popcol <- c("#fdae61", "#d73027", "#abd9e9", "#2c7bb6", "#313695")
names(popcol)=c("Kansas", "Swan_Creek", "Dauphin_River", "Matheson", "Red_River")
popcol


c <- "NC_041338.1"
start <- 15260000
end <- 15900000
mat <- read.table(paste0(c,"_lostruct_input_matrix.txt"),header=T)
colnames(mat)=popmap$sample
p <- loci %>% 
    filter(chromosome==c & position >=start & position <=end) %>% 
    pull(position)
mat_sub <- mat[rownames(mat)%in%p,]
  
  
  out <- cov_pca(mat_sub, k=2) # out is a numeric vector
  PC1_perc <- out[2]/out[1]
  PC2_perc <- out[3]/out[1]
  samples <- colnames(mat_sub)
  matrix.out <- t(matrix(out[4:length(out)], ncol=length(samples), byrow=T))
  out <- as_tibble(cbind(samples, matrix.out)) %>% 
    rename(name=samples, PC1="V2", PC2="V3") %>% 
    mutate(PC1=as.double(PC1), PC2=as.double(PC2))
  
  #kmeans to cluster along PC1 and PCA
  try_3_clusters <-try(kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1]))))
  if("try-error" %in% class(try_3_clusters)){
    kmeans_cluster <-kmeans(matrix.out[,1], 2, centers=c(min(matrix.out[,1]), max(matrix.out[,1])))
  }else{
    kmeans_cluster <- kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1])))
  }
  out$cluster <- kmeans_cluster$cluster - 1
  out$cluster <- as.character(out$cluster)
  betweenSS_perc <- kmeans_cluster$betweenss / kmeans_cluster$totss
  
  mat_sub %>% as.tibble() -> snps
  snps_geno <- snps %>% gather("name","genotype",1:ncol(snps)) %>% 
    group_by(name, genotype) %>%
    summarize(count=n()) %>%
    spread(genotype, count)
  
  snps_geno[is.na(snps_geno)]=0
  
  #heterozygosity <- snps_geno %>%   
  #  summarize(het=`1`/(`0` + `1` + `2`)) #this would cause errors if there is no '2' genotype
  
  heterozygosity <- snps_geno %>% 
    ungroup() %>% 
    mutate(total=rowSums(dplyr::select(., -name))) %>% 
    group_by(name) %>% 
    summarize(het= `1`/ total) 
  
  het <- inner_join(out, heterozygosity)
  het.cluster <- het %>% 
    group_by(cluster) %>% 
    summarize(mean=mean(het)) 
  
  out <- left_join(out, popmap, by=c("name"="sample"))
  
  ###PCA plot by pop
  
  pca_plot_pop <- out %>%
    ggplot(., aes(x=PC1, y=PC2, col=pop)) + 
    geom_point(size = 4, alpha =0.7)+
    theme_bw(base_size = 20) +
    scale_color_manual(values=popcol, name="Population", 
                       labels=c("Cedar Bluff", "Swan Creek", 
                                "Dauphin River", "Matheson Island", "Red River"))+
    xlab(paste0("PC",1," (",round(100*PC1_perc,2),"%)")) +
    ylab(paste0("PC",2," (",round(100*PC2_perc,2),"%)")) +
    theme(legend.position = "bottom")+
    guides(colour=guide_legend(nrow=2))
  

  #LD heatmap zoom
  ld<- read.table("NC_041338.1.ld", header=T)
  
  ld_plot <- ld %>% 
    ggplot(aes(x=BP_A/1e6,y=BP_B/1e6,color=R2))+
    geom_point()+
    xlab("SNP 1 (Position in Mb)")+
    ylab("SNP 2 (Position in Mb)")+
    geom_vline(xintercept=start/1e6, size=1, linetype="dotted")+
    geom_vline(xintercept=end/1e6, size=1, linetype="dotted")+
    scale_color_gradient2(low="snow", high="red", name=expression(italic(r^2)), 
                          limits=c(0,1),
                          breaks=seq(0,1,0.25))+
    theme_classic(base_size = 20)+
    coord_cartesian(xlim=c(start/1e6-1, end/1e6+1), ylim=c(start/1e6-1, end/1e6+1))+
    theme(legend.position = c(0.9,0.4))
  

  #PCA: 3 clusters or not
  pca_plot_cluster <- out %>%
    ggplot(., aes(x=PC1, y=PC2, col=cluster)) + 
    geom_point(size = 4, alpha =0.5)+
    theme_bw(base_size = 20) +
    scale_color_manual(name="Cluster", values=c("red","purple","blue")) +
    xlab(paste0("PC",1," (",round(100*PC1_perc,2),"%)")) +
    ylab(paste0("PC",2," (",round(100*PC2_perc,2),"%)")) +
    annotate(geom="text", -Inf, Inf,
             label=paste0("Discreteness=",round(betweenSS_perc, digits=4)), hjust=-0.2, vjust=2, size=6)+
    theme(legend.position = "none")
  
  
  #boxplot to compare heterozygosity among clusters
  het_plot <- het %>% 
    ggplot(., aes(x=as.character(cluster), y=het, fill=as.character(cluster))) + 
    geom_boxplot() + 
    theme_bw(base_size = 20) +
    scale_fill_manual(name="Cluster", values=c("red","purple","blue")) +
    xlab("Cluster") + ylab("Heterozygosity") +
    geom_signif(comparisons = list(c("1", "0"), c("1","2"), c("0","2")), 
                map_signif_level=TRUE)+
    theme(legend.position = "none")
  
  ### Genotype frequency
popcluster <- out %>% 
    group_by(pop, cluster) %>% 
    summarise(n=n()) %>% 
    group_by(pop) %>% 
    mutate(total=sum(n)) %>% 
    mutate(freq=n / total)
  popcluster$cluster <- as.character(popcluster$cluster)
  
  freq_plot <-  popcluster %>% 
    ggplot(aes(x=fct_rev(pop), y=freq, fill=cluster)) +
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(name="cluster", values=c("red","purple","blue")) +
    xlab("")+
    ylab("Frequency")+
    scale_x_discrete(expand=c(0,0), labels=c("Red River", "Matheson Island", "Dauphin River",
                                             "Swan Creek", "Cedar Bluff"))+
    scale_y_continuous(expand=c(0,0))+
    theme_classic(base_size=20)+
    theme(legend.position = "none")+
    coord_flip()
  
  
#write.table(out, "./chr8_inv_pca_cluster.txt", quote=F, row.names = F)
  
  
  #output figures
  
  pdf("chr8_inv_summary.pdf", width = 16, height = 12)
  
  p1 <- pca_plot_pop + ld_plot +
    plot_layout(widths = c(2,1))

  p <- p1 / (pca_plot_cluster | het_plot | freq_plot ) + 
    plot_annotation(tag_levels = "A",
                    title=paste0("Chromosome 8 (Position ", start, " - ", end, ")"),
                    theme=theme(plot.title = element_text(size=20)))
  p
  dev.off()

  
  





