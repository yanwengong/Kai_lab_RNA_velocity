## generate pdf plot for BRCA vs normal breast tissue, epi cells
## DDRtree cell color by group, arrow by velocity
## Ddrtree pesudotime 
## Ddrtree + scenery 
## Tsne + velocity , 1k from each individual 

## libraries
library(monocle)
library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(velocyto.R)
library(pagoda2)
library(Matrix)
library(igraph)
library(rlist)

## functions


## set up path
input_path <- "/Users/yanwengong/Documents/kai_lab/rna_velocity_2019_summer/data/brca_nomal_epi"
plot_output_path <- "/Users/yanwengong/Documents/kai_lab/rna_velocity_2019_summer/plot/brca_nomal_epi"
## read in files

## load RData
load('normal_brca_6samples.RData')

## load seurat object
brca_normal_seurat_object<-readRDS(paste(input_path, "brca_normal_seurat_object", sep = "/")) #1290

## get the 1k cell for each 6 samples 
set.seed(20190715)
sample_cell_each1k <- brca_normal_seurat_object@meta.data %>%tibble::rownames_to_column("cell_barcode") %>% 
  filter(assig.label != "UnClass") %>%
  group_by(group) %>% sample_n(1000) %>% data.frame() %>% arrange(cell_barcode) %>% select(cell_barcode) 

## Ddrtree pesudotime
nature_comm_list <- read.delim(paste(input_path, "Nat_Comm_Droplet_Ordering_Genes.txt", sep = "/"), header = TRUE) #214

# Extract data, phenotype data, and feature data from the SeuratObject
data <- brca_normal_seurat_object@raw.data[, c(sample_cell_each1k$cell_barcode)]
meta_sub <- brca_normal_seurat_object@meta.data %>% subset(rownames(brca_normal_seurat_object@meta.data)%in%c(sample_cell_each1k$cell_barcode))
meta_sub_order<- meta_sub[order(row.names(meta_sub)),]

pd <- new('AnnotatedDataFrame', data =meta_sub_order)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

# Construct monocle cds
monocle_cds_each1k <- newCellDataSet(data,
                                     phenoData = pd,
                                     featureData = fd,
                                     lowerDetectionLimit = 0.5,
                                     expressionFamily = negbinomial.size())
# estimate size factors and dispersions
monocle_cds_each1k <- estimateSizeFactors(monocle_cds_each1k)
monocle_cds_each1k <- estimateDispersions(monocle_cds_each1k) #this functions use way too much memory

# select genes
nature_comm_genes <- nature_comm_list$gene

# each 1k, nature comm genes

monocle_cds_each1k_ddrtree_nc <- setOrderingFilter(monocle_cds_each1k, ordering_genes = unique(as.character(nature_comm_genes)))
monocle_cds_each1k_ddrtree_nc <- reduceDimension(monocle_cds_each1k_ddrtree_nc, method = 'DDRTree')
monocle_cds_each1k_ddrtree_nc <- orderCells(monocle_cds_each1k_ddrtree_nc)
#png(paste(plot_output_path, "ddrtree_each1k_naturecom_colorLabel.png", sep = "/"), width = 500, height = 350)
#plot_cell_trajectory(monocle_cds_each1k_ddrtree_nc, color_by = "assig.label", cell_size = 0.5)
#dev.off()

# colored by pseudotime
plot_cell_trajectory(monocle_cds_each1k_ddrtree_nc, color_by = "State")
monocle_cds_ddrtree <- orderCells(monocle_cds_each1k_ddrtree_nc, root_state = 3)
pdf(paste(plot_output_path, "ddrtree_each1k_naturecom_pseudotime.pdf", sep = "/"), width = 5, height = 4)
plot_cell_trajectory(monocle_cds_ddrtree, color_by = "Pseudotime", cell_size = 0.5)
dev.off()

# save the ddrtree location
ddrtree_loc_1kCell_eachSample<- monocle_cds_each1k_ddrtree_nc@reducedDimS %>% data.frame() %>% t() %>% data.frame() 
colnames(ddrtree_loc_1kCell_eachSample) <- c("ddrtree1", "ddrtree2")
ddrtree_loc_1kCell_eachSample$Row <- rownames(ddrtree_loc_1kCell_eachSample)
rownames(ddrtree_loc_1kCell_eachSample) <- c()
ddrtree_loc_1kCell_eachSample$Row <- gsub("_", ":", gsub("..$", "x", ddrtree_loc_1kCell_eachSample$Row))
ddrtree_loc_1kCell_eachSample <- ddrtree_loc_1kCell_eachSample[, c(3,1,2)]
#write.table(brca_normal_1keach_ddrtree_nat_com, file = "/Users/yanwengong/Documents/kai_lab/rna_velocity_2019_summer/data/loc_files/brca_normal_1keach_ddrtree_nat_com.txt", sep = "\t",
#            row.names = TRUE, col.names = TRUE)

## Ddrtree + scenery
scenergy <- read.delim(paste(input_path, "scEnergy_NormalBrca_Quy.txt",  sep = "/"))
# rename the cell name in scenergy
scenergy$Row <- gsub("_", ":", gsub("..$", "x", scenergy$Row))
# filter the cells based on subset file 
scenergy_ddrtree_loc <- scenergy%>% inner_join(ddrtree_loc_1kCell_eachSample, by = "Row")
# plot ddrtree color code with scEnergy
mid<-median(scenergy_ddrtree_loc$scEnergy)
pdf(paste(plot_output_path, "ddrtree_each1k_scEnergeryColor_blueredDark.pdf", sep = "/"), width = 5, height = 4)
scenergy_ddrtree_loc %>% ggplot(aes(x = ddrtree1, y = ddrtree2, color = scEnergy)) + 
  geom_point(size = 0.5) + scale_color_gradient(low="blue", high="red")+ theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

#scenergy_ddrtree_loc <- scenergy_ddrtree_loc %>% arrange(scEnergy)
#scnergy_color_rampallette <- colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(scenergy_ddrtree_loc$scEnergy))
#names(scnergy_color_rampallette) <- as.character(scenergy_sort$Row)
#scenergy_ddrtree_loc %>% ggplot(aes(x = ddrtree1, y = ddrtree2, color = scEnergy)) + 
#  geom_point(size = 0.5) + 
#  scale_colour_gradientn(colours=rainbow(2))+ theme_bw() 

pdf(paste(plot_output_path, "ddrtree_each1k_scEnergeryColor_bluered.pdf", sep = "/"), width = 5, height = 4)
scenergy_ddrtree_loc %>% ggplot(aes(x = ddrtree1, y = ddrtree2)) + 
  geom_point(aes(colour = scEnergy), size = 0.5, alpha = 1) + scale_colour_gradient2(low = "blue",
                                                high = "red", mid = "white", midpoint = mid)+ theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

####### involce velocity ####################

## functions

## this is function that build pagoda object from 
## input: the filtered spliced data
## output: the pagoda object r with pca reduction 

build_pagoda <- function (input){
  ## make gene names unique by appending sequence numbers to duplicates
  rownames(input) <- make.unique(rownames((input)))
  
  ## NOTE: i need to do this step here not in filter_cell_gene function to avoid error in velocity calculation 
  input <- input %>% as.matrix() %>% as("dgCMatrix")
  
  ## create pagoda2 object; 
  r <- Pagoda2$new(input,modelType='plain',trim=10,log.scale=T)
  ## adjust the variance, to normalize the extent to which genes with very
  ## different expression will contribute to the downstream analysos
  r$adjustVariance(plot=T,do.par=T,gam.k=10)
  
  ## reduce the dataset dimentsions by running PCA
  ## nPcs number of PCs\n
  ## n.odgenes whether a certain number of top overdispersed genes should be used
  ## maxit maximum number of irlba iterations to use
  ## used r package irlba
  set.seed(0)
  r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
  
  ## make a knn graph
  ## k use k neighbors
  ## use cosine distance A*B/ (|A|*|B|)
  set.seed(1)
  r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine') ## optional for PCA?
  
  ## calculate clusters based on the KNN graph for PCA
  set.seed(2)
  #r$getKnnClusters(method=infomap.community,type='PCA') ## does not work again
  # object 'infomap.community' not found
  
  return(r)
}

assign_cell_group <- function(cell_list){
  for (i in 1:nrow(cell_list)){
    cell_i <- cell_list[i,1];
    group_i <- cell_list[i,2];
    names(group_i)<-cell_i
    
    if (i == 1){
      group <- group_i;
    } else {
      group <- c(group, group_i);
    }
  }
  group <- as.factor(group)
  
  return(group)
}

calculate_velocity <- function(emat, nmat, r_object, cluster_label, fit_quantile){
  ## keep cells that in pagoda object
  emat <- emat[,rownames(r_object$counts)]; 
  nmat <- nmat[,rownames(r_object$counts)]; 
  
  ## calculate cell-cell distance 
  cell_dist <- as.dist(1-armaCor(t(r_object$reductions$PCA))) 
  
  ## filter gene ## cluster label
  emat <- filter.genes.by.cluster.expression(emat,cluster_label,min.max.cluster.average = 0.5)
  nmat <- filter.genes.by.cluster.expression(nmat,cluster_label,min.max.cluster.average = 0.05)
  intersect_genes <- length(intersect(rownames(emat),rownames(nmat))) #2333
  print(paste("number of intersect genes between splice and unsplice: ", intersect_genes, sep = ""))
  
  ## calculate velocity
  rvel <- velocyto.R::gene.relative.velocity.estimates(emat,nmat,kCells=25,deltaT=1,cell.dist=cell_dist,fit.quantile=fit_quantile)
  return(rvel)
}

generate_loc_file <- function(low_dim_file){
  #low_dim_file$X <- gsub("..$", "x", gsub("_", ":", low_dim_file$X)) 
  #low_dim_file <- low_dim_file%>% dplyr::filter(X %in% cell_group_file$X)
  ## match with cell group
  loc <- as.matrix(apply(low_dim_file[,-1],2,as.numeric))
  row.names(loc) <- low_dim_file[,1]
  return(loc)
} 

## read in cell group info
cell_iden <- read.csv(file=paste(input_path, "cell_barcode.csv", sep = "/"), header=TRUE, sep=",")
## match the row name in loom file
cell_iden$X <- gsub("..$", "x", gsub("_", ":", cell_iden$X))


## filter emat based on cell_list and build pangoda object
# nat_com_ddrtree_loc <-  read.delim(paste(input_path, "brca_normal_1keach_ddrtree_nat_com.txt", sep="/"), header = TRUE) 
# ddrtree_loc_1kCell_eachSample<-nat_com_ddrtree_loc
# ddrtree_loc_1kCell_eachSample$Row <- rownames(ddrtree_loc_1kCell_eachSample)
# rownames(ddrtree_loc_1kCell_eachSample) <- c()
# ddrtree_loc_1kCell_eachSample$Row <- gsub("_", ":", gsub("..$", "x", ddrtree_loc_1kCell_eachSample$Row))
# ddrtree_loc_1kCell_eachSample <- ddrtree_loc_1kCell_eachSample[, c(3,1,2)]


emat_c_filtered_1keach <- emat_c%>%as.matrix()%>%as.data.frame()%>%select(one_of(ddrtree_loc_1kCell_eachSample$Row)) #33538 genes x 6k cells
r_1keach <- build_pagoda(emat_c_filtered_1keach)

ddrtree_loc_1kCell_eachSample_loc <- generate_loc_file(ddrtree_loc_1kCell_eachSample)
cell_group_1keach <- cell_iden %>% filter(X %in% ddrtree_loc_1kCell_eachSample$Row) %>%
  select(X, assig.label) %>% assign_cell_group()

## velocity calculation
rvel_1keach <- calculate_velocity(emat_c, nmat_c, r_1keach, cell_group_1keach, fit_quantile = 0.02)

## generate cell color
cell.colors_1keach <- pagoda2:::fac2col(cell_group_1keach)
##TODO: update loc_file 

par(mfrow=c(1,1))
pdf(paste(plot_output_path, "velocity_ddrtree_1keach_nat_com.pdf", sep = "/"), width = 5, height = 4.5)
show.velocity.on.embedding.cor(ddrtree_loc_1kCell_eachSample_loc,rvel_1keach,n=200,scale='sqrt',cell.colors=ac(cell.colors_1keach,alpha=0.8),cex=0.4,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.001, xlab = "Component 1", ylab = "Component 2")
dev.off()


## velocity on t-SNE 
## read in cell location file
tsne_loc <- read.csv(file=paste(input_path, "tsne_df.csv", sep = "/"), header=TRUE, sep=",")
tsne_loc$X <- gsub("_", ":", gsub("..$", "x", tsne_loc$X))
tsne_loc <- generate_loc_file(tsne_loc)

par(mfrow=c(1,1))
pdf(paste(plot_output_path, "velocity_tsne_1keach_cellGroupColor.pdf", sep = "/"), width = 5, height = 4.5)
show.velocity.on.embedding.cor(tsne_loc,rvel_1keach,n=200,scale='sqrt',cell.colors=ac(cell.colors_1keach,alpha=0.8),cex=0.4,arrow.scale=7,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.001, xlab = "Component 1", ylab = "Component 2")
dev.off()

## plot velocity with scEnergy as color
scnergy_color_rampallette <- colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(scenergy_ddrtree_loc$scEnergy))
scenergy_sort <- scenergy_ddrtree_loc %>% arrange(scEnergy)
names(scnergy_color_rampallette) <- as.character(scenergy_sort$Row)

pdf(paste(plot_output_path, "velocity_tsne_1keach_scEnergyColor.pdf", sep = "/"), width = 5, height = 4)
show.velocity.on.embedding.cor(tsne_loc,rvel_1keach,n=200,scale='sqrt',cell.colors=ac(scnergy_color_rampallette,alpha=0.8),cex=0.4,arrow.scale=7,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.001, xlab = "Component 1", ylab = "Component 2")
dev.off()

