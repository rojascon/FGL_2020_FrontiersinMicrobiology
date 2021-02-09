#################################################################################
#
#               Anaerobic Microbial Communities in a Stratified Sulfidic Lake
#                      
#              Rojas et al 2021.Organic electron donors and terminal electron 
#       acceptors structure anaerobic microbial communities and interactions in 
#                     a permanently stratified sulfidic lake
#
#                               By: Connie Rojas
#                               Created: 4 Aug 2020
#                            Last updated: 8 Feb 2021
################################################################################

## CODE FOR: running constrained correspondence analysis (CCA) on 
# microbial communities from watercolumn to determine which geochemical
# parameters are important (WATER COLUMN data)

source(file="scripts/00_background.R"); #load necessary packages and specifications

################################################################################
#             1.  Load OTU table and sample metadata
#           Chloroplast, Mitochondria, and Unknown were removed
################################################################################
taxa=read.table("data/fgl_otu.txt", sep="\t", header=T);
meta_data=read.table("data/fgl_meta.txt", sep="\t", header=T);


################################################################################
#             1.  Prepare data for CCA analysis 
################################################################################
#subset metadata to parameters of interest
ecopar=meta_data[,c(1,11,13,16,19)];
rownames(ecopar)=ecopar$sample; 
ecopar$sample=NULL;
ecopar=ecopar[order(row.names(ecopar)),];

#calculate OTU abundances from raw counts
otu_cca<-taxa[, c(6, 7:ncol(taxa))];
colnames(otu_cca)[1]="Taxa";
otu_cca=aggregate(.~Taxa, otu_cca, sum);
otu_cca[,-1] <- lapply(otu_cca[,-1], function(x) (x/sum(x))*100);

#transpose table for analysis
rownames(otu_cca)=otu_cca$Taxa;
otu_cca$Taxa=NULL;
otu_cca_t=as.data.frame(t(otu_cca));
otu_cca_t=otu_cca_t[order(row.names(otu_cca_t)),];


################################################################################
#             2.  Run constrained correspondence analysis (CCA)
################################################################################
fgl.cca <- cca(otu_cca_t~ sulfide_z+DOC+ammonia+methane_x, 
               data=ecopar)
print(summary(fgl.cca))
print(anova.cca(fgl.cca, by="terms"));
print(anova.cca(fgl.cca, by="axis"));


################################################################################
#             3.  visualize output from CCA and save plot
################################################################################
pdf("images/03_WC_cca_plot.pdf", width =7, height = 6)
plot(fgl.cca, 
     xlim=c(-3,3), ylim=c(-3,3),
     xlab="CCA1 (49%)",
     ylab="CCA2 (10%)",
     display=c("ob","cn","wa")); 

#save the plot
dev.off();

