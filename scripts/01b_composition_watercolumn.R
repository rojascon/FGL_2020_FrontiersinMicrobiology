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

## CODE FOR: generating bubbleplots of microbial community composition for 
# WATER COLUMN data at the bacterial class taxonomic level

source(file="scripts/00_background.R"); #load necessary packages and specifications

################################################################################
#             1.  Load OTU table and sample metadata
#           Chloroplast, Mitochondria, and Unknown were removed
################################################################################
taxa=read.table("data/fgl_otu.txt", sep="\t", header=T);
meta_data=read.table("data/fgl_meta.txt", sep="\t", header=T);


################################################################################
#             2. Create Class level composition bubbleplots                 
################################################################################
#select your taxonomic level
#1-Kingdom 2-Phylum 3-Class 4-Order 5-Family 6-Genus
taxa2<-taxa[, c(3, 7:ncol(taxa))];
colnames(taxa2)[1]="Taxa";

#calculate taxa relative abundances 
taxa_rel=aggregate(.~Taxa, taxa2, sum);
taxa_rel[,-1] <- lapply(taxa_rel[,-1], function(x) (x/sum(x))*100); 
#print(colSums(taxa_rel[-1])); 

#select taxa >0.5% relative abundance
taxa_rel$AVG=rowMeans(taxa_rel[-1]);
taxa_rel=taxa_rel[taxa_rel$AVG>0.5,]; taxa_rel$AVG=NULL;
rownames(taxa_rel)=taxa_rel$Taxa; taxa_rel$Taxa=NULL;

#rearrange data for ggplot
taxon=rownames(taxa_rel); rtaxon=rep(taxon,ncol(taxa_rel));
taxa_rel<-taxa_rel%>% gather(sample, abun,
                             CR2:CR7);
taxa_rel$taxon=rtaxon;

#merge taxa abundance table with metadata
taxa_meta=merge(taxa_rel, meta_data, by="sample");
taxa_meta$depthf=factor(taxa_meta$depth);
taxa_meta$depthf=forcats::fct_rev(factor(taxa_meta$depthf));

#rename archaea order so that it appears first in plot
taxa_meta$taxon[taxa_meta$taxon=="Nanoarchaeia"]="1Nanoarchaeia";

#generate bubble plot
bubble=ggplot(data=taxa_meta)+ 
  geom_point(mapping=aes(x = depthf, y = taxon, 
                         size = abun),
             stat="identity")+
  coord_flip()+
  labs(y="",x="",colour="",
       size="Relative Abundance (%)")+
  theme_bw()+
  theme(panel.background = element_rect(colour = "black", size=1),
        legend.position="bottom",
        legend.title=element_text(size=10, color="black",face="bold"),
        legend.text=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10, angle = 90, vjust = 0.66),
        axis.ticks.x=element_blank());

#plot(bubble);

##save image 
ggsave(filename="01_WC_bubble_order.pdf",
       device="pdf",path="./images",
       plot=bubble,
       width=4,
       height=4.5,
       units="in",
       dpi=500);