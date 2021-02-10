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
# ENRICHMENTS data at the bacterial family taxonomic level

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load OTU table and sample metadata
#           Chloroplast, Mitochondria, and Unknown were removed
################################################################################
taxa=read.table("data/enrichments_otu.txt", sep="\t", header=T);
meta_data=read.table("data/enrichments_meta.txt", sep="\t", header=T);

meta_data$CarbonSource<-factor(meta_data$CarbonSource, 
                               levels=c("NoC", "methane","acetate",
                                        "propionate","butyrate","lignin",
                                        "chitin","cellulose"));


################################################################################
#         2. Calculate abundances and format data for composition bubbleplots                 
################################################################################
#select your taxonomic level
#1-Kingdom 2-Phylum 3-Class 4-Order 5-Family 6-Genus
taxa2<-taxa[, c(5, 7:ncol(taxa))];
colnames(taxa2)[1]="Taxa";

#calculate taxa relative abundances 
taxa_rel=aggregate(.~Taxa, taxa2, sum);
taxa_rel[,-1] <- lapply(taxa_rel[,-1], function(x) (x/sum(x))*100); 
#print(colSums(taxa_rel[-1])); 

#select taxa >0.23% relative abundance
taxa_rel$AVG=rowMeans(taxa_rel[-1]);
taxa_rel=taxa_rel[taxa_rel$AVG>0.23,]; taxa_rel$AVG=NULL;
rownames(taxa_rel)=taxa_rel$Taxa; taxa_rel$Taxa=NULL;

#rearrange data for ggplot
taxon=rownames(taxa_rel); rtaxon=rep(taxon,ncol(taxa_rel));
taxa_rel<-taxa_rel%>% gather(sample, abun,
                             SP.Fe.Su.a.1:SP.SO4.Su.P.2);
taxa_rel$taxon=rtaxon;

#merge taxa abundance table with metadata
taxa_meta=merge(taxa_rel, meta_data, by="sample"); 
phylum_labels=taxa[, c(2,5)];
colnames(phylum_labels)[2]="taxon";

#sort data frame by phylum, so bacterial families are organized by their phylum
taxa_meta2=merge(taxa_meta, phylum_labels, by="taxon");
taxa_meta2=taxa_meta2[order(taxa_meta2$Phylum, decreasing=TRUE),];
my.order=unique(taxa_meta2$taxon);
taxa_meta2$taxon=factor(taxa_meta2$taxon, levels=my.order);


################################################################################
#             3. Plot Bacterial Family bubbleplot and save                 
################################################################################
#select color-palette for bubbles
new_col=c("#66c2a5", "#e78ac3", "#fdc086");

#generate bubble plot
#samples are grouped by carbon source and color-coded by electron source
bubble=ggplot(data=taxa_meta2)+ 
  geom_point(mapping=aes(x = ElectronAcceptor, y = taxon, 
                         size = abun, colour=ElectronAcceptor),
             stat="identity")+
  facet_grid(facets=.~CarbonSource,
             scales="free_x")+ 
  scale_colour_manual(values=new_col)+
  labs(y="",x="",colour="Electron Acceptor",
       size="Relative Abundance (%)")+
  theme_bw()+
  theme(legend.position="bottom",
        legend.title=element_text(size=12, color="black",face="bold"),
        legend.text=element_text(size=11),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=10),
        axis.ticks.x=element_blank(),
        strip.text = element_text(size =11),
        panel.background = element_rect(size=2));

#plot(bubble);

##save image 
ggsave(filename="05_ENRICH_bubble_family.pdf",
       device="pdf",path="./images",
       plot=bubble,
       width=11,
       height=7.5,
       units="in",
       dpi=500);
