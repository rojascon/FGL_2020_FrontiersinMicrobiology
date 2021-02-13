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

## CODE FOR: generating stacked bar plots of microbial community composition for 
# ENRICHMENTS data at the bacterial phylum, ando order taxonomic levels

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
#             2. Create Phylum level composition barplots                 
################################################################################
#select your taxonomic level
#1-Kingdom 2-Phylum 3-Class 4-Order 5-Family 6-Genus
taxa2<-taxa[, c(2, 7:ncol(taxa))];
colnames(taxa2)[1]="Taxa";

#calculate taxa relative abundances 
taxa_rel=aggregate(.~Taxa, taxa2, sum);
taxa_rel[,-1] <- lapply(taxa_rel[,-1], function(x) (x/sum(x))*100); 
#print(colSums(taxa_rel[-1]));   

#keep taxa >1% relative abundance across samples
taxa_rel$AVG=rowMeans(taxa_rel[,-1]);
taxa_rel=taxa_rel[taxa_rel$AVG>1,];
taxa_rel$AVG=NULL;

#denote the rest of taxa as "Other"
newrow=c(NA, 100-colSums(taxa_rel[2:ncol(taxa_rel)])); 
taxa_rel=rbind(taxa_rel, newrow); 
taxa_rel$Taxa=as.character(taxa_rel$Taxa);
taxa_rel[nrow(taxa_rel),1]="Other";

#melt data frame for ggplot
pbar<-reshape2::melt(taxa_rel, id.vars="Taxa",value.name = "abun");
colnames(pbar)[2]="sample";
pbar=merge(pbar, meta_data, by="sample");

#set color palette
phy_col=c('#66c2a5','#fc8d62','#7570b3','#e78ac3','#a6d854','#ffd92f',
          "#6baed6","#fc9272","#238b45","#525252","#e7298a","#08589e",
          "#bf812d");

#create plot (group samples by carbon source)
barphy=ggplot(data=pbar, 
              mapping=aes(x=sample,y=abun, fill=Taxa))+
  geom_bar(stat = "identity")+
  facet_grid(~CarbonSource, scales="free_x")+ 
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="Bacterial Phylum")+
  scale_fill_manual(values=phy_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_text(size =9, face="bold"),
        panel.background = element_rect(size=2))+
  guides(fill=guide_legend(ncol=1));

#plot(barphy);

##save image 
ggsave(filename="05_ENRICH_stacked_phyla.pdf",
       device="pdf",path="./images",
       plot=barphy,
       width=8.5,
       height=5,
       units="in",
       dpi=500);


################################################################################
#             3. Create Order level composition barplots                 
################################################################################
#select your taxonomic level
#1-Kingdom 2-Phylum 3-Class 4-Order 5-Family 6-Genus
taxa2<-taxa[, c(4, 7:ncol(taxa))];
colnames(taxa2)[1]="Taxa";

#calculate taxa relative abundances 
taxa_rel=aggregate(.~Taxa, taxa2, sum);
taxa_rel[,-1] <- lapply(taxa_rel[,-1], function(x) (x/sum(x))*100);
#print(colSums(taxa_rel[-1]));

#keep taxa >1% relative abundance across samples
taxa_rel$AVG=rowMeans(taxa_rel[,-1]);
taxa_rel=taxa_rel[taxa_rel$AVG>1,];
taxa_rel$AVG=NULL;

#denote the rest of taxa as "Other"
newrow=c(NA, 100-colSums(taxa_rel[2:ncol(taxa_rel)])); 
taxa_rel=rbind(taxa_rel, newrow); 
taxa_rel$Taxa=as.character(taxa_rel$Taxa);
taxa_rel[nrow(taxa_rel),1]="Other";

#melt data frame for ggplot
obar<-reshape2::melt(taxa_rel, id.vars="Taxa",value.name = "abun");
colnames(obar)[2]="sample";
obar=merge(obar, meta_data, by="sample");

#set color palette
ord_col=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
          "#117777", "#44AAAA", "#117744", "#88CCAA", 
          "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", 
          "#771122", "grey","#AA4455", "#DD7788", "#636363");

##create plot (group samples by carbon source)
barord=ggplot(data=obar, 
              mapping=aes(x=sample,y=abun, fill=Taxa))+
  geom_bar(stat = "identity")+
  facet_grid(~CarbonSource, scales="free_x")+ 
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="Bacterial Order")+
  scale_fill_manual(values=ord_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=9),
        legend.title = element_text(size=10, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=11, face="bold"),
        axis.text.y = element_text(size=11),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_text(size =9, face="bold"),
        panel.background = element_rect(size=2))+
  guides(fill=guide_legend(ncol=1));

#plot(barord);

##save image 
ggsave(filename="05_ENRICH_stacked_order.pdf",
       device="pdf",path="./images",
       plot=barord,
       width=9.5,
       height=6,
       units="in",
       dpi=500);


################################################################################
#             5. Download table of taxa abundances                 
################################################################################
#calculate abundances
taxa2<-taxa[, c(4, 7:ncol(taxa))];
colnames(taxa2)[1]="Taxa";
taxa_rel=aggregate(.~Taxa, taxa2, sum);
taxa_rel[,-1] <- lapply(taxa_rel[,-1], function(x) (x/sum(x))*100);

#add sample depth labels
meta_data$newname=paste(meta_data$CarbonSource, 
                        meta_data$ElectronAcceptor, 
                        sep="_")
colnames(taxa_rel)=meta_data$newname[match(colnames(taxa_rel),
                                         meta_data$sample)];
colnames(taxa_rel)[1]=1;
taxa_rel=taxa_rel[, order(names(taxa_rel))];
colnames(taxa_rel)[1]="Taxa";

#write file
#write.csv(taxa_rel, file="data/05_ENRICH_order.csv",row.names=F);
