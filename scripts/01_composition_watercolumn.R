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
# WATER COLUMN data at the bacterial phylum, order, and genus taxonomic levels

source(file="scripts/00_background.R"); #load necessary packages and specifications

################################################################################
#             1.  Load OTU table and sample metadata
#           Chloroplast, Mitochondria, and Unknown were removed
################################################################################
taxa=read.table("data/fgl_otu.txt", sep="\t", header=T);
meta_data=read.table("data/fgl_meta.txt", sep="\t", header=T);


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

#create plot (group samples by water depth)
barphy=ggplot(data=pbar, 
              mapping=aes(x=sample,y=abun, fill=Taxa))+
  geom_bar(stat = "identity")+
  facet_grid(~depth, scales="free_x")+ 
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="Bacterial Phylum")+
  scale_fill_manual(values=phy_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_text(size =14, face="bold"),
        panel.background = element_rect(size=2))+
  guides(fill=guide_legend(ncol=1));

#plot(barphy);

##save image 
ggsave(filename="01_WC_stacked_phyla.pdf",
       device="pdf",path="./images",
       plot=barphy,
       width=8,
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
          "#117777", "#737373", "#44AA77", "#88CCAA", "#777711", "#AAAA44", 
          "#DDDD77","#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455");

#create plot (group samples by water depth)
barord=ggplot(data=obar, 
              mapping=aes(x=sample,y=abun, fill=Taxa))+
  geom_bar(stat = "identity")+
  facet_grid(~depth, scales="free_x")+ 
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="Bacterial Order")+
  scale_fill_manual(values=ord_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_text(size =14, face="bold"),
        panel.background = element_rect(size=2))+
  guides(fill=guide_legend(ncol=1));

#plot(barord);

##save image 
ggsave(filename="01_WC_stacked_order.pdf",
       device="pdf",path="./images",
       plot=barord,
       width=8,
       height=5,
       units="in",
       dpi=500);


################################################################################
#             4. Create Genus level composition barplots                 
################################################################################
#select your taxonomic level
#1-Kingdom 2-Phylum 3-Class 4-Order 5-Family 6-Genus
taxa2<-taxa[, c(6, 7:ncol(taxa))];
colnames(taxa2)[1]="Taxa";

#calculate taxa relative abundances 
taxa_rel=aggregate(.~Taxa, taxa2, sum);
taxa_rel[,-1] <- lapply(taxa_rel[,-1], function(x) (x/sum(x))*100);
#print(colSums(taxa_rel[-1]));

#keep taxa >0.7% relative abundance across samples
taxa_rel$AVG=rowMeans(taxa_rel[,-1]);
taxa_rel=taxa_rel[taxa_rel$AVG>0.7,];
taxa_rel$AVG=NULL;

#denote the rest of taxa as "Other"
newrow=c(NA, 100-colSums(taxa_rel[2:ncol(taxa_rel)])); 
taxa_rel=rbind(taxa_rel, newrow); 
taxa_rel$Taxa=as.character(taxa_rel$Taxa);
taxa_rel[nrow(taxa_rel),1]="Other";

#melt data frame for ggplot
gbar<-reshape2::melt(taxa_rel, id.vars="Taxa",value.name = "abun");
colnames(gbar)[2]="sample";
gbar=merge(gbar, meta_data, by="sample");

#set color palette
gen_col=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", 
          "#737373","#117744","#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", 
          "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","black","grey")

#create plot (group samples by water depth)
bargen=ggplot(data=gbar, 
              mapping=aes(x=sample,y=abun, fill=Taxa))+
  geom_bar(stat = "identity")+
  facet_grid(~depth, scales="free_x")+ 
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="Bacterial Genus")+
  scale_fill_manual(values=gen_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_text(size =14, face="bold"),
        panel.background = element_rect(size=2))+
  guides(fill=guide_legend(ncol=1));

#plot(bargen);

##save image 
ggsave(filename="01_WC_stacked_gen.pdf",
       device="pdf",path="./images",
       plot=bargen,
       width=6,
       height=6,
       units="in",
       dpi=500);


################################################################################
#             5. Download table of taxa abundances                 
################################################################################
#calculate abundances
taxa2<-taxa[, c(6, 7:ncol(taxa))];
colnames(taxa2)[1]="Taxa";
taxa_rel=aggregate(.~Taxa, taxa2, sum);
taxa_rel[,-1] <- lapply(taxa_rel[,-1], function(x) (x/sum(x))*100);

#add sample depth labels
colnames(taxa_rel)=meta_data$depth[match(colnames(taxa_rel),
        meta_data$sample)];
colnames(taxa_rel)[1]=1;
taxa_rel=taxa_rel[, order(names(taxa_rel))];
colnames(taxa_rel)[1]="Taxa";

#write file
#write.csv(taxa_rel, file="data/01_wc_genus.csv",row.names=F);
