#Fayettevielle Green Lake Manuscript
#FrontiersinMicrobiology (Rojas et al 2020)
#Script to generate stacked bar plots of microbial community composition for ENRICHMENTS

source(file="scripts/background.R") #load necessary packages and specifications

#read in your data
taxa=read.table("data/enrichments_otu_2021.txt", sep="\t", header=T)
meta_data=read.table("data/enrichments_meta.txt", sep="\t", header=T)

meta_data$CarbonSource<-factor(meta_data$CarbonSource, 
                        levels=c("NoC", "methane","acetate",
                                 "propionate","butyrate","lignin",
                                 "chitin","cellulose"))

#select your taxonomic level
#1-Kingdom 2-Phylum 3-Class 4-Order 5-Family 6-Genus
taxa<-taxa[, c(6, 7:ncol(taxa))] 
print(colnames(taxa[1])); colnames(taxa)[1]="Taxa"

#calculate proportions of bacterial taxa
taxa_rel=aggregate(.~Taxa, taxa, sum)
taxa_rel[,-1] <- lapply(taxa_rel[,-1], function(x) (x/sum(x))*100) 
print(colSums(taxa_rel[-1]))   

#keep families >1% relative abundance across samples (1% for order; 0.7 for genera)
taxa_rel$AVG=rowMeans(taxa_rel[,-1]);
taxa_rel=taxa_rel[taxa_rel$AVG>1,];
taxa_rel$AVG=NULL;

#denote the rest of phyla as "Other"
newrow=c(NA, 100-colSums(taxa_rel[2:ncol(taxa_rel)])); 
taxa_rel=rbind(taxa_rel, newrow); 
taxa_rel$Taxa=as.character(taxa_rel$Taxa);
taxa_rel[nrow(taxa_rel),1]="Other";

#melt data frame for ggplot
pbar<-reshape2::melt(taxa_rel, id.vars="Taxa",value.name = "abun");
colnames(pbar)[2]="sample";
pbar=merge(pbar, meta_data, by="sample");

#set color palette
new_col=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
          "#117777", "#44AAAA", "#117744", "#88CCAA", 
          "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", 
          "#771122", "grey","#AA4455", "#DD7788", "#636363")

#generate stacked bar plot
stacked=ggplot(data=pbar, 
               mapping=aes(x=sample,y=abun, fill=Taxa))+
  geom_bar(stat = "identity")+
  facet_grid(~CarbonSource, scales="free_x")+ 
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="Bacterial Order")+
  scale_fill_manual(values=new_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_text(size =10, face="bold"),
        panel.background = element_rect(size=2))+
  guides(fill=guide_legend(ncol=1));

plot(stacked);

##save image 
ggsave(filename="enrich_stackedbar_order.pdf",
       device="pdf",path="./images",
       plot=stacked,
       width=9.3,
       height=7,
       units="in",
       dpi=400)


#trash
#select taxa >1% relative abundance
taxa_rel$AVG=rowMeans(taxa_rel[-1])
taxa_rel=taxa_rel[taxa_rel$AVG>1,]; taxa_rel$AVG=NULL  
rownames(taxa_rel)=taxa_rel$Taxa; taxa_rel$Taxa=NULL  

taxon=rownames(taxa_rel); rtaxon=rep(taxon,ncol(taxa_rel))
taxa_rel<-taxa_rel%>% gather(sample, abun,
                             SP.Fe.Su.a.1:SP.SO4.Su.P.2)  
taxa_rel$taxon=rtaxon

#merge taxa abundance table with metadata
taxa_meta=merge(taxa_rel, meta_data, by="sample") 

stacked=ggplot(data=taxa_meta)+ 
  geom_bar(mapping=aes(x = sample, y = abun, fill = taxon),
           stat="identity")+
  facet_grid(facets=.~CarbonSource,
             scales="free_x") +
  scale_fill_manual(values=new_col)+
  labs(y="Relative Abundance (%)",x="Carbon Source",fill="Bacterial Family")+
  theme_bw()+
  theme(text = element_text(size=12),
        legend.position="right",
        legend.title=element_text(size=12, color="black",face="bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=14, color="black",face="bold"),
        axis.title.x=element_text(size=14, color="black",face="bold"),
        panel.background = element_rect(size=2))

plot(stacked)

##write csv
gin=taxa_rel
gin$AVG=rowMeans(gin[,-1]);
gin=gin[order(-gin$AVG),]
meta_data$newname=paste(meta_data$CarbonSource, meta_data$ElectronAcceptor, sep="_")
colnames(gin)=meta_data$newname[match(colnames(gin),meta_data$sample)]
colnames(gin)[1]=1;
gin=gin[, order(names(gin))]
write.csv(gin, file="enrich_genus.csv",row.names=F)