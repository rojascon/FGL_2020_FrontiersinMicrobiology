#Fayettevielle Green Lake Manuscript
#FrontiersinMicrobiology (Rojas et al 2020)
#Script to generate stacked bar plots of microbial community composition for ENRICHMENTS

source(file="scripts/background.R") #load necessary packages and specifications

#read in your data
taxa=read.table("data/enrichments_otu.txt", sep="\t", header=T)
meta_data=read.table("data/enrichments_meta.txt", sep="\t", header=T)

meta_data$CarbonSource<-factor(meta_data$CarbonSource, 
                        levels=c("NoC", "methane","acetate",
                                 "propionate","butyrate","lignin",
                                 "chitin","cellulose"))

#select your taxonomic level
#1-Kingdom 2-Phylum 3-Class 4-Order 5-Family 6-Genus
taxa<-taxa[, c(5, 7:ncol(taxa))] 
print(colnames(taxa[1])); colnames(taxa)[1]="Taxa"

#calculate proportions of bacterial taxa
taxa_rel=aggregate(.~Taxa, taxa, sum)
taxa_rel[,-1] <- lapply(taxa_rel[,-1], function(x) (x/sum(x))*100) 
print(colSums(taxa_rel[-1]))   

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

#generate stacked barplot
new_col=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
          "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", 
          "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", 
          "#771122", "#AA4455", "#DD7788","grey", "black")

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

##save image 
ggsave(filename="enrich_stackedbar.pdf",
       device="pdf",path="./images",
       plot=stacked,
       width=13,
       height=7,
       units="in",
       dpi=400)