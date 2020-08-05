#Fayettevielle Green Lake Manuscript
#FrontiersinMicrobiology (Rojas et al 2020)
#Script to generate stacked bar plots of microbial community composition for WATER COLUMN

source(file="scripts/background.R") #load necessary packages and specifications

#read in your data
taxa=read.table("data/fgl_otu.txt", sep="\t", header=T)
meta_data=read.table("data/fgl_meta.txt", sep="\t", header=T)

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
taxa_rel<-taxa_rel%>% gather(sample, abun, CR2:CR7)  
taxa_rel$taxon=rtaxon

#merge taxa abundance table with metadata
taxa_meta=merge(taxa_rel, meta_data, by="sample") 

#generate stacked bar plot
new_col=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", 
          "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", 
          "#F6C141", "#F1932D", "#E8601C", "#DC050C")

stacked=ggplot(data=taxa_meta)+ 
  geom_bar(mapping=aes(x = sample, y = abun, fill = taxon),
           stat="identity")+
  facet_grid(facets=.~depth,
             scales="free_x") +
  scale_fill_manual(values=new_col)+
  labs(y="Relative Abundance (%)",x="Depth (m)",fill="Bacterial Family")+
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
ggsave(filename="wc_stackedbar.pdf",
       device="pdf",path="./images",
       plot=stacked,
       width=10,
       height=5,
       units="in",
       dpi=400)





