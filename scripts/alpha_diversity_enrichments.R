#Fayettevielle Green Lake Manuscript
#FrontiersinMicrobiology (Rojas et al 2020)
#Script to generate alpha diversity boxplots and kruskal-wallis test for ENRICHMENTS

source(file="scripts/background.R") #load necessary packages and specifications

#read in your data
div=read.table("data/enrichments_diversity_2021.txt", sep="\t", header=T)
#colnames(div)[1]= "sample"
meta_data=read.table("data/enrichments_meta.txt", sep="\t", header=T)

#merge alpha-diversity with metadata
am=alpha_met=merge(div, meta_data, by="sample")
#am=final[complete.cases(am), ]
am$CarbonSource<-factor(am$CarbonSource, 
                        levels=c("NoC", "methane","acetate",
                                 "propionate","butyrate","lignin",
                                 "chitin","cellulose"))
View(am)

#obtain mean diversity values for each carbon source
aggregate(shannon~CarbonSource,am, mean )

#run kruskal wallis test and post-hoc test
print(
  kruskal.test(am$shannon~as.factor(am$CarbonSource))
)
print(
  dunnTest(am$shannon~as.factor(am$CarbonSource),data=am,method="bh")
)

mod=lm(shannon~CarbonSource*ElectronAcceptor, data=am)

summary(glht(mod, 
             linfct=mcp(ElectronAcceptor="Tukey")),
        test=adjusted("BH"))

#generate alpha-diversity boxplots
box_col=c("#66c2a5", "#e78ac3", "#fdc086")

alphabox=ggplot(data=am)+
  geom_boxplot(mapping=aes(x=ElectronAcceptor, 
                           y=chao1, 
                           fill=ElectronAcceptor))+
  facet_wrap(~CarbonSource, scales="free_x",nrow = 1)+
  theme_bw()+ 
  labs(y="Chao 1 Richness", 
       x=" ")+
  scale_fill_manual(values=box_col)+
  theme(axis.title.x=element_text(size=11, color="black",face="bold"),
        axis.title.y=element_text(size=11, color="black",face="bold"),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=11),
        legend.position="right",
        legend.title=element_text(size=11, face="bold"),
        legend.text=element_text(size=11),
        panel.background = element_rect(size=1.3, color="black"),
        strip.text = element_text(size=10, face="bold"))

plot(alphabox)

##save image 
ggsave(filename="enrich_boxplot.pdf", 
       device="pdf",path="./images",
       plot=alphabox, 
       width=8.5, 
       height=3, 
       units="in", 
       dpi=400)

