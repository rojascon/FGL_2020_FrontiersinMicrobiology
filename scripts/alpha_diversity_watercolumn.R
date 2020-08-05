#Fayettevielle Green Lake Manuscript
#FrontiersinMicrobiology (Rojas et al 2020)
#Script to generate alpha diversity boxplots and kruskal-wallis test for WATER COLUMN

source(file="scripts/background.R") #load necessary packages and specifications

#read in your data
div=read.table("data/fgl_diversity.txt", sep="\t", header=T)
colnames(div)[1]= "sample"
meta_data=read.table("data/fgl_meta.txt", sep="\t", header=T)

#merge alpha-diversity with metadata
am=alpha_met=merge(div, meta_data, by="sample")
View(am)

#run kruskal wallis test and post-hoc test
print(
  kruskal.test(am$chao1~as.factor(am$depth))
  )
print(
  dunnTest(am$chao1~as.factor(am$depth),data=am,method="bh")
)

#generate alpha-diversity boxplots
box_col=c("#bd0026","#08519c","#31a354","goldenrod2","#8856a7")

alphabox=ggplot(data=am)+
  geom_boxplot(mapping=aes(x=as.factor(depth), 
                           y=chao1, 
                           fill=as.factor(depth)))+
  theme_bw()+ 
  labs(y="Chao 1 Richness", 
        x="Depth (m)")+
  scale_fill_manual(values=box_col)+
  theme(axis.title.x=element_text(size=14, color="black",face="bold"),
        axis.title.y=element_text(size=14, color="black",face="bold"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        legend.position="none",
        panel.background = element_rect(size=2))

plot(alphabox)

# #save image 
ggsave(filename="wc_boxplot.pdf",
       device="pdf",path="./images",
       plot=alphabox,
       width=3.5,
       height=3,
       units="in",
       dpi=400)
