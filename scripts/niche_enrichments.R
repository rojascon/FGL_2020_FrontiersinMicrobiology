#Fayettevielle Green Lake Manuscript
#FrontiersinMicrobiology (Rojas et al 2020)
#Script to generate niche breadth plots of microbial community composition for ENRICHMENTS
#niche breadth was calculated using Levin's (1968) niche breadth index

source(file="scripts/background.R") #load necessary packages and specifications

#read in your data
niche=read.table("data/enrichments_niche.txt", sep="\t", header=T)

#set up plot color palette
new_col=c("#7fc97f","#beaed4","#fdc086")

#plot ASV niche breadth vs. ASV average relative abundance across enrichment samples
#each sample was color-coded by niche breadth category:
#0-specialist #1-other (i.e. average niche breadth) 2-generalist

niche_enrich=ggplot(data=niche)+ 
  geom_point(mapping=aes(x = MeanRelAbund, y = niche_breadth, colour = category))+
  scale_colour_manual(values=new_col)+
  scale_y_continuous(breaks=seq(from=1, to=18, by=3))+
  scale_x_continuous(breaks=seq(from=-1, to=-6, by=-1))+
  labs(y="Levin's niche breadth",
       x="ASV Average Relative Abundance (log)",
       colour="")+
  theme_bw()+  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.position="right",
        legend.text=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12, face="bold"), 
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12, face="bold"))

plot(niche_enrich)

##save image 
ggsave(filename="enrich_niche.pdf",
       device="pdf",path="./images",
       plot=niche_enrich,
       width=6,
       height=4,
       units="in",
       dpi=400)

