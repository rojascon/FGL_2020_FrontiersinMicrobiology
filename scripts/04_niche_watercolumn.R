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

## CODE FOR: generating niche breadth ~ abundance for each bacterial taxon
#niche breadth was calculated using Levin's (1968) niche breadth index 
#(WATER COLUMN data)

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Read in your table of Levin's niche breadth values
################################################################################
niche=read.table("data/fgl_niche.txt", sep="\t", header=T);

#calculate +/- 1 SD from the mean niche breadth
MeanNiche=mean(niche$niche_breadth);
SDNiche=sd(niche$niche_breadth);
upper=MeanNiche+SDNiche;
lower=MeanNiche-SDNiche;


################################################################################
#             2.  Scatter of Levin's niche breadth vs. taxa abunance
################################################################################
#set up color palette
new_col=c("#1b9e77", "#e6ab02", "darkorchid2", "#e7298a",
          "steelblue1","#d9d9d9");

#plot! specialist and generalist bacterial taxa of interest were shown in color
#all other taxa were in gray
#dashed lines denote boundaries for generalists or specialists

niche_wc=ggplot(data=niche)+ 
  geom_point(mapping=aes(x = MeanRelAbund, 
                         y = niche_breadth, 
                         colour = colorcode),
             size=2)+
  scale_colour_manual(values=new_col)+
  scale_y_continuous(breaks=seq(from=1, to=5, by=1))+
  labs(y="Levin's niche breadth",
       x="Average Relative Abundance (log)",
       colour="")+
  geom_hline(yintercept=upper, linetype='dashed', col = 'black')+
  geom_hline(yintercept=lower, linetype='dashed', col = 'black')+
  theme_bw()+  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.position="right",
        legend.text=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12, face="bold"), 
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12, face="bold"));

#plot(niche_wc);

##save image 
ggsave(filename="04_WC_niche.pdf",
       device="pdf",path="./images",
       plot=niche_wc,
       width=6,
       height=4,
       units="in",
       dpi=400);