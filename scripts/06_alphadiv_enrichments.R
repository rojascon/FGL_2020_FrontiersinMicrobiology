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

## CODE FOR: generating alpha diversity boxplots and 
# testing whether microbial community alpha diversity varies with Carbon source
#or electron acceptor using linear models (ENRICHMENTS data)

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load table of alphadiversity values; 
#                       load sample metadata
################################################################################
div=read.table("data/enrichments_diversity.txt", sep="\t", header=T);
meta_data=read.table("data/enrichments_meta.txt", sep="\t", header=T);


################################################################################
#             2.  Run linear model of alpha-diversity metric ~ 
#                       CarbonSource * Electron Acceptor
#                               Chao1 Richness
################################################################################
#merge alpha-diversity with metadata
am=alpha_met=merge(div, meta_data, by="sample")

am$CarbonSource<-factor(am$CarbonSource, 
                        levels=c("NoC", "methane","acetate",
                                 "propionate","butyrate","lignin",
                                 "chitin","cellulose"))

am$ElectronAcceptor<-factor(am$ElectronAcceptor)

#obtain mean alpha-diversity values for each depth
aggregate(chao1~CarbonSource,am, mean);
aggregate(chao1~CarbonSource,am, sd);

#run linear model and assess significance of predictor variable
mod=lm(chao1~CarbonSource*ElectronAcceptor, data=am)
#summary(mod);
Anova(mod);

##run Tukey HSD posthoc tests
summary(glht(mod, 
             linfct=mcp(CarbonSource="Tukey")),
        test=adjusted("BH"));


################################################################################
#             2.  Run linear model of alpha-diversity metric ~ 
#                       CarbonSource * Electron Acceptor
#                              Shannon Diversity
################################################################################
#obtain mean alpha-diversity values for each depth
aggregate(shannon~CarbonSource,am, mean);
aggregate(shannon~CarbonSource,am, sd);

#run linear model and assess significance of predictor variable
mod=lm(shannon~CarbonSource*ElectronAcceptor, data=am)
#summary(mod);
Anova(mod);

##run Tukey HSD posthoc tests
summary(glht(mod, 
             linfct=mcp(CarbonSource="Tukey")),
        test=adjusted("BH"));

summary(glht(mod, 
             linfct=mcp(ElectronAcceptor="Tukey")),
        test=adjusted("BH"));


################################################################################
#             4.  Generate boxplots of microbial community 
#           alpha diversity ~ Carbon Source * Electron Acceptor
################################################################################
#select color-palette for Electron Acceptor 
box_col=c("#66c2a5", "#e78ac3", "#fdc086")

#build plot!
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
ggsave(filename="06_ENRICH_box_chao1.pdf", 
       device="pdf",path="./images",
       plot=alphabox, 
       width=9.5, 
       height=3, 
       units="in", 
       dpi=400)

