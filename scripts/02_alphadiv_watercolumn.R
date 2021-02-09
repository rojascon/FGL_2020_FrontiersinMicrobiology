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
# testing whether microbial community alpha diversity varies with sample depth
# using linear models (WATER COLUMN data)


################################################################################
#             1.  Load table of alphadiversity values; 
#                       load sample metadata
################################################################################
div=read.table("data/fgl_diversity.txt", sep="\t", header=T);
meta_data=read.table("data/fgl_meta.txt", sep="\t", header=T);


################################################################################
#             3.  Run linear model of alpha-diversity metric ~ sample depth
#                               Chao1 Richness
################################################################################
#merge alpha-diversity with metadata
am=alpha_met=merge(div, meta_data, by="sample");

#obtain mean alpha-diversity values for each depth
aggregate(chao1~depth,am, mean);
aggregate(chao1~depth,am, sd);

#run linear model and assess significance of predictor variable
am$depth=as.factor(am$depth);
mod=lm(chao1~depth, data=am);
#summary(mod);
Anova(mod);


################################################################################
#             3.  Run linear model of alpha-diversity metric ~ sample depth
#                                 Shannon Diversity
################################################################################
#merge alpha-diversity with metadata
am=alpha_met=merge(div, meta_data, by="sample");
View(am);

#obtain mean alpha-diversity values for each depth
aggregate(shannon~depth,am, mean);
aggregate(shannon~depth,am, sd);

#run linear model and assess significance of predictor variable
am$depth=as.factor(am$depth);
mod=lm(shannon~depth, data=am);
#summary(mod);
Anova(mod);

##run Tukey HSD posthoc tests
summary(glht(mod, 
             linfct=mcp(depth="Tukey")),
        test=adjusted("BH"));


################################################################################
#             4.  Generate boxplots of microbial community 
#                       alpha diversity ~ sample depth
################################################################################
#select color-palette
box_col=c("#8856a7", "goldenrod2","#31a354","#08519c","#bd0026" )
am$depth=forcats::fct_rev(factor(am$depth));

#build plot!
alphabox=ggplot(data=am)+
  geom_boxplot(mapping=aes(x=as.factor(depth), 
                           y=shannon, 
                           fill=as.factor(depth)))+
  coord_flip()+
  theme_bw()+ 
  labs(y="Shannon Diversity", 
        x="Depth (m)")+
  scale_fill_manual(values=box_col)+
  theme(axis.title.x=element_text(size=14, color="black",face="bold"),
        axis.title.y=element_text(size=14, color="black",face="bold"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        legend.position="none",
        panel.background = element_rect(size=2));

plot(alphabox);

##save image 
ggsave(filename="02_WC_box_shannon.pdf",
       device="pdf",path="./images",
       plot=alphabox,
       width=3.5,
       height=3,
       units="in",
       dpi=400);


