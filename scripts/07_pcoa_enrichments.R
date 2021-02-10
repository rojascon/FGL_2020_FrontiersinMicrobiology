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

## CODE FOR: generating PCoA plots based on microbial commmunity dissimilarity matrices
## Bray-Curtis & Jaccard (ENRICHMENTS data)
## and testing significance of patterns using PERMANOVAs
## and creating boxplots of average Bray-Curtis dissimilarity

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load OTU table and sample metadata
#           Chloroplast, Mitochondria, and Unknown were removed
################################################################################
taxa=read.table("data/enrichments_otu.txt", sep="\t", header=T);
meta_data=read.table("data/enrichments_meta.txt", sep="\t", header=T);
meta_data=meta_data[order(meta_data$sample),];

meta_data$CarbonSource<-factor(meta_data$CarbonSource, 
                               levels=c("NoC", "methane","acetate",
                                        "propionate","butyrate","lignin",
                                        "chitin","cellulose"));


################################################################################
#             2. Generate the 2 types of distance matrices
################################################################################
#transpose OTU table so OTUs are columns
taxadf=taxa[,7:ncol(taxa)];
taxadf=t(taxadf);
taxadf=taxadf[order(row.names(taxadf)),];

###BRAY-CURTIS distance
bray<-apply(taxadf, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
#print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

###JACCARD distance
jac=(taxadf>0)*1;
#print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");


################################################################################
#             3. Conduct PERMANOVAs
################################################################################
#does microbial community structure vary with sample depth?
#Bray-Curtis
print(adonis(bray.dist~CarbonSource*ElectronAcceptor,             
             data=meta_data, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~CarbonSource*ElectronAcceptor,             
             data=meta_data, method = "jaccard",         
             permutations = 999));


################################################################################
#             3. Make PCoA ordination plot 
################################################################################
pcoa_dec=cmdscale(bray.dist, eig=TRUE);  #bray.dist or jac.dist
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "sample"); 
pcoa_met=merge(pcoa,meta_data,by="sample");

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

#set color-palette for Carbon source
box_col=c("#66c2a5","#fc8d62", "#8da0cb", "#e78ac3","#a6d854","#ffd92f",
          "#e5c494","#b3b3b3") 

#plot the PCoA-color coded by Carbon Source and 
#with shape indicating Electron Acceptor
enrich_pcoa=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(colour=CarbonSource,
                         shape=ElectronAcceptor), 
             size = 5)+
  scale_colour_manual(values=box_col)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       colour="Carbon Source", 
       shape="Electron Acceptor")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=11),
        legend.title=element_text(size=13, face="bold"),
        axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12, face="bold"), 
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12, face="bold"));

plot(enrich_pcoa)

##save image 
ggsave(filename="07_ENRICH_pcoa.pdf",
       device="pdf",path="./images",
       plot=enrich_pcoa,
       width=7,
       height=4.5,
       units="in",
       dpi=400)


################################################################################
#             4. Create boxplots showing Bray-Curtis similarity for 
#           microbial community samples WITHIN VS BETWEEN Electron Acceptors
#                         for each Carbon Sourcee
################################################################################
#melt distance matrix obtained from vegan
bsim=as.matrix(bray.dist);
bsim=reshape2::melt(bsim);
bsim$sim=1-bsim$value;
bsim=bsim[bsim$sim!=1,]; #eliminate comparisons of sample against itself

#attach metadata
bsim$carbon_1=gsub(x=bsim$Var1,  #extract carbon source name
                   pattern="SP|Su|Fe|S|SO4|[0123456789.]", replacement="");
bsim$carbon_2=gsub(x=bsim$Var2,
                   pattern="SP|Su|Fe|S|SO4|[0123456789.]", replacement="");

bsim$electron_1=gsub(x=bsim$Var1, #extract electron source name
                     pattern="SP.|Su.|[^FeSO4]", replacement="");

bsim$electron_2=gsub(x=bsim$Var2,
                     pattern="SP.|Su.|[^FeSO4]", replacement="");

colnames(bsim)=c("sample_1","sample_2","value","sim",
                 "carbon_1","carbon_2","electron_1","electron_2");

#remove same sample comparisons 
#(e.g. SP.Fe.Su.a.2 vs SP.Fe.Su.a.2)
bsim=bsim[!duplicated(bsim$sim),];

#restrict data frame to WITHIN OR BETWEEN categories
#"within"= same carbon source, same e- acceptor
#"between"= same carbon source, different e- acceptor 

df=bsim[bsim$carbon_1==bsim$carbon_2,];

same_electron=df[df$electron_1==df$electron_2,];
same_electron$category="W"; #W= within

diff_electron=df[df$electron_1!=df$electron_2,];
diff_electron$category="B";  #B= between

df2=rbind(same_electron, diff_electron);
df2$category<-factor(df2$category, levels=c("W","B"));

#make carbon source labels informative
df2$carbon_1<-factor(df2$carbon_1, 
                     levels=c("NC","m","a","P","B","Li","Ch","C"));

carbon.labs=c("NoC", "methane","acetate",
                "propionate","butyrate","lignin",
                "chitin","cellulose");
names(carbon.labs) <- levels(df2$carbon_1);
              
#create Bray-Curtis boxplots
bray_enrich=ggplot(data=df2)+
  geom_boxplot(mapping=aes(x=as.factor(category), 
                           y=sim))+
  facet_grid(.~carbon_1, 
             labeller = labeller(carbon_1=carbon.labs))+
  theme_bw()+ 
  labs(y="Bray-Curtis Similarity (0-1)", 
       x="")+
  theme(axis.title.y=element_text(size=12, color="black",face="bold"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        legend.position="none",
        panel.background = element_rect(size=2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=10, face="bold"));

#plot(bray_enrich);

# #save image 
ggsave(filename="07_ENRICH_braycurtis_plot.pdf",
       device="pdf",path="./images",
       plot=bray_enrich,
       width=8,
       height=3,
       units="in",
       dpi=400)

