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
## Bray-Curtis & Jaccard (WATER COLUMN data)
## and testing significance of patterns using PERMANOVAs

source(file="scripts/00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load OTU table and sample metadata
#           Chloroplast, Mitochondria, and Unknown were removed
################################################################################
taxa=read.table("data/fgl_otu.txt", sep="\t", header=T);
meta_data=read.table("data/fgl_meta.txt", sep="\t", header=T);
meta_data=meta_data[order(meta_data$sample),];


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
print(adonis(bray.dist~depth,             
             data=meta_data, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~depth,             
             data=meta_data, method = "jaccard",         
             permutations = 999));


################################################################################
#             3. Make PCoA ordination plot 
################################################################################
#obtain coordinates for ordination plot
pcoa_dec=cmdscale(bray.dist, eig=TRUE);  #bray.dist or jac.dist
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "sample"); 
pcoa_met=merge(pcoa,meta_data,by="sample");

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

#set color-palette
box_col=c("#bd0026","#08519c","#31a354","goldenrod2","#8856a7");

#plot the PCoA-color coded by sample depth
pcoa=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=as.factor(depth)),
             size = 4,
             shape=21)+
  scale_fill_manual(values=box_col)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="Depth (m)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold"));

#plot(pcoa);

# #save image 
ggsave(filename="03_WC_pcoa.pdf",
       device="pdf",path="./images",
       plot=pcoa,
       width=6,
       height=4,
       units="in",
       dpi=400);


################################################################################
#             4. Create boxplots showing Bray-Curtis dissimilarity of 
#           microbial community samples at 21m vs. every other depth
################################################################################
#melt distance matrix obtained from vegan
bsim=as.matrix(bray.dist);
bsim=reshape2::melt(bsim);
bsim=bsim[bsim$value!=0,]; #eliminate comparisons of sample against itself

#attach metadata
colnames(bsim)[1]="sample"; 
md=meta_data[,1:2];
bsim=merge(bsim,md, by="sample");
colnames(bsim)=c("sample_1","sample","dissim","depth_1");
bsim=merge(bsim,md, by="sample");
colnames(bsim)=c("sample_2","sample_1", "dissim","depth_1","depth_2");

#remove same sample comparisons (e.g. CR1 vs CR1)
bsim=bsim[!duplicated(bsim$dissim),];

#restrict distance matrix to only comparisons from 21m vs every depth
df1=bsim[(bsim$depth_1==21)|(bsim$depth_2==21),]; 

#add informative labels to comparisons
df1$depths=paste(df1$depth_1, df1$depth_2)
df1$vs = gsub(x=df1$depths,pattern="21", replacement="")
df1$vs=trimws(df1$vs)
df1$category=paste("21 m vs.",df1$vs,"m")
df1$category[df1$vs==""]=" 21 m vs. 21 m"

#create Bray-Curtis boxplots
braybox=ggplot(data=df1)+
  geom_boxplot(mapping=aes(x=as.factor(category), 
                           y=dissim))+
  theme_bw()+ 
  labs(y="Bray-Curtis Dissimilarity (0-1)", 
       x="21m vs")+
  scale_x_discrete(labels=c("21m", "24m","30m","45m","52.5m"))+
  theme(axis.title.y=element_text(size=10, color="black",face="bold"),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        legend.position="none",
        panel.background = element_rect(size=2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank());

plot(braybox);

# #save image 
ggsave(filename="03_WC_braycurtis_plot.pdf",
       device="pdf",path="./images",
       plot=braybox,
       width=5,
       height=4,
       units="in",
       dpi=400);