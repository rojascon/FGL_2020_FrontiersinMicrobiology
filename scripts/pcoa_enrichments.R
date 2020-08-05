#Fayettevielle Green Lake Manuscript
#FrontiersinMicrobiology (Rojas et al 2020)
#Script to generate PCoA plots and Bray-Curtis plots 
#for ENRICHMENT data

source(file="scripts/background.R") #load necessary packages and specifications

#read in your data
taxa=read.table("data/enrichments_otu.txt", sep="\t", header=T)
meta_data=read.table("data/enrichments_meta.txt", sep="\t", header=T)

meta_data$CarbonSource<-factor(meta_data$CarbonSource, 
                               levels=c("NoC", "methane","acetate",
                                        "propionate","butyrate","lignin",
                                        "chitin","cellulose"))

######################Part 1: PCoA by Carbon source and Electron Acceptor ########
############################################################

#prepare otu data for PCoA analyses
otu_pcoa=taxa
otu_pcoa<-otu_pcoa[, c(6, 7:ncol(otu_pcoa))]    
colnames(otu_pcoa)[1]="Taxa"
otu_pcoa[,-1] <- lapply(otu_pcoa[,-1], 
                       function(x) (x/sum(x))*100) 
colSums(otu_pcoa[-1])  

otu_pcoa$Taxa=NULL 
otu_pcoa_t=as.data.frame(t(otu_pcoa))
otu_pcoa_t=otu_pcoa_t[order(row.names(otu_pcoa_t)),]

#calculate distance matrix for ordination
dist.mat=vegdist(otu_pcoa_t, method="bray")

#conduct PERMANOVA statistical analyses
print(adonis(dist.mat~CarbonSource*ElectronAcceptor,             
             data=meta_data, method = "bray",         
             permutations = 999)) 

#obtain coordinates for ordination plot
pcoa_dec=cmdscale(dist.mat, eig=TRUE)
pcoa=as.data.frame(pcoa_dec$points)
colnames(pcoa)=c("Axis1","Axis2") 
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "sample"); 
pcoa_met=merge(pcoa,meta_data,by="sample") 

#calculate axis % explained
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; sum(pcoa_per) 
pcoa_per[1]; pcoa_per[2] #for axis 1 and 2 respectively

#visualize PCoA plot
box_col=c("#66c2a5","#fc8d62", "#8da0cb", "#e78ac3","#a6d854","#ffd92f","#e5c494","#b3b3b3") 

enrich_pcoa=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(colour=CarbonSource,
                         shape=ElectronAcceptor), size = 5)+
  scale_colour_manual(values=box_col)+
  labs(y="PC2 (10.93%)",x="PC1 (24.23%)",
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
        axis.title.y=element_text(size=12, face="bold"))

plot(enrich_pcoa)

# #save image 
ggsave(filename="enrich_pcoa_plot.pdf",
       device="pdf",path="./images",
       plot=enrich_pcoa,
       width=7,
       height=4.5,
       units="in",
       dpi=400)

########################Part 3: Bray-Curtis boxplots ########################
######################## Within vs. Between Carbon sources ########################
############################################################

#format distance matrix obtained from vegan
bsim=as.matrix(dist.mat)
bsim=reshape2::melt(bsim)
bsim$sim=1-bsim$value
bsim=bsim[bsim$sim!=1,] #eliminate comparisons of sample against itself

bsim$carbon_1=gsub(x=bsim$Var1,  #extract carbon source name
                   pattern="SP|Su|Fe|S|SO4|[0123456789.]", replacement="")
bsim$carbon_2=gsub(x=bsim$Var2,
                   pattern="SP|Su|Fe|S|SO4|[0123456789.]", replacement="")

bsim$electron_1=gsub(x=bsim$Var1, #extract electron source name
                     pattern="SP.|Su.|[^FeSO4]", replacement="")

bsim$electron_2=gsub(x=bsim$Var2,
                     pattern="SP.|Su.|[^FeSO4]", replacement="")

colnames(bsim)=c("sample_1","sample_2","value","sim",
                 "carbon_1","carbon_2","electron_1","electron_2")

bsim=bsim[!duplicated(bsim$sim),]

#restrict data frame to only within carbon source comparisons
#will compare samples of same carbon source and same e- acceptor ="within" category
#will also compare samples of same carbon source and different e- acceptor ="between category"

df=bsim[bsim$carbon_1==bsim$carbon_2,]

same_electron=df[df$electron_1==df$electron_2,]
same_electron$category="W"  #W= within

diff_electron=df[df$electron_1!=df$electron_2,]
diff_electron$category="B"  #B= between

df2=rbind(same_electron, diff_electron)
df2$category<-factor(df2$category, levels=c("W","B"))

df2$carbon_1<-factor(df2$carbon_1, levels=c("NC","m","a","P","B","Li","Ch","C"))

carbon.labs=c("NoC", "methane","acetate",
                "propionate","butyrate","lignin",
                "chitin","cellulose")
names(carbon.labs) <- levels(df2$carbon_1)
              
#create Bray-Curtis boxplots
bray_enrich=ggplot(data=df2)+
  geom_boxplot(mapping=aes(x=as.factor(category), 
                           y=sim))+
  facet_grid(.~carbon_1, labeller = labeller(carbon_1=carbon.labs))+
  theme_bw()+ 
  labs(y="Bray-Curtis Distance", 
       x="")+
  theme(axis.title.y=element_text(size=12, color="black",face="bold"),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        legend.position="none",
        panel.background = element_rect(size=2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=11, face="bold"))

plot(bray_enrich)

# #save image 
ggsave(filename="enrich_braycurtis_plot.pdf",
       device="pdf",path="./images",
       plot=bray_enrich,
       width=8,
       height=3,
       units="in",
       dpi=400)

