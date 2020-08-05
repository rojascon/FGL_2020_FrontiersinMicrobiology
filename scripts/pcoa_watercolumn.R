#Fayettevielle Green Lake Manuscript
#FrontiersinMicrobiology (Rojas et al 2020)
#Script to generate PCoA plots of microbial community composition for WATER COLUMN

source(file="scripts/background.R") #load necessary packages and specifications

#read in your data
taxa=read.table("data/fgl_otu.txt", sep="\t", header=T)
meta_data=read.table("data/fgl_meta.txt", sep="\t", header=T)

########################Part 1: CCA ########################
############################################################

#prepare metadata for constrained correspondence analysis (CCA)
ecopar=meta_data
ecopar=ecopar[order(ecopar$sample),]
rownames(ecopar)=ecopar$sample; ecopar$sample=NULL
ecopar=data.frame(ecopar[, (colnames(ecopar) 
                            %in% c("sulfide_z","DOC", "ammonia", "methane_x"))])

#prepare otu data for constrained correspondence analysis (CCA)
otu_cca=taxa
otu_cca<-otu_cca[, c(6, 7:ncol(otu_cca))]    
colnames(otu_cca)[1]="Taxa"
otu_cca[,-1] <- lapply(otu_cca[,-1], 
                       function(x) (x/sum(x))*100) 
colSums(otu_cca[-1])  

otu_cca$Taxa=NULL 
otu_cca_t=as.data.frame(t(otu_cca))
otu_cca_t=otu_cca_t[order(row.names(otu_cca_t)),]

#run constrained correspondence analysis (CCA)
fgl.cca <- cca(otu_cca_t~ sulfide_z+DOC+ammonia+methane_x, data=ecopar)
print(summary(fgl.cca))
print(anova.cca(fgl.cca, by="terms"))
print(anova.cca(fgl.cca, by="axis"))

#visualize results of constrained correspondence analysis (CCA)
pdf("images/wc_cca_plot.pdf", width =7, height = 6)
plot(fgl.cca, 
     xlim=c(-3,3), ylim=c(-3,3),
     xlab="CCA1 (49%)",
     ylab="CCA2 (10%)",
     display=c("ob","cn","wa")) 

#save the plot
dev.off()

########################Part 2: PCoA by depth ########################
############################################################

#calculate distance matrix for ordination
dist.mat=vegdist(otu_cca_t, method="bray")

#conduct PERMANOVA statistical analyses
print(adonis(dist.mat~depth,             
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
box_col=c("red","blue","green","yellow","purple")
wc_pcoa=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(aes(fill=as.character(depth)), size = 4, shape=21)+
  scale_fill_manual(values=box_col)+
  labs(y="PC2 (20.45%)",x="PC1 (24.74%)",fill="depth (m)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",legend.text=element_text(size=11),
        legend.title=element_text(size=13, face="bold"),
        axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12, face="bold"), 
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12, face="bold"))

plot(wc_pcoa)

# #save image 
ggsave(filename="wc_pcoa_plot.pdf",
       device="pdf",path="./images",
       plot=wc_pcoa,
       width=6,
       height=4,
       units="in",
       dpi=400)

########################Part 3: Bray-Curtis boxplots ########################
######################## 21m vs. every depth ########################
############################################################

#format distance matrix obtained from vegan
bsim=as.matrix(dist.mat)
bsim=reshape2::melt(bsim)
bsim=bsim[bsim$value!=0,] #eliminate comparisons of sample against itself

colnames(bsim)[1]="sample"; 
md=meta_data[,1:2]
bsim=merge(bsim,md, by="sample")
colnames(bsim)=c("sample_1","sample","dissim","depth_1")
bsim=merge(bsim,md, by="sample")
colnames(bsim)=c("sample_2","sample_1", "dissim","depth_1","depth_2")
bsim=bsim[!duplicated(bsim$dissim),]

#restrict distance matrix to only comparisons from 21m vs every depth
df1=bsim[(bsim$depth_1==21)|(bsim$depth_2==21),]; 

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
  labs(y="Bray-Curtis Distance", 
       x="")+
  theme(axis.title.y=element_text(size=10, color="black",face="bold"),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=8),
        legend.position="none",
        panel.background = element_rect(size=2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot(braybox)

# #save image 
ggsave(filename="wc_braycurtis_plot.pdf",
       device="pdf",path="./images",
       plot=braybox,
       width=5,
       height=4,
       units="in",
       dpi=400)