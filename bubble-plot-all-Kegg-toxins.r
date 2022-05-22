library(dplyr)
library(strex)
library(ggplot2)

all<-read.table("ALL_Genes_w_Phylo_Bacteria.txt",sep="\t",header=TRUE)

#10 kegg toxins can be found in the whole data

K16619 <- all[grep("KO:K16619", all$source),]
K11068 <- all[grep("KO:K11068", all$source),]
K11041 <- all[grep("KO:K11041", all$source),]
K11031 <- all[grep("KO:K11031", all$source),]
K11030 <- all[grep("KO:K11030", all$source),]
K10953 <- all[grep("KO:K10953", all$source),]
K01387 <- all[grep("KO:K01387", all$source),]
K01114 <- all[grep("KO:K01114", all$source),]
K01197 <- all[grep("KO:K01197", all$source),]
K03699 <- all[grep("KO:K03699", all$source),]


keggtoxinsall<- rbind(K16619,K11068,K11041,K11031,K11030,K10953,K03699,K01387,K01197,K01114)

write.table(keggtoxinsall,"Kegg-toxins_Genes_w_Phylo_Bacteria.txt",sep="\t",row.names=FALSE)

keggtoxinsall <- read.table("Kegg-toxins_Genes_w_Phylo_Bacteria.txt",sep="\t",header=TRUE)



keggtoxinsall <- keggtoxinsall %>% group_by(Metagenome, Type, Colony, Experiment, TotalGeneCount, class) %>% count(product_name)

keggtoxinsall$Proportion<- (keggtoxinsall$n / keggtoxinsall$TotalGeneCount)*10000

keggtoxinsall$Type<-factor(keggtoxinsall$Type, levels=c("Diseased","Newly Infected","Not Infected","Healthy"))

keggtoxinsall$Metagenome<-factor(keggtoxinsall$Metagenome, levels=c("McAH","McBH","McEH","McGH","McLGH","McD1","McD2","McD3","McD4","McD6","McAI","McBI","McEI","McGNI","McLGI"))

mycolors<-c("Diseased"="red3","Newly Infected"="orange3","Not Infected"="green3","Healthy"="blue3")



zz <- ggplot(keggtoxinsall,aes(Metagenome,product_name))+
  geom_point(aes(size=Proportion,color=Type),alpha=0.8)+ #the alpha value here is the transparency of the color
  scale_size(range = c(1, 17), name="Relative Abundance x 10^-4")+
  scale_y_discrete(limits = rev(levels(keggtoxinsall$product_name)))+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle=90,hjust=1,vjust=0))+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(legend.text = element_text(size=10))+
  scale_colour_manual(values=mycolors)+
  scale_x_discrete(drop = FALSE)
  
  

ggsave("keggtoxins-2.tiff", units="in", width=12, height=8, dpi=300,compression = 'lzw')



## plot this data as box and whiskers instead, showing range of points by disease state and phylogeny

xx <- ggplot(keggtoxinsall,aes(Proportion,class))+
  geom_boxplot()+
  geom_point(aes(color=Type),size=3)+
  theme_bw() +
  labs(x = "Relative Abundance x 10^-4",title = paste("Kegg-toxins")) + 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.y=element_blank())+
  theme(legend.text = element_text(size=12))+
  theme(legend.position="bottom")+
  scale_colour_manual(values=mycolors,name="Tissue Type")


ggsave("keggtoxins-boxplot.tiff", units="in", width=12, height=8, dpi=300,compression = 'lzw')


## plot standar box plot -- in which the type of tissue is used as seprate 

zz <-ggplot(keggtoxinsall, aes(x=Proportion, y=class, color=Type)) +
  geom_boxplot() + theme_bw() + labs(x = "Relative Abundance x 10^-4",title = paste("Kegg-toxins")) + 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  theme(legend.position="bottom")+
  scale_colour_manual(values=mycolors,name="Tissue Type")

  
 

ggsave("keggtoxins-boxplot-3.tiff", units="in", width=12, height=8, dpi=300,compression = 'lzw')

