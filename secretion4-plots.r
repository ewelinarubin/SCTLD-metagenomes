library(dplyr)
library(strex)
library(ggplot2)

all<-read.table("all-keggs-w-taxa.txt",sep="\t",header=TRUE)



K03194 <- all[grep("K03194", all$kegg_id),]
K03195 <- all[grep("K03195", all$kegg_id),]
K03196 <- all[grep("K03196", all$kegg_id),]
K03197 <- all[grep("K03197", all$kegg_id),]
K03198 <- all[grep("K03198", all$kegg_id),]
K03199 <- all[grep("K03199", all$kegg_id),]
K03200 <- all[grep("K03200", all$kegg_id),]
K03201 <- all[grep("K03201", all$kegg_id),]
K03202 <- all[grep("K03202", all$kegg_id),]
K03203 <- all[grep("K03203", all$kegg_id),]
K03204 <- all[grep("K03204", all$kegg_id),]
K03205 <- all[grep("K03205", all$kegg_id),]


sec4 <- rbind(K03194,K03195,K03196,K03197,K03198,K03199,K03200,K03201,K03202,K03203,K03204,K03205)

write.table(sec4,"secretion4-allbacteria.txt",sep="\t",row.names=FALSE)

sec4 <- read.table("secretion4-allbacteria.txt",sep="\t",header=TRUE)



sec4 <- sec4 %>% group_by(msample, type, totalgenecount, class) %>% count(kegg_id)

sec4$proportion <- (sec4$n / sec4$totalgenecount)*10000

sec4$type <- factor(sec4$type, levels=c("diseased","newly infected","not infected","healthy"))

sec4$msample <- factor(sec4$msample, levels=c("McAH","McBH","McEH","McGH","McLGH","McD1","McD2","McD3","McD4","McD6","McAI","McBI","McEI","McGNI","McLGI"))

mycolors <- c("diseased"="red3","newly infected"="orange3","not infected"="green3","healthy"="blue3")



zz <- ggplot(sec4,aes(msample,kegg_id))+
  geom_point(aes(size=proportion,color=type),alpha=0.8)+ #the alpha value here is the transparency of the color
  scale_size(range = c(1, 17), name="Relative Abundance x 10^-4")+
  scale_y_discrete(limits = rev(levels(sec4$kegg_id)))+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle=90,hjust=1,vjust=0))+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(legend.text = element_text(size=10))+
  scale_colour_manual(values=mycolors)+
  scale_x_discrete(drop = FALSE)
  
  

ggsave("secretion-4-bubble-plot.tiff", units="in", width=12, height=8, dpi=300,compression = 'lzw')



## plot this data as box and whiskers instead, showing range of points by disease state and phylogeny

xx <- ggplot(sec4,aes(proportion,class))+
  geom_boxplot()+
  geom_point(aes(color=type),size=3)+
  theme_bw() +
  labs(x = "Relative Abundance x 10^-4",title = paste("Secretion-type-IV")) + 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.y=element_blank())+
  theme(legend.text = element_text(size=12))+
  theme(legend.position="bottom")+
  scale_colour_manual(values=mycolors,name="Tissue Type")


ggsave("secretion4-boxplot.tiff", units="in", width=12, height=8, dpi=300,compression = 'lzw')

xy <- ggplot(sec4,aes(proportion,class,color=type))+
  geom_boxplot()+
  theme_bw() +
  labs(x = "Relative Abundance x 10^-4",title = paste("Secretion-type-IV")) + 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.y=element_blank())+
  theme(legend.text = element_text(size=12))+
  theme(legend.position="bottom")+
  scale_colour_manual(values=mycolors,name="Tissue Type")
  
  
ggsave("secretion4-boxplot-2.tiff", units="in", width=12, height=8, dpi=300,compression = 'lzw')


sec4 <- read.table("secretion4-allbacteria.txt",sep="\t",header=TRUE)
sec4 <- sec4 %>% group_by(msample, type, totalgenecount, family) %>% count(kegg_id)
sec4$proportion <- (sec4$n / sec4$totalgenecount)*10000
sec4$type <- factor(sec4$type, levels=c("diseased","newly infected","not infected","healthy"))
sec4$msample <- factor(sec4$msample, levels=c("McAH","McBH","McEH","McGH","McLGH","McD1","McD2","McD3","McD4","McD6","McAI","McBI","McEI","McGNI","McLGI"))
mycolors <- c("diseased"="red3","newly infected"="orange3","not infected"="green3","healthy"="blue3")


xy <- ggplot(sec4,aes(proportion,family))+
  geom_boxplot()+
  geom_point(aes(color=type),size=3)+
  theme_bw() +
  labs(x = "Relative Abundance x 10^-4",title = paste("Secretion-type-IV")) + 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.y=element_blank())+
  theme(legend.text = element_text(size=12))+
  theme(legend.position="bottom")+
  scale_colour_manual(values=mycolors,name="Tissue Type")


ggsave("secretion5-boxplot.tiff", units="in", width=12, height=8, dpi=300,compression = 'lzw')



yy <- ggplot(sec4,aes(proportion,family,color=type))+
  geom_boxplot()+
  theme_bw() +
  labs(x = "Relative Abundance x 10^-4",title = paste("Secretion-type-IV")) + 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.y=element_blank())+
  theme(legend.text = element_text(size=12))+
  theme(legend.position="bottom")+
  scale_colour_manual(values=mycolors,name="Tissue Type")
  
  
ggsave("secretion6-boxplot-2.tiff", units="in", width=12, height=8, dpi=300,compression = 'lzw')
