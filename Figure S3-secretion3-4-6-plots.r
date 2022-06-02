library(dplyr)
library(strex)
library(ggplot2)

all<-read.table("all-keggs-w-taxa.txt",sep="\t",header=TRUE)

#KEGGs for secretion 3 
#secretion-III-keggs K03219, K03221, K03222, K03223, K03224, K03225, K03226, K03227, K03228, K03229, K03230, K04056 K04057, K04058, K04059


K03219 <- all[grep("K03219", all$kegg_id),]
K03221 <- all[grep("K03221", all$kegg_id),]
K03222 <- all[grep("K03222", all$kegg_id),]
K03223 <- all[grep("K03223", all$kegg_id),]
K03224 <- all[grep("K03224", all$kegg_id),]
K03225 <- all[grep("K03225", all$kegg_id),]
K03226 <- all[grep("K03226", all$kegg_id),]

K03227 <- all[grep("K03227", all$kegg_id),]
K03228 <- all[grep("K03228", all$kegg_id),]
K03229 <- all[grep("K03229", all$kegg_id),]
K03230 <- all[grep("K03230", all$kegg_id),]
K04056 <- all[grep("K04056", all$kegg_id),]
K04057 <- all[grep("K04057", all$kegg_id),]
K04058 <- all[grep("K04058", all$kegg_id),]
K04059 <- all[grep("K04059", all$kegg_id),]


sec3<- rbind(K03219,K03221,K03222,K03223,K03224,K03225,K03226,K03227,K03228,K03229,K03230,K04056,K04057,K04058,K04059)

write.table(sec3,"secretion3-allbacteria.txt",sep="\t",row.names=FALSE)

sec3 <- read.table("secretion3-allbacteria.txt",sep="\t",header=TRUE)

sec3 <- sec3 %>% group_by(msample, type, totalgenecount, class) %>% count(kegg_id)

sec3$proportion <- (sec3$n / sec3$totalgenecount)*10000

sec3$type <- factor(sec3$type, levels=c("diseased","newly infected","not infected","healthy"))

sec3$msample <- factor(sec3$msample, levels=c("McAH","McBH","McEH","McGH","McLGH","McD1","McD2","McD3","McD4","McD6","McAI","McBI","McEI","McGNI","McLGI"))

mycolors <- c("diseased"="red3","newly infected"="orange3","not infected"="green3","healthy"="blue3")

xx <- ggplot(sec3,aes(proportion,class))+
  geom_boxplot()+
  geom_point(aes(color=type),size=3)+
  theme_bw() +
  labs(x = "Relative Abundance x 10^-4",title = paste("Secretion-type-III")) + 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.y=element_blank())+
  theme(legend.text = element_text(size=12))+
  theme(legend.position="bottom")+
  scale_colour_manual(values=mycolors,name="Tissue Type")


ggsave("secretion3-boxplot.tiff", units="in", width=6, height=4, dpi=300,compression = 'lzw')



#keggs for secretion type IV 

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

xy <- ggplot(sec4,aes(proportion,class))+
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

#keggs for secretion type VI 

K11891 <- all[grep("K11891", all$kegg_id),]
K11892 <- all[grep("K11892", all$kegg_id),]
K11903 <- all[grep("K11903", all$kegg_id),]
K11904 <- all[grep("K11904", all$kegg_id),]
K11906 <- all[grep("K11906", all$kegg_id),]
K11907 <- all[grep("K11907", all$kegg_id),]
K11912 <- all[grep("K11912", all$kegg_id),]
K11913 <- all[grep("K11913", all$kegg_id),]
K11915 <- all[grep("K11915", all$kegg_id),]


sec6 <- rbind(K11891,K11892,K11903,K11904,K11906,K11907,K11912,K11913,K11915)

write.table(sec6,"secretion6-allbacteria.txt",sep="\t",row.names=FALSE)

sec6 <- read.table("secretion6-allbacteria.txt",sep="\t",header=TRUE)

sec6 <- sec6 %>% group_by(msample, type, totalgenecount, class) %>% count(kegg_id)

sec6$proportion <- (sec6$n / sec6$totalgenecount)*10000

sec6$type <- factor(sec6$type, levels=c("diseased","newly infected","not infected","healthy"))

sec6$msample <- factor(sec6$msample, levels=c("McAH","McBH","McEH","McGH","McLGH","McD1","McD2","McD3","McD4","McD6","McAI","McBI","McEI","McGNI","McLGI"))

mycolors <- c("diseased"="red3","newly infected"="orange3","not infected"="green3","healthy"="blue3")

xz <- ggplot(sec6,aes(proportion,class))+
  geom_boxplot()+
  geom_point(aes(color=type),size=3)+
  theme_bw() +
  labs(x = "Relative Abundance x 10^-4",title = paste("Secretion-type-VI")) + 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.y=element_blank())+
  theme(legend.text = element_text(size=12))+
  theme(legend.position="bottom")+
  scale_colour_manual(values=mycolors,name="Tissue Type")


ggsave("secretion6-boxplot.tiff", units="in", width=12, height=8, dpi=300,compression = 'lzw')











