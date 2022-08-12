
rm(list=ls())
library(ggplot2)

LCM <- read.table("LCM_WGCNA.txt",header = T, check.names = F, sep="\t")
head(LCM)
names(LCM)[3]<-paste("mGM")
names(LCM)[4]<-paste("3C")
names(LCM)[5]<-paste("4BS+V")
names(LCM)[6]<-paste("5BS+V")
names(LCM)[7]<-paste("6BS+V")
names(LCM)[8]<-paste("PM")
names(LCM)[9]<-paste("3PM")
names(LCM)[10]<-paste("4PM")
names(LCM)[11]<-paste("5/6PM")
names(LCM)[12]<-paste("1M")
names(LCM)[13]<-paste("2M")
head(LCM)
#scale data
LCMsc <- scale(LCM[,3:13])
head(LCM[3:13])
head(LCMsc)
# X<-LCM([,1:2])
data <- t(apply(LCM[,3:13], 1, scale))
# Identity name
rownames(data) <- paste(LCM$Gene)
colnames(data) <- letters[1:11]
head(data) 
names(data)[1]<-paste("mGM")
names(data)[2]<-paste("3C")
names(data)[3]<-paste("4BS+V")
names(data)[4]<-paste("5BS+V")
names(data)[5]<-paste("6BS+V")
names(data)[6]<-paste("PM")
names(data)[7]<-paste("3PM")
names(data)[8]<-paste("4PM")
names(data)[9]<-paste("5/6PM")
names(data)[10]<-paste("1M")
names(data)[11]<-paste("2M")
write.csv(data, file = "./LCMz0627.csv")

#loading data
LCMz <- read.table("LCMz.txt",header = T, check.names = F, sep="\t")
head(LCMz)
names(LCMz)[3]<-paste("mGM")
names(LCMz)[4]<-paste("3C")
names(LCMz)[5]<-paste("4BS+V")
names(LCMz)[6]<-paste("5BS+V")
names(LCMz)[7]<-paste("6BS+V")
names(LCMz)[8]<-paste("PM")
names(LCMz)[9]<-paste("3PM")
names(LCMz)[10]<-paste("4PM")
names(LCMz)[11]<-paste("5/6PM")
names(LCMz)[12]<-paste("1M")
names(LCMz)[13]<-paste("2M")
head(LCMz)
#
library(reshape2)
# transfer data 
LCMz = melt(LCMz)
## Using Gene, Group as id variables
head(LCMz)

names(LCMz) <- c("Gene","Group","Stage","Expression")
library(ggplot2)
tiff('WGCNA_0812.tiff', width = 1920, height = 1920, res = 300)
ggplot(LCMz,aes(x=Stage, y=Expression, group=Gene)) + geom_line(color="gray90",size=0.8) + 
  geom_hline(yintercept =0,linetype=2) +
  stat_summary(aes(group=1),fun.y=mean, geom="line", size=1, color="#c51b7d") + 
  facet_wrap(.~Group) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=8,angle = 90, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold"),
        strip.text = element_text(size = 8, face = "bold"))+
  labs(y="Z-score")
dev.off()
