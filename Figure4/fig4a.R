library(dplyr)
library(tidyselect)
library(ggplot2)
library(reshape2)

dt = read.table("../00.data/vOTU.reads_norm.profile.relative.family",sep="\t", header=T, row.names=1,check.names=F)
sample_map = read.table("../00.data/sample.info.20220103.tsv",sep="\t", header=T,check.names=F)
#proportion = read.table("../00.data/phylum_proportion.tsv",sep="\t")

colors = c(
  "#a6cee3","#1f78b4","#b2df8a","#33a02c",
  "#d9d9d9","#bc80bd","#ccebc5","#ffed6f",
  "#cab2d6","#6a3d9a","#b15928",
  "#fb9a99","#e31a1c","#fdbf6f","#ff7f00",
  "#8dd3c7","#fb8072","#bebada","#fdb462","#b3de69","#fccde5"
)

sample_map <- sample_map %>% filter(Project_2 == "Discover") %>% dplyr::select(Sample, Project_1) %>% unique()
#y <- proportion %>% filter(value > 10) %>% rownames()
#cutoff = 1e-5
#x = names(rowMeans(dt)>cutoff)
#select_names = intersect(x, y)
dt1 = dt[,sample_map$Sample]

#---------
# sort
dtf = dt1[order(rowSums(dt1),decreasing = T),]
rowMeans(dtf)[1:18]
top=18
if(nrow(dtf) > top){
  dt_f = dtf[1:(top-1),]
  dtff = rbind(dt_f, 1-colSums(dt_f))
  rownames(dtff)[top] = "other"
  dtf = dtff
  rm(dtff)
  rm(dt_f)
}


#dtf = dt1[order(rowSums(dt1),decreasing = T),]
dtf = dtf[,order(dtf[1,], decreasing = T)]

dl = melt(as.matrix(dtf))
dm = merge(sample_map, dl, by.x='Sample', by.y='Var2')
data <- dm %>%
  select(Var1,Sample,value, Project_1)

data$Var1 = factor(data$Var1, levels=unique(rownames(dtf)))
data$Sample = factor(data$Sample, levels=rev(unique(colnames(dtf))))
data1 = unique(data)
p1 <- ggplot(data1) +
  geom_bar(aes(x=Sample, y=value, fill=Var1), stat='identity', position='stack')+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  scale_fill_manual(values=colors)+
  scale_y_continuous(expand = c(0,0))
p1
#---------
# pie plot
project_order = read.table("../00.data/project_order.list",sep="\t", header=T)
data2 = data1 %>%
  group_by(Project_1,Var1) %>%
  summarise(value=mean(value))
data2$Project_1 = factor(data2$Project_1, levels=unique(project_order$project_1))
p2 <- ggplot(data2, aes(x="",y=value, fill=Var1))+
  geom_bar(stat="identity",width=1, color='white')+
  coord_polar(theta="y",direction=1)+
  facet_wrap(Project_1~.,nrow = 4)+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=colors)
  #scale_fill_brewer(palette="Set3")
# scale_fill_futurama() # 好像颜色很多
p2

library(ggpubr)
ggarrange(plotlist=list(p1,p2), ncol=1,common.legend=T, widths=c(1,2), heights=c(1,1))
