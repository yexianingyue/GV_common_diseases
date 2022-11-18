library(gridExtra) # grid.arrange

sample_map = read.table("sample_group_for_age_gender.tsv", sep="\t", header=T, check.names=F)
load("pcoa.Rdata")

eigs = round(pcoa$eig/sum(pcoa$eig)*100, digits = 2)
point = pcoa$points
data = merge(point, sample_map, by.x='row.names', by.y='ID')

xlab = paste("PCoA 1 ( ", eigs[1], "% )", sep="")
ylab = paste("PCoA 2 ( ", eigs[2], "% )", sep="")



x_min = min(data$V1)
x_max = max(data$V1)
y_min = min(data$V2)
y_max = max(data$V2)

age_color = structure(c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6"),
                         names=c("40~49 years","50~59 years","60~69 years","20~29 years","30~39 years","10~19 years","70~79 years","80~89 years","0~9 years"))

p <- ggplot(data,aes(x=V1, y=V2, color=age1))+
  geom_point()+
  # geom_text_repel(aes(label=`Row.names`))+
  xlim(x_min, x_max)+
  ylim(y_min, y_max)+
  theme_bw()+
  theme(legend.position = 'none')+
  scale_color_manual(values=age_color)+
  ylab(ylab)+
  xlab(xlab)
p

use_pro = unique(data$age1)

sigFunc = function(x){
  if(x < 0.05){ifelse(x<(1e-3), format(x, scientific=T),signif(x,3))} 
  else{NA}}

yplot <- 
  data %>%
  ggplot(.,aes(y=V2,x=age1, fill=age1))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        legend.position = "none")+
  geom_signif(comparisons =combn(use_pro,2,list),test='wilcox.test',step_increase = 0.1,map_signif_level=sigFunc)+
  # ylim(y_min, y_max)+
  scale_fill_manual(values=age_color)+
  scale_y_continuous(breaks=pretty(data$V2))
yplot

xplot <- 
  data %>%
  ggplot(.,aes(x=V1,y=age1, fill=age1))+
  geom_boxplot() +
  geom_signif(comparisons =combn(use_pro,2,list),test='wilcox.test',step_increase = 0.1,map_signif_level=sigFunc)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.position = "none")+
  # xlim(x_min, x_max)
  scale_fill_manual(values=age_color)+
  scale_x_continuous(breaks=pretty(data$V1))
xplot
# pcoa_age_with_boxplot.pdf
grid.arrange(xplot,p,yplot,
             ncol=3,
             layout_matrix = rbind(c(1,1,NA),
                                   c(3,3,4),
                                   c(3,3,4))
)
