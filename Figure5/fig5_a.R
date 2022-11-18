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

gender_color = structure(c("#da552f", "#276ce6"), names=c("Female","Male"))

p <- ggplot(data,aes(x=V1, y=V2, color=Gender))+
  geom_point()+
  # geom_text_repel(aes(label=`Row.names`))+
  xlim(x_min, x_max)+
  ylim(y_min, y_max)+
  theme_bw()+
  theme(legend.position = 'none')+
  scale_color_manual(values=gender_color)+
  xlab(xlab)+
  ylab(ylab)
p

use_pro = c("Male","Female")
yplot <- 
  data %>%
  ggplot(.,aes(y=V2,x=Gender, fill=Gender))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        legend.position = "none")+
  geom_signif(comparisons =combn(use_pro,2,list),test='wilcox.test',step_increase = 0.1)+
  #ylim(y_min, y_max)+
  scale_fill_manual(values=gender_color)

xplot <- 
  data %>%
  ggplot(.,aes(x=V1,y=Gender, fill=Gender))+
  geom_boxplot() +
  geom_signif(comparisons =combn(use_pro,2,list),test='wilcox.test',step_increase = 0.1)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.position = "none")+
  #xlim(x_min, x_max)+
  scale_fill_manual(values=gender_color)

grid.arrange(xplot,p,yplot,
             ncol=3,
             layout_matrix = rbind(c(1,1,NA),
                                   c(3,3,4),
                                   c(3,3,4))
)
