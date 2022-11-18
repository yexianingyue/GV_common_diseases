library(dplyr)
library(tidyselect)
library(reshape2)
library(ggplot2)
library(vegan)
library(ggsignif)

# 生成画图数据，后面不需要再运行
# dt = read.table("../00.data/case_control.profile.tpm", sep="\t", header=T, row.names=1, check.names=F)
# #dt[dt<1e-5] = 0
# 
# diver = data.frame(alpha=diversity(t(dt), index='shannon'))
# 
# 
# 
# sample_map = read.table("../00.data/sample.info.20220323.tsv",sep="\t", header=T)
# project_map = read.table("../00.data/project.group",sep="\t", header=T)
# sample_map <- sample_map %>%
#   filter(Project_2 == "Discover") %>%
#   dplyr::select(Sample, Group,Project) %>%
#   unique()
# 
# dm = merge(diver, sample_map, by.x='row.names', by.y='Sample')
# 
# ggplot(dm, aes(x=Group, fill=Group, y=alpha))+
#   geom_boxplot()+
#   geom_signif(comparisons=list(c("Disease", "Control")),
#               test='wilcox.test')+
#   facet_wrap(Project~., scale='free')
# 
# 
# data <- dm %>%
#   group_by(Project, Group) %>%
#   summarise(value = mean(alpha)) %>%
#   dcast(Project~Group) %>%
#   # 相对于健康人，病人的多样性减少还是增加
#   mutate(fold_change=ifelse(Control>Disease, -(Control/Disease-1), (Disease/Control-1)),
#          enriched=ifelse(Control>Disease,"Control", "Disease"))
# 
# 
# ps = unique(dm$Project)
# unique(dm$Group)
# for(p in ps){
#   x = subset(sample_map, Project==p)
#   x_con = subset(x, Group=='Control')$Sample
#   x_dis = subset(x, Group == "Disease")$Sample
#   cons = dm %>% filter(`Row.names` %in% x_con) %>% select(alpha) %>% unlist()
#   diss = dm %>% filter(`Row.names` %in% x_dis) %>% select(alpha) %>% unlist()
#   pval = wilcox.test(cons,diss,exact = F)$p.value
#   data[data$Project == p, 'pval'] = pval
# }
# 
# data = merge(data, unique(project_map[,-2]), by='Project', all.y = F, all.x = T)
# data$shape = ifelse(data$pval < 0.01, "p<0.01", ifelse(data$pval < 0.05, "p<0.05", "p≥0.05"))
# write.table(data, "vOTU.shannon.tsv", row.names=FALSE, sep="\t", quote=FALSE)

#--------------
#   plot

data = read.table("vOTU.shannon.tsv",sep="\t", header=T, check.names=F)
orders = read.table("../00.data/project_order.list", sep="\t", header=T)
data$Project = factor(data$Project, levels=rev(orders$project))
load("../00.data/colors.map")
data$shape = ifelse(data$pval < 0.001, 
      "***", ifelse(data$pval < 0.01, 
      "**",  ifelse(data$pval < 0.05,
      "*",   NA)))
data$text_pos = ifelse(data$fold_change < 0, data$fold_change-0.02, data$fold_change+0.02)

shannon_foldchange <- ggplot(data, aes(y=Project, x=fold_change, color=group, fill=group))+
  geom_bar(stat='identity',width=0.8)+
  geom_text(aes(label=shape, x=text_pos), nudge_y=-0.35,size=8,color="black")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  ggtitle("Shannon index")
shannon_foldchange

# save(shannon_foldchange,file="ggplot_shannon_foldchange.RData")
