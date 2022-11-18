library(tidyverse)

source("../scripts/zy_plot_ROC.R")
dt = read.table("./fig8.tsv", sep="\t", header=T)
dt$group = paste(dt$nspecies, dt$all, sep="_")

x <- calc_auc(dt, pred="Disease", true="Group",group="group")
head(x$table)
data = x$table
data$name = rownames(data)
data = separate(data, col=name, into=c("nspecies","risk"), sep="_")
data$nspecies = as.numeric(data$nspecies)

p <- ggplot(data, aes(x=nspecies, y=auc, color=risk))+
  geom_point()+
  geom_line()+
  scale_x_continuous(trans="log10")+
  theme_bw()

x <- plot_roc(dt %>% filter(nspecies == 7161), pred="Disease", true="Group",group="all")
x$plot

library(patchwork)
p+x$plot+
  plot_layout(widths=c(2,1))
