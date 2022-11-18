library(ComplexHeatmap)
library(dplyr)
library(tidyselect)
library(reshape2)
library(circlize)
rm(list=ls())

# define function for heatmap
sigFun <- function(x){
  if(is.na(x)){NA}
  if(x<0.01){"+"
  }else if(x<0.05){
    "*"
  }
}

load("colors.map")

marker = read.table("enriched_vOTU_family.tsv",sep="\t", check.names=F, header=T)
dt = read.table("wilcox_vOTU_family.tsv",sep="\t", check.names=F, header=T)
project_map = read.table("project.group",sep="\t", header=T, row.names=1, check.names=F)
dt$log2_fd = ifelse(dt$enriched == "Disease", 0-log2(dt$fold_change), log2(dt$fold_change))
project_ord = read.table("project_order.list", header=T)

pq <- data.frame(dcast(dt, name~Project, value.var="qval"), row.names=1, check.names=F)
pv <- data.frame(dcast(dt, name~Project, value.var="pvalue"), row.names=1, check.names=F)
fd <- data.frame(dcast(dt, name~Project,value.var='log2_fd'), row.names=1, check.names=F)

pq[is.na(pq)] = 1
pv[is.na(pv)] = 1

#filter
# marker = marker %>% filter(known != "unknown")
#dt <- dt %>% filter(name %in% marker$name)




# plot data
fd = fd[marker$name,]
pq = pq[marker$name,colnames(fd)]

# 聚类
rod = hclust(dist(fd,method = "canb"),method = "complete")$order
# cod = hclust(dist(t(fd),method = "maximum"),method = "complete")$order
cod = project_ord$project

fd = fd[rod, cod]
pq = pq[rod, cod]

# Heatmap color map
col_legend = colorRamp2(c(-1,-0.5,0, 0.5,1),
                        c("#8E0152","#D5589D", "#ffffff","#6EAE36","#276419")
)
# significant barplot
barplot_dt = matrix(0,ncol=2, nrow=nrow(fd),
                    dimnames = list(
                      rownames(fd),
                      c("Control","Disease")
                    ))

for( i in rownames(pq)){
  cc = sum((pq[i,fd[i,]>0]<0.05)+0) # count -> Control
  cd = sum((pq[i,fd[i,]<0]<0.05)+0) # count -> Disease
  barplot_dt[i,"Control"] = cc
  barplot_dt[i,"Disease"] = cd
}

# left annotation
row_ha = rowAnnotation(
  bar=anno_barplot(barplot_dt)
)

# 分割
rownames(marker) = marker$name
marker = marker[rownames(fd),]
rsp = data.frame(enriched=marker$enriched)
rsp = factor(rsp$enriched, levels=c("Control","Disease"))

# heatmap
cm = colnames(fd)
Heatmap(as.matrix(fd),
        row_split = rsp, border=T,
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = T,
        row_names_gp = gpar(fongsize=1),
        column_names_gp = gpar(col = colors[project_map[cm,'group']]),
        use_raster=F,
        col = col_legend,
        right_annotation = row_ha,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", sigFun(pq[i, j])), x, y=ifelse(pq[i, j]<0.01,y,y+unit(-0.015,"npc")), 
                    gp = gpar(fontsize = ifelse(pq[i, j]<0.01,15,20))
          )
          # unit color
          grid.rect(x = x, y = y, width = width, height = height, 
                    gp = gpar(col = "grey", fill = NA))
        }
)
# 5 11
#heatmap.enriched_genus.pdf
