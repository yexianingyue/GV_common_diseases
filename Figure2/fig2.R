library(UpSetR)

vtax = read.csv('04.tax/all.tax_family2',header = F,sep='\t',stringsAsFactors = F)

hqckv  = read.csv('02.dbCom/dbCom.ckv',sep='\t',header = F,stringsAsFactors = F)
hqckv = hqckv[hqckv$V8=='Complete',]

map = read.csv('02.dbCom/dbcom.list',sep='\t',header = F,stringsAsFactors = F)
prof = read.csv('02.dbCom/dbcom.i95_c70.uniq.list',sep='\t',header = F,stringsAsFactors = F)
colnames(prof) = c('clstr','V1')
#prof = prof[prof$V1%in%hqckv$V1,]
prof$Group = map$V1[match(prof$V1,map$V2)]
prof = prof[!is.na(prof$Group),]
#prof = merge(prof,dat[,c('V1','Group')],by='V1',all.x = T)
prof$value = 1

#dat = prof[prof$clstr%in%prof$clstr[prof$Group=='vOTU'],]
dat = prof
aa = plyr::count(dat$clstr)
dat$Uniq = aa$freq[match(dat$clstr,aa$x)]
dat = dat[dat$Group=='vOTU',]
dat$Uniq=ifelse(dat$Uniq==1,'Uniq','Share')
#write.table(dat,'aa',sep='\t',quote = F,row.names = F)

dat = dcast(data = prof[,c('clstr','Group')],clstr~Group)
row.names(dat) = dat$clstr
dat = dat[,-1]
#write.table(dat,'aa',sep='\t',quote = F)
dat[dat!=0] = 1


gg=c("vOTU","CHVD","GPD","GVD","MGV")
aa = c()
for(x in gg){
    aa[[x]] = row.names(dat)[dat[,x]!=0]
}
fig = venn.diagram(
                   aa,
                   #list(votu=aa$votu,GPD=aa$GPD,GVD=aa$GVD,MGV=aa$MGV),
                   #fill=colorRampPalette(brewer.pal(8,'Spectral'))(5), 
                   alpha=c(0.5,0.5,0.5,0.5,0.5), 
                   col = colorRampPalette(brewer.pal(8,'Spectral'))(5),
                   cex=0.9,
                   cat.pos = 6,
                   cat.dist = 0.05,
                   filename = NULL)

grid.draw(fig)

upset(dat, sets = gg, mb.ratio = c(0.6, 0.4), order.by = "freq", 
      number.angles = 90,point.size = 2, line.size = 1, mainbar.y.label = "title",
      sets.x.label = "title", text.scale = c(1, 1, 1, 1, 1, 1))
colSums(dat)

#统计我们的数据库中不同科水平的高质量基因组共有情况=====
hq = read.csv('02.dbCom/hq.list',sep='\t',header = F)
vtax = read.csv('04.tax/all.tax_family2',header = F,sep='\t',stringsAsFactors = F)
#dat = prof[prof$clstr%in%prof$clstr[prof$Group=='vOTU'],]
dat = prof[prof$V1%in%hq$V1 & prof$Group!='RefSeq',]
#dat = prof[prof$V1%in%hq$V1,]
dat$tax = vtax$V4[match(dat$V1,vtax$V1)]
aa = unique(dat[!is.na(dat$tax),c('clstr','tax')])
dat = dat[dat$clstr%in%aa$clstr,]

aa =  plyr::count(dat[!is.na(dat$tax),c('clstr','tax')])

aa = aa[order(-aa$freq),]
aa = aa[!duplicated(aa$clstr),]

bb = plyr::count(aa$tax)
bb = bb[order(-bb$freq),]
bb$tax = 'Others'
bb$tax[1:15] = bb$x[1:15]
#bb$tax = bb$x

aa$tax2 = bb$tax[match(aa$tax,bb$x)]

# for(x in unique(aa$tax2)){
#   write.table(aa$clstr[aa$tax2==x],paste('08.db_PD_new/hq.',x,'.list',sep=''),quote=F,row.names = F,col.names = F)
# }

# bb=aa$clstr[aa$tax2=='Siphoviridae']
# for(x in 1:10){
#  write.table(sample(bb,size = 1500,replace = F),paste('08.db_PD_new/hq.Siphoviridae/list.',x,sep=''),quote=F,row.names = F,col.names = F)
# }

dat$tax = aa$tax2[match(dat$clstr,aa$clstr)]
sum(dat$Group=='vOTU')
datf = dcast(data=dat[,c('clstr','tax','Group')],clstr+tax~Group)
pd_prof = datf[,3:ncol(datf)]
row.names(pd_prof) = gsub('\\|','_',datf$clstr)

datf$Uniq = 'Uniq'

#datf$Uniq[apply(datf[,c('GPD','GVD','MGV','RefSeq','CHVD')],1,sum)!=0] = 'Share'
datf$Uniq[apply(datf[,c('GPD','GVD','MGV','CHVD')],1,sum)!=0] = 'Share'
datf$Uniq[datf$Uniq == 'Share' & datf$vOTU == 0] = 'pre'

Matrix_overlap_otu = datf

#datf = aggregate(datf$vOTU,datf[,c('tax','Uniq')],sum)
datf = plyr::count(datf[,c('Uniq','tax')])


aa = aggregate(datf$freq,list(tax=datf$tax),sum)
aa = aa[order(-aa$x),]
gg = unique(bb$tax)
#aa$tax2[1:20] = aa$tax[1:20]
#datf$Rate = datf$freq/aa$x[match(datf$tax,aa$tax)]
datf$tax[!datf$tax%in%gg] = 'Others'
datf = aggregate(datf$freq,datf[,c('Uniq','tax')],sum)
datf$tax =factor(datf$tax,gg)

datf$Group = 'G1'
#datf$Group[datf$tax%in%gg[4:12]] = 'G2'
datf$Group[datf$tax%in%gg[5:length(gg)]] = 'G3'
datf$Uniq =factor(datf$Uniq,c('Uniq','Share','pre'))

Lab = datf[datf$Uniq=='Uniq',]
datf = datf#[datf$tax%in%Lab$tax,]
f2=ggplot(datf)+
geom_bar(aes(x=tax,y=x,fill=Uniq),position = 'stack',stat='identity',width = 0.70)+
#geom_text(data=Lab,aes(x=tax,y=x,label=x),color='red',size=3,angle=90)+
scale_y_continuous(limits =c(0,1500))+
#scale_fill_manual(values = c('black','white'))+
facet_grid(~Group,scales = 'free',space = 'free')+
theme_classic()+
theme(
      axis.text.x = element_text(
                                 size=7,
                                 angle=90,
                                 hjust=1
                                 ),
      legend.key.size = unit(0.5, "cm"),
      #axis.line = element_blank(),
      axis.text.y = element_text(
                                 size=7,
                                 angle=90
                                 )
      )
f1
ggarrange(f1,f2,nrow = 2,align = 'hv')
