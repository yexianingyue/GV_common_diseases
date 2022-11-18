#marker的前50的ko功能比较============================================================================
Fun = read.csv('aa.txt',sep='\t',stringsAsFactors = F)
mfun = read.csv('bb.txt',sep='\t',header = F,stringsAsFactors = F)

dd = melt(Fun[,c('Kos','X..Control','X..Disease')])
colnames(dd) = c('Anno','Enrich','Rate')
dd$Anno = factor(dd$Anno,rev(Fun$Kos))

Fun$lab = ''
Fun$lab[Fun$qvalue<0.05] = '*'
Fun$lab[Fun$qvalue<0.01] = '**'
Fun$lab[Fun$qvalue<0.001] = '***'
Fun$Pos = as.numeric(apply(Fun[,4:5],1,max))+2
f1=ggplot(dd)+
geom_bar(aes(x=Rate,y=Anno,fill=Enrich),stat='identity',position = 'dodge',width = 0.6)+
geom_text(data=Fun,aes(x=Pos,y=Kos,label=lab),size=3)+
theme_bw()+
theme(
      axis.text.x = element_text(
                                 angle=90,
                                 hjust=1,
                                 size=7,
                                 vjust=0.5
                                 ),
      axis.text.y= element_text(
                                size=7
                                )
      )


mfunf = mfun[mfun$V1%in%Fun$Kos,]
mfunf$val = 1
mfunf$V1 = factor(mfunf$V1,rev(Fun$Kos))
f2=ggplot(mfunf)+
geom_bar(aes(x=val,y=V1,fill=V3),stat='identity',position = 'stack',width = 0.6)+
scale_fill_manual(values = c(brewer.pal(12,'Paired'),brewer.pal(8,'Set2'),brewer.pal(9,'Set1')))+
theme_transparent()+
theme(
      axis.text = element_blank()
      )

ggarrange(f1,f2,align = 'hv',widths = c(1,0.9))
#疾病marker的分布============================================================================
markerlist = read.csv('aa.txt',sep='\t',stringsAsFactors = F,)

p = read.csv('06.host/vOTU.host.tax',sep='\t',header=F,check.names = F)
vtax = read.csv('../mapping.file/votu.tax.mapping20220613',sep='\t',stringsAsFactors = F)
p = p[,c('V1','V3','V6','V7')]
colnames(p)[2:4] = c('phylum','family','genus')
p$phylum = sub('_[A-Z]$','',p$phylum)
p = merge(vtax[,c('ID','Family')],p,by.x = 'ID',by.y='V1',all.x = T)
p[is.na(p)] = 'Unknown'

p = p[p$ID%in%markerlist[,1],]
p$Group = markerlist[match(p$ID,markerlist[,1]),2]

lab = 'family'
dat = plyr::count(p[,c('ID','Family','Group',lab)])
ff = plyr::count(dat$ID)
ff = ff$x[ff$freq>1]
a1 = unique(dat[dat$ID%in%ff,c('ID','Family','Group')])
a1[,lab] = 'Multiple'
a2 = dat[!dat$ID%in%ff,c('ID','Family','Group',lab)]
dat = rbind(a1,a2)
#write.table(dat,'marker.tax_host.list',sep='\t',quote = F,row.names = F)
dat = plyr::count(dat[,2:4])

n = plyr::count(markerlist[,2])
#dat$freq = dat$freq/n$freq[match(dat$Group,n$x)]

aa = aggregate(dat$freq,list(dat$Family),sum)
aa = aa$Group.1[head(order(-aa$x),n=6)]
dat$Family[!dat$Family%in%aa] = 'Others'
aa = aggregate(dat$freq,list(dat[,lab]),sum)
aa = aa$Group.1[head(order(-aa$x),n=15)]
dat[!dat[,lab]%in%aa,lab] = 'Others'
dat = aggregate(dat$freq,dat[,c('Family','Group',lab)],sum)
dat = dat[order(-dat$x),]

aa=dcast(family+Group~Family,data=dat)
aa[is.na(aa)]=0
#write.table(aa,'bb',sep='\t',quote = F,row.names = F)

Or = aggregate(dat$x,list(dat$Family),sum)
Or = Or$Group.1[order(-Or$x)]
dat$Family = factor(dat$Family,Or)
Or = aggregate(dat$x,list(dat[,lab]),sum)
Or = Or$Group.1[order(-Or$x)]
dat[,lab] = factor(dat[,lab],Or)

Col = data.frame(ID=sort(unique(dat[,lab])),Col=c(brewer.pal(12,'Paired'),brewer.pal(12,'Set2'),
                                                  brewer.pal(9,'Set1'))[1:length(unique(dat[,lab]))])
dat$Col = Col$Col[match(dat[,lab],Col$ID)]
dat = dat[order(-dat$x),]
dat$Fill = paste(dat$Family,dat$Group,dat[,lab],sep='_')
dat$Fill = factor(dat$Fill,dat$Fill)


f1=ggplot(dat)+
geom_bar(aes(x=Family,y=x,fill=Fill),position = 'stack',
         stat = 'identity',width = 0.6,color=NA)+
#scale_y_continuous(limits =c(0,3000))+
#scale_y_continuous(limits =c(0,130))+
scale_fill_manual(values = dat$Col)+
facet_grid(~Group,scales = 'free',space = 'free')+
guides(fill=F)+
theme_bw()+
theme(
      axis.text.x = element_text(
                                 #size=5,
                                 angle=45,
                                 hjust=1
                                 ),
      legend.key.size = unit(0.5, "cm"),
      axis.line = element_blank(),
      axis.text.y = element_text(
                                 #size=5,
                                 angle=45
                                 )
      )

f2 = ggplot(Col)+
geom_tile(aes(x=1,y=ID,fill=ID))+
#facet_grid(ph~.,scales = 'free',space = 'free')+
scale_fill_manual(values = Col$Col)

ggarrange(f1,f2,f3,f4,align = 'h')
