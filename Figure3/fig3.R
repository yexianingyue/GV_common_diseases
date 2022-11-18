#基因稀释曲线================================================================================
prof = read.csv('01.geneclstr/all.virus.AA.i90_c80_cluster.tsv',sep='\t',header = F,stringsAsFactors = F,colClasses = 'character')
prof$V1 = as.numeric(as.factor(prof$V1)) 
prof$V2 = sub('_\\d+$','',prof$V2)
prof$V2 = as.numeric(as.factor(prof$V2)) 
prof = unique(prof)
colnames(prof) = c('geneClstr','contigID')
ff = plyr::count(prof$geneClstr)
ff = ff[ff$freq>1,]

n = max(prof$contigID)
dat = data.frame(N=c(rep(c(seq(4e2,3e4,3e3),seq(4e4,n,2e4)),each=3)),all=NA,nonsingle=NA,stringsAsFactors = F)
for(x in 1:nrow(dat)){
    Sa = sample(1:n,size = dat$N[x],replace = F)
    dat$all[x]=length(unique(prof$geneClstr[prof$contigID%in%Sa]))
    dat$nonsingle[x]=length(unique(prof$geneClstr[prof$contigID%in%Sa & prof$geneClstr%in%ff$x]))
}

#rm(prof)
genRaf = aggregate(dat[,c('all','nonsingle')],list(Group=dat$N),mean)
f1=ggplot(genRaf)+
geom_line(aes(x=Group,y=all))+
geom_line(aes(x=Group,y=nonsingle),color='#888888')+
#geom_errorbar(aes(x=Group,ymin=a_low,ymax=a_up))+
#geom_errorbar(aes(x=Group,ymin=s_low,ymax=s_up))+
scale_y_continuous(breaks = seq(0,2e6,4e5))+
theme_test()
f1

#高质量病毒的科水平kegg功能比较============================================================================
# gg = c('Siphoviridae','Myoviridae','Microviridae','Podoviridae','Retroviridae',
#        'Quimbyviridae','Salasmaviridae','Autographiviridae','Inoviridae','Metaviridae',
#        'Podoviridae_crAss-like','Others')


mkegg =  read.csv('/share/data2/guorc/Database/KEGG/20201208/all_ko.mapping',sep='\t',stringsAsFactors = F)
mkegg = mkegg[mkegg$LevelA=="Metabolism",]

ckv = read.csv('00.Data/vOTU.ckv',header = F,sep='\t',stringsAsFactors = F)
prof = read.csv('05.fun/votu_kegg.m8',header = F,sep='\t',stringsAsFactors = F)
vtax = read.csv('../mapping.file/votu.tax.mapping20220613',sep='\t',stringsAsFactors = F)
#vtax = vtax[vtax$ID%in%ckv$V1[ckv$V8%in%c('Complete','High-quality')],]
vtax$ngene = ckv$V5[match(vtax$ID,ckv$V1)]
allgene = aggregate(vtax$ngene,list(tax=vtax$Family),sum)
allsp = plyr::count(vtax$Family)
allsp = allsp[order(-allsp$freq),]

#gg = c(Annogene$x[head(order(-Annogene$freq),n=10)],'Others')
gg = c(allsp$x[1:11],'Others')
allsp$x[12:nrow(allsp)] = 'Others'
allsp = aggregate(allsp$freq,list(tax=allsp$x),sum)
allgene$tax[!allgene$tax%in%gg] = 'Others'
allgene = aggregate(allgene$x,list(allgene$tax),sum)


prof$Contig = sub('_\\d+$','',prof$V1)
prof = prof[prof$Contig%in%vtax$ID,]
prof$tax = vtax$Family[match(prof$Contig,vtax$ID)]
prof$ko = sub('.*\\|','',prof$V2)

allanno = plyr::count(prof$tax)
allanno$x[!allanno$x%in%gg] = 'Others'
allanno = aggregate(allanno$freq,list(allanno$x),sum)

lab='myLevelA'
mkeggf = read.csv('../mapping.file/metabolism_ko.mapping3',skip=1,sep='\t',stringsAsFactors = F)
mkeggf = mkeggf[!mkeggf[,lab]%in%c('Protein families: metabolism',"Others","Unclassified: metabolism"),]
mkeggf$myLevelA[mkeggf$myLevelB=='Peptidoglycan biosynthesis and degradation proteins'] = 'Peptidoglycan biosynthesis and degradation'
p = merge(prof[,c('V1','Contig','tax','ko')],unique(mkeggf[,c(lab,'KO')]),by.x='ko',by.y='KO')
colnames(p)[5] = 'Group' 
aa = plyr::count(p$V1)
#p$val = 1/aa$freq[match(p$V1,aa$x)]
p$val = 1
pf = aggregate(p$val, list(Group=p$Group), sum)
pf$rate = pf$x/sum(allanno$x)
#pf$rate = pf$x/sum(allgene$x)
pf= pf[order(pf$rate),]
Or = pf$Group
pf$Group = factor(pf$Group,Or) 
f1=ggplot(pf)+
geom_bar(aes(x=rate,y=Group),fill='#aaaaaa',color='black',position='dodge',stat='identity',width = .75)+
#scale_fill_manual(values = colorRampPalette(brewer.pal(11,'Spectral'))(11))+
guides(fill=F)+
theme_test()


pf = p
pf$tax[!pf$tax%in%gg] = 'Others'
pf = aggregate(pf$val, list(Group=pf$Group,tax=pf$tax), sum)
pf$rate = pf$x/allanno$x[match(pf$tax,allanno$Group.1)]
#pf$rate = pf$x/allgene$x[match(pf$tax,allgene$Group.1)]
pf$Group = factor(pf$Group,Or)
pf$tax = factor(pf$tax,gg)
#write.table(pf,'../../AMG.list',sep='\t',col.names = F,row.names = F,quote = F)
f2=ggplot(pf)+
geom_tile(aes(x=tax,y=Group,fill=rate))+
scale_fill_gradientn(trans='log2',colors = c('white',brewer.pal(11,'Spectral')[6:1]),
                     breaks=c(0.001,0.01,0.1,0.5))+
#breaks=c(0.0001,0.001,0.003,0.01))+
#guides(fill=F)+
theme_test()+
theme(
      axis.text.x = element_text(
                                 angle=45,
                                 hjust = 1,
                                 vjust = 1
                                 )
      )

ggarrange(f1,f2,nrow=1,align = 'hv',widths = c(1.3,2))

#高质量病毒的科水平cazy功能比较============================================================================
Top = 20

ckv = read.csv('00.Data/vOTU.ckv',header = F,sep='\t',stringsAsFactors = F)
prof = read.csv('05.fun/votu_cazy.m8',header = F,sep='\t',stringsAsFactors = F)
vtax = read.csv('../mapping.file/votu.tax.mapping20220613',sep='\t',stringsAsFactors = F)
#vtax = vtax[vtax$ID%in%ckv$V1[ckv$V8%in%c('Complete','High-quality')],]
vtax$ngene = ckv$V5[match(vtax$ID,ckv$V1)]
allgene = aggregate(vtax$ngene,list(tax=vtax$Family),sum)
allsp = plyr::count(vtax$Family)
allsp = allsp[order(-allsp$freq),]


prof$Contig = sub('_\\d+$','',prof$V1)
prof = prof[prof$Contig%in%vtax$ID,]
prof$tax = vtax$Family[match(prof$Contig,vtax$ID)]
prof$cazy = sub('\\|$|\\|\\d+.*$','',prof$V2)
prof$cazy = sub('.*\\|','',prof$cazy)
prof$cazyid = prof$cazy
prof$cazy = gsub('_|\\d+','',prof$cazy)


allanno = plyr::count(prof$tax)
allanno = allanno[order(-allanno$freq),]
gg = c(allanno$x[1:Top],'Others')
allanno$x[!allanno$x%in%gg] = 'Others'
allanno = aggregate(allanno$freq,list(allanno$x),sum)

#gg = c(allsp$x[1:Top],'Others')
#allsp$x[(Top+1):nrow(allsp)] = 'Others'
allsp$x[!allsp$x%in%gg] = 'Others'
allsp = aggregate(allsp$freq,list(tax=allsp$x),sum)
allgene$tax[!allgene$tax%in%gg] = 'Others'
allgene = aggregate(allgene$x,list(allgene$tax),sum)



Or = plyr::count(prof$cazyid)
Or = Or[tail(order(Or$freq),n=30),]
pf = prof[prof$cazyid%in%Or$x,c('tax','cazyid')]
pf$tax[!pf$tax%in%gg] = 'Others'
pf = plyr::count(pf)
pf$cazyid = factor(pf$cazyid,Or$x)
pf$tax = factor(pf$tax,c(setdiff(gg,'Unknown'),'Unknown'))
#write.table(p,'../../antibio.list',sep='\t',col.names = F,row.names = F,quote = F)
ggplot(pf)+
geom_bar(aes(x=freq,y=cazyid,fill=tax),position='stack',stat='identity',width = .75)+
scale_fill_manual(values = brewer.pal(12,'Paired'))+
#guides(fill=F)+
theme_test()
