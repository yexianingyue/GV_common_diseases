world <- ne_countries(scale = "medium", returnclass = "sf")
sites <- data.frame(longitude = c(-80.144005, -80.109), latitude = c(26.479005, 26.83))

ggplot(data = world) +
geom_sf() +
geom_point(data = sites, aes(x = longitude, y = latitude), size = 4, 
           shape = 23, fill = "darkred") +
#coord_sf(xlim = c(-88, -78), ylim = c(24.5, 33), expand = FALSE)
coord_sf(xlim = c(73.33, 135.05), ylim = c(3.51, 53.33), expand = FALSE)


sf::sf_use_s2(FALSE)


china = sf::read_sf("china.geojson")
info = read.csv('map.txt',sep='\t',stringsAsFactors = F)
aa = aggregate(info$sample,info[,c('province','lon_province','lat_province','proj')],sum)
area <- subset(china, china$name%in%aa$province)
area$Group = aa$x[match(area$name,aa$province)]

ggplot(world) +
geom_sf(fill = "antiquewhite1",color='#a2a3a5',size=0.1) +
geom_sf(data = fortify(china),fill = 'antiquewhite1', color = 'grey',size=0.1) +
geom_sf(data = area,aes(fill = Group), color = NA) +
#scale_fill_viridis_c(trans = "log10", alpha = .4) +
#scale_fill_gradientn(colors=brewer.pal(10,'Spectral'),trans = "log10") +
scale_fill_viridis_c(trans = "log10", alpha = .4) +
#scale_fill_gradient(low = 'grey90',high = 'steelblue',trans = "log10") +
#geom_point(data = aa, aes(x =lat_province , y = lon_province,size=x),color='red',alpha=.4)+
geom_label_repel(data = aa, aes(x =lat_province , y = lon_province, label = proj), size=3,force = 0.3,force_pull = 0.5)+
coord_sf(xlim = c(75, 135), ylim = c(18, 55))+
#coord_sf()+
annotation_scale(location = "bl", width_hint = 0.2) +
xlab("Longitude") + ylab("Latitude")+
theme(
      # panel.grid.major = element_line(
      #   color = gray(0.5),
      #   linetype = "dashed",
      #   size = 0.5
      # ),
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "aliceblue")
      )



china = sf::read_sf("china.geojson")
info = read.csv('map2.txt',sep='\t',stringsAsFactors = F)
aa = aggregate(info$sample,info[,c('lab','longitude','latitude')],sum)

ggplot(data=world) +
geom_sf(fill = "#dddddd",color='#dddddd',size=0.1) +
geom_sf(data = fortify(china),fill = 'antiquewhite1', color = 'antiquewhite1',size=0.1) +
#geom_sf(data = area,aes(fill = Group), color = NA) +
#scale_fill_viridis_c(trans = "log10", alpha = .4) +
#scale_fill_gradientn(colors=brewer.pal(10,'Spectral'),trans = "log10") +
scale_fill_gradient(low = 'grey90',high = 'steelblue',trans = "log10") +
geom_point(
           data = aa,
           aes(x = longitude, y = latitude, size = x),
           color = '#d53e4f',
           #color = 'black',
           alpha = .6
           )+
#geom_text_repel(data = aa, aes(x =lat_province , y = lon_province, label = proj), 
#                fontface = "bold", nudge_x = c(1, -1.5, 2, 2, -1), size=3,
#                nudge_y = c(0.25,-0.25, 0.5, 0.5, -0.5))+
coord_sf(xlim = c(75, 135), ylim = c(18, 55))+
annotation_scale(location = "bl", width_hint = 0.2) +
#coord_sf()+
xlab("Longitude") + ylab("Latitude")+
theme(
      # panel.grid.major = element_line(
      #   color = gray(0.5),
      #   linetype = "dashed",
      #   size = 0.5
      # ),
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "aliceblue")
      )

#物种稀释曲线================================================================================
prof = read.csv('vOTU.clstr',sep='\t',header = F,stringsAsFactors = F,colClasses = 'character')
prof = tidyr::separate_rows(prof,'V3', sep = ",")[,c(1,3)]
colnames(prof) = c('speClstr','contigID')
prof$speClstr = as.numeric(as.factor(prof$speClstr)) 
prof$contigID = as.numeric(as.factor(prof$contigID)) 
ff = plyr::count(prof$speClstr)
ff = ff[ff$freq>1,]

n = max(prof$contigID)
dat = data.frame(N=c(rep(c(seq(1,n,1e3)),each=3)),all=NA,nonsingle=NA,stringsAsFactors = F)
for(x in 1:nrow(dat)){
    Sa = sample(1:n,size = dat$N[x],replace = F)
    dat$all[x]=length(unique(prof$speClstr[prof$contigID%in%Sa]))
    dat$nonsingle[x]=length(unique(prof$speClstr[prof$contigID%in%Sa & prof$speClstr%in%ff$x]))
}
speRaf = aggregate(dat[,c('all','nonsingle')],list(Group=dat$N),mean)
speRaf2 = aggregate(dat[,c('all','nonsingle')],list(Group=dat$N),sd)
speRaf$a_up = speRaf$all+speRaf2$all
speRaf$a_low = speRaf$all-speRaf2$all
speRaf$s_up = speRaf$nonsingle+speRaf2$nonsingle
speRaf$s_low = speRaf$nonsingle-speRaf2$nonsingle
f2=ggplot(speRaf)+
geom_line(aes(x=Group,y=all))+
geom_line(aes(x=Group,y=nonsingle),color='#888888')+
#geom_errorbar(aes(x=Group,ymin=a_low,ymax=a_up))+
#geom_errorbar(aes(x=Group,ymin=s_low,ymax=s_up))+
scale_y_continuous(breaks = seq(0,7e4,1.5e4))+
theme_test()
f2
ggarrange(f1,f2,ncol=2)
d.csv('votu.tax.mapping20220613',sep='\t',stringsAsFactors = F)
p1 = read.csv('21.vOTU/vOTU.fna.len',sep='\t',stringsAsFactors = F,header = F)
p2 = read.csv('20.new_catalog/final.vfa.len',sep='\t',stringsAsFactors = F,header = F)
p1$G = 'votu'
p2$G = 'all'


p = rbind(p1,p2)

cutoff = nrow(p2)/nrow(p1)
f1=ggplot()+
geom_histogram(data=p1,aes(x=V2,y=..count..*cutoff),fill='red',bins=50,alpha=0.3)+
geom_histogram(data=p2,aes(x=V2,y=..count..),fill='blue',bins=50,alpha=0.3)+
scale_x_continuous(trans = 'sqrt',breaks=c(1e4,1e5,2e5,3e5,4e5,5e5))+
scale_y_continuous(
                   sec.axis = sec_axis( ~./cutoff, #对次坐标轴刻度标签的二次映射（极差范围指定真实极差即可）  
                                       name = "Dentity")           #次坐标名
                   )+
theme_test()



aa = plyr::count(vtax$Family)
aa = aa$x[head(order(-aa$freq),n=15)]
vtax$Family[!vtax$Family%in%aa] = 'Others'
vtax$Family = factor(vtax$Family,c(aa,'Others'))

f2=ggplot(vtax)+
geom_density(aes(x=ContigSize,y=..density..,color=Family),alpha=0.3)+
scale_x_continuous(trans='sqrt',breaks=c(1e4,1e5,2e5,3e5,4e5,5e5))+
scale_color_manual(values = c(brewer.pal(12,'Paired'),brewer.pal(8,'Set2')))+
guides(color=F)+
theme_test()

ggarrange(f1,f2,nrow = 2,align = 'hv')


vtax$Comp = vtax$ContigSize/(ckv$V10[match(vtax$ID,ckv$V1)]/100)
f3=ggplot(vtax)+
geom_density(aes(x=Comp,y=..density..,color=Family),alpha=0.3)+
scale_x_continuous(trans='sqrt',breaks=c(1e4,1e5,2e5,3e5,4e5,5e5))+
scale_color_manual(values = c(brewer.pal(12,'Paired'),brewer.pal(8,'Set2')))+
guides(color=F)+
theme_test()
ggarrange(f1,f3,nrow = 2,align = 'hv')

