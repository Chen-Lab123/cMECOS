##Figure3A
specie_count_all_uhgg=read.delim("species_count_unculture_uhgg",header = F)
colnames(specie_count_all_uhgg)=c("Species","Count_all","Count_origin","Count_uhgg")
specie_count_all_uhgg=cbind(specie_count_all_uhgg,t(matrix(unlist(strsplit(specie_count_all_uhgg$Species,split = ";")),7)))
colnames(specie_count_all_uhgg)[5:11]=c("d","p","c","o","f","g","s")

specie_count_all_uhgg[specie_count_all_uhgg$Count_uhgg==0,"Count_uhgg"]=1
ggplot(specie_count_all_uhgg,aes(x=count_ori2,y=Count_uhgg,size=Count_all,col=p))+geom_jitter()+theme_classic()+geom_hline(yintercept = 20)+geom_vline(xintercept = 0.5)+scale_y_log10()+scale_color_manual(values=brewer.pal(7, "Set1"))

specie_count_all_uhgg$type="C"
specie_count_all_uhgg[specie_count_all_uhgg$Count_origin==0,"type"]="CE"

specie_count_all_uhgg$type_r="C"
specie_count_all_uhgg[specie_count_all_uhgg$Count_uhgg <= 20,"type_r"]="Rare"

specie_count_all_uhgg$type_a=paste(specie_count_all_uhgg$type,specie_count_all_uhgg$type_r,sep = "_")

ggplot(data.frame(table(specie_count_all_uhgg[,c("p","type_a")])),aes(fill=p,y=Freq,x=1))+geom_bar(stat="identity",position = "fill")+theme_classic()+scale_fill_manual(values=brewer.pal(7, "Set1"))+coord_polar(theta = "y")+facet_wrap(.~type_a)

##Figure3B
all_anno=read.delim("all_anno",header=F)
ac_nr_mat=acast(all_anno,V3~V5)
ac_nr_mat[ac_nr_mat>1]=1

ac_nr_mat=data.frame(ac_nr_mat)
ac_nr_mat$gene=rownames(ac_nr_mat)
upset(ac_nr_mat,nsets = 4,order.by = "freq")

##Figure3C
sp_gene_count=read.delim("sp_count.txt")
ggplot(sp_gene_count,aes(y=spieces ,x=count,fill=type))+geom_bar(stat="identity")+theme_classic()+ylim(subset(sp_gene_count,type=="cer_uniq")[order(subset(sp_gene_count,type=="cer_uniq")$count),]$spieces )

##Figure3D
count_cog_cer_uniq=read.delim("count_cog_rare_uniq",header=F)
count_cog_rm=read.delim("count_cog_rm",header=F)
count_cog=merge(count_cog_rm,count_cog_cer_uniq,by="V2")
colnames(count_cog)=c("cog_id","gene_all_rm","gene_cer_uniq")
count_cog$gene_all_rm_all=339179
count_cog$gene_cer_uniq_all=38762
count_cog$gene_all_rm_all_l=count_cog[,4]-count_cog[,2]
count_cog$gene_cer_uniq_all_l=count_cog[,5]-count_cog[,3]
count_cog$fd=(count_cog$gene_cer_uniq*count_cog$gene_all_rm_all_l)/(count_cog$gene_all_rm*count_cog$gene_cer_uniq_all_l )
count_cog$pvalue=apply(count_cog,1,function(a){fisher.test(matrix(as.numeric(c(a[2],a[3],a[6],a[7])),nrow=2))$p.value})
count_cog=merge(count_cog,cog_anno,by.x="cog_id",by.y="COG_id")
count_cog=count_cog[order(count_cog$fd),]
ggplot(count_cog,aes(y=cog_anno ,x=log2(fd),fill=class))+geom_bar(stat="identity")+ylim(as.character(count_cog$cog_id))+theme_classic()

##Figure3E
count_rm=read.delim("count_all_rm",header=F)
count_cer=read.delim("count_CER_uniq",header=F)
count_map=merge(count_rm,count_cer,by="V1",all = T)
count_map[is.na(count_map$V2.y),3]=0
colnames(count_map)=c("map_id","gene_all_rm","gene_cer_uniq")

count_map$gene_all_rm_all=339179
count_map$gene_cer_uniq_all=38762
count_map$gene_all_rm_all_l=count_map[,4]-count_map[,2]
count_map$gene_cer_uniq_all_l=count_map[,5]-count_map[,3]
count_map$pvalue=apply(count_map,1,function(a){fisher.test(matrix(as.numeric(c(a[2],a[3],a[6],a[7])),nrow=2))$p.value})
count_map$qvalue=p.adjust(count_map$pvalue,method = "fdr")
count_map$fd=(count_map$gene_cer_uniq/count_map$gene_cer_uniq_all)/(count_map$gene_all_rm/count_map$gene_all_rm_all )
count_map$log_qval=-log10(count_map$qvalue)
count_map[count_map$log_qval>5,"log_qval"]=5
count_map=merge(count_map,map_anno,by.x="map_id",by.y="kegg_pathway_id")
ggplot(subset(count_map,class2=="Metabolism"),aes(x=log2(fd),y=log_qval,col=class1))+geom_point()+geom_hline(yintercept = -log10(0.05))+theme_classic()+geom_vline(xintercept = c(log2(1.5),-log2(1.5)))+scale_color_manual(values=c(brewer.pal(6, "Set2"),brewer.pal(7, "Set3")))
##Figure3F
map_561=read.delim("map0561_cre_uniq_ko",header = F)
#map_561=read.delim("map0561_cre_uniq",header = F)
ac_map_561=acast(map_561,V2~V3,value.var = "V1",fill = 0)

ac_map_561[ac_map_561>1]=1
ce_rare_1genome=read.tree("ce_rare_msa_1_genome.tree")

plot(ce_rare_1genome,align.tip.label = T)
ce_rare_ge=read.delim("ce_rare_1_genome",header = F,row.names = "V2")


ac_map_561=rbind(ac_map_561,matrix(0,nrow=1,ncol=27))
rownames(ac_map_561)[23]=ce_rare_ge[ce_rare_1genome$tip.label,"V1"][!ce_rare_ge[ce_rare_1genome$tip.label,"V1"]%in% rownames(ac_map_561)  ]
pheatmap(ac_map_561[ce_rare_ge[ce_rare_1genome$tip.label,"V1"],],cluster_rows = F)


###Figure S10
##species for reads
all_sp13=read.csv("../all_13_sp_nozero_absolute.csv")
sample_sheet=read.delim("../sampleSheet_a.txt")
sample_sheet=subset(sample_sheet,seq_id %in% colnames(all_sp13))

##species for contig
cds_count_zym=read.delim("../202205/cds_count_zym_new",header = F)
matrix_contig_sp=acast(cds_count_zym,V2~V1,fill=0)

ac_bin_sp_pat=acast(bin_hq_gtdb,Species~Patient  )
#colnames(ac_bin_sp_pat)=colnames(matrix_contig_sp)

colnames(matrix_contig_sp)=colnames(ac_bin_sp_pat)

count_sp_pa_rarefraction=data.frame(patient_count=rep(1:13,each=10),rep_exp=rep(1:10,13))
count_sp_pa_rarefraction$count_sp=0
count_sp_pa_rarefraction$count_contig=0
count_sp_pa_rarefraction$count_bin=0
for (i in 1:130){
  select_pa= sample(unique(sample_sheet$PatientID),count_sp_pa_rarefraction[i,1])
  sp_count_pa=apply(all_sp13[,subset(sample_sheet,PatientID %in% select_pa)$seq_id],1,sum)
  count_sp_pa_rarefraction[i,3]=length(sp_count_pa[sp_count_pa>0])
  if(count_sp_pa_rarefraction[i,1]>1){
    contig_count_pa=apply(matrix_contig_sp[,select_pa],1,sum)
    bin_count_pa=apply(ac_bin_sp_pat[,select_pa],1,sum)
  }else{
    contig_count_pa=matrix_contig_sp[,select_pa]
    bin_count_pa=ac_bin_sp_pat[,select_pa]
  }
  count_sp_pa_rarefraction[i,4]=length(contig_count_pa[contig_count_pa>0])
  count_sp_pa_rarefraction[i,5]=length(bin_count_pa[bin_count_pa>0])
}
ggplot(count_sp_pa_rarefraction,aes(x=as.factor(patient_count),y=count_sp))+geom_boxplot()+theme_classic()

###Figure S12
count_mge=read.delim("count_mge",header = F)
count_res=read.delim("count_resfinder",header = F)
count_vfdb=read.delim("count_vfdb",header = F)
colnames(count_mge)=c("bin_ID","MGE_count")
colnames(count_res)=c("bin_ID","RES_count")
colnames(count_vfdb)=c("bin_ID","VF_count")

bin_hq_gtdb=merge(bin_hq_gtdb,count_mge)
bin_hq_gtdb=merge(bin_hq_gtdb,count_res)
bin_hq_gtdb=merge(bin_hq_gtdb,count_vfdb)

bin_hq_gtdb2=merge(bin_hq_gtdb,specie_count_all_uhgg,by.x=c("Phylum","Class","Order","Family","Genus","Species"),by.y=c("p","c","o","f","g","s"))

p1=ggplot(subset(bin_hq_gtdb2,type_a=="CE_Rare"),aes(y=Species,x=MGE_count))+geom_bar(stat="identity")+theme_classic()
p2=ggplot(subset(bin_hq_gtdb2,type_a=="CE_Rare"),aes(y=Species,x=RES_count))+geom_bar(stat="identity")+theme_classic()
p3=ggplot(subset(bin_hq_gtdb2,type_a=="CE_Rare"),aes(y=Species,x=VF_count))+geom_bar(stat="identity")+theme_classic()
grid.arrange(p1,p2,p3,nrow=1)