
if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
}
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)
# source
source('z:/projects/codes/mg_base.R')
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(text = element_text(family = 'Times'),panel.grid = element_blank())
}
my_volcano=function(dat,p_cutoff=0.05,fc_cutoff=1,col=c("red","blue","black"),
                    ylab='-log10 (adj.PVal)',xlab='log2 (FoldChange)',leg.pos='right'){
  degs_dat=dat$DEG
  degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                              ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
  p=ggplot(degs_dat,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
    geom_point()+
    scale_color_manual(values=col)
    # geom_text_repel(
    #   data = tcga.diff$DEG[tcga.diff$DEG$adj.P.Val<p_fit & abs(tcga.diff$DEG$logFC)>fc_fit,],
    #   #aes(label = Gene),
    #   size = 3,
    #   segment.color = "black", show.legend = FALSE )
    theme_bw()
    theme(
      legend.title = element_blank(),
      legend.position = leg.pos,
    )+
    ylab(ylab)
    xlab(xlab)
    geom_vline(xintercept=c(-fc_cutoff,fc_cutoff),lty=3,col="black",lwd=0.5)|FoldChange|>2
    geom_hline(yintercept = -log10(p_cutoff),lty=3,col="black",lwd=0.5)padj<0.05
  return(p)
}

bioForest=function(rt=null,col){

  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  

  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  

  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  

  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
  
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  

  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) 
  return(geneList)
}
mg_nomogram=function(clinical_riskscore,os,status,title='Nomogram',
                     quick=T,mks = c(1,3,5)){
  #clinical_riskscore=dat1[,3:5]
  #os=dat1[,1]
  #status=dat1[,2]
  #sum(is.na(norm.stat.al[,3]))
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  #summary(cox2)
  #surv=Survival(cox2)
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  #cindex.orig=1-rcorr.cens(fp,Surv(norm.stat.al$time,norm.stat.al$status))
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=mks
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*mks[1],12*mks[2],12*mks[3])
  }else{
    cut.time=c(365*mks[1],365*mks[2],365*mks[3])
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['100%'])]
  print(cut.time)
  #regplot(cox2)
  #  print(regplot(cox3
  #,observation=pbc[2,] 

  #              ,title=title
  #              ,failtime = cut.time
  #              ,prfail = TRUE #cox
  #              ,showP = T 
  #              ,droplines = F
  #,colors = mg_colors[1:3] 
  #,rank="decreasing") 
  #,interval="confidence"
  #,rank="decreasing"
  #,clickable=T
  #              ,points=TRUE)) 
  
  #  plot(nom)
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
    #    surv0 <- function(x)surv(cut.time[i],lp=x) 
    #    survs=c(survs,list(surv0))
  }
  #f2<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[2]) 
  #cal3<-calibrate(f2, cmethod="KM", method="boot",u=cut.time[2],m=100,B=200) 
  #f3<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[3]) 
  #cal5<-calibrate(f3, cmethod="KM", method="boot",u=cut.time[3],m=100,B=200) 
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
    #cal3=get_best_calibrate(cox2,cut.time[2])
    #cal5=get_best_calibrate(cox2,cut.time[3])
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  #customLayout::lay_show(cl)
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  #lines(cal_all[[1]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[1], pch = 16)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
      #lines(cal_all[[i]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[i], pch = 16)
    }
    #plot(cal3,lwd = 2,lty = 0,errbar.col = mg_colors[3],xlim = c(0,1),ylim= c(0,1),col = mg_colors[3],add = T)
    #lines(cal3[,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[3], pch = 16)
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=as.data.frame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  #plot(cal1,xlim=c(0,1),ylim=c(0,1))
  #plot(cal3,xlim=c(0,1),ylim=c(0,1))
  #plot(cal5,xlim=c(0,1),ylim=c(0,1))
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2024)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Partial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}


tcga.cli1<-read.delim('origin_datas/TCGA/Merge_LUAD_clinical.txt',sep='\t',header = T)
colnames(tcga.cli1)[1:20]
tcga.cli1$number_pack_years_smoked
table(tcga.cli1$tobacco_smoking_history)
table(tcga.cli1$diagnosis)
tcga.cli1=data.frame(Samples=tcga.cli1$A0_Samples,
                     Age=tcga.cli1$A17_Age,
                     Gender=tcga.cli1$A18_Sex,
                     Smoking_history=tcga.cli1$tobacco_smoking_history,
                     T.stage=tcga.cli1$A3_T,
                     N.stage=tcga.cli1$A4_N,
                     M.stage=tcga.cli1$A5_M,
                     Stage=tcga.cli1$A6_Stage)
tcga.cli1$Samples=paste0(tcga.cli1$Samples,'-01')
rownames(tcga.cli1)=tcga.cli1$Samples
head(tcga.cli1)
table(tcga.cli1$Smoking_history)

table(tcga.cli1$T.stage)
tcga.cli1$T.stage=gsub('[ab]','',tcga.cli1$T.stage)
tcga.cli1$T.stage[tcga.cli1$T.stage=='TX']<-NA

table(tcga.cli1$N.stage)
tcga.cli1$N.stage[tcga.cli1$N.stage=='NX'|tcga.cli1$N.stage=='']<-NA

table(tcga.cli1$M.stage)
tcga.cli1$M.stage=gsub('[ab]','',tcga.cli1$M.stage)
tcga.cli1$M.stage[tcga.cli1$M.stage=='MX'|tcga.cli1$M.stage=='']<-NA

table(tcga.cli1$Stage)
tcga.cli1$Stage=gsub('[AB]','',tcga.cli1$Stage)
tcga.cli1$Stage[tcga.cli1$Stage=='']<-NA
tcga.cli1$Stage=gsub('Stage ','',tcga.cli1$Stage)


tcga.pancancer.cli=read.xlsx
head(tcga.pancancer.cli)
tcga.cli2=tcga.pancancer.cli[which(tcga.pancancer.cli$type=='LUAD'),]
head(tcga.cli2)
tcga.cli2=data.frame(Samples=paste0(tcga.cli2$bcr_patient_barcode,'-01'),
                     tcga.cli2[,c('OS','OS.time','DSS','DSS.time','DFI','DFI.time','PFI','PFI.time')])
head(tcga.cli2)
tcga.cli2$OS.time
tcga.cli2=tcga.cli2 %>% drop_na(OS.time)
tcga.cli2=tcga.cli2[tcga.cli2$OS.time>0,]
dim(tcga.cli2)

tcga.cli=merge(tcga.cli1,tcga.cli2,by='Samples')
rownames(tcga.cli)=tcga.cli$Samples
tcga.cli=as.data.frame(tcga.cli)
fivenum(as.numeric(tcga.cli$Age))
tcga.cli$Age1=ifelse(as.numeric(tcga.cli$Age)>66,'>66','<=66')
dim(tcga.cli)
head(tcga.cli)
#509


tcga.data<-read.delim('origin_datas/TCGA/LUAD_TPM.txt',sep='\t',header = T,row.names = 1,check.names = F)
tcga.data[1:4,1:4]
table(substr(colnames(tcga.data),14,15))


sample_T=colnames(tcga.data)[which(as.numeric(substr(colnames(tcga.data),14,15))==1)]#肿瘤样本
sample_N=colnames(tcga.data)[which(as.numeric(substr(colnames(tcga.data),14,15))==11)]#正常样本
length(sample_N)
tcga_type=data.frame(Samples=c(sample_T,sample_N),type=rep(c('Tumor','Normal'),c(length(sample_T),length(sample_N))))
rownames(tcga_type)=tcga_type$Samples
table(tcga_type$type)
# Normal  Tumor 
# 59    513

genecode=read.delim
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]


range(tcga.data)
tcga.data=log2(tcga.data[intersect(rownames(tcga.data),mrna_genecode$SYMBOL),tcga_type$Samples]+1)
range(tcga.data)
tcga.exp=tcga.data[,intersect(tcga.cli$Samples,sample_T)]
dim(tcga.exp)
# 19503   500
tcga.cli=tcga.cli[intersect(tcga.cli$Samples,sample_T),]
dim(tcga.cli)

##GSE31210#############
load('origin_datas/GEO/GSE31210.RData')
GSE31210.cli=GSE31210$Sample
GSE31210.cli=data.frame(Samples=GSE31210.cli$Acc,
                        Age=GSE31210.cli$`age (years)`,
                        Gender=GSE31210.cli$gender,
                        Status=GSE31210.cli$death,
                        OS.time=GSE31210.cli$`days before death/censor`)
rownames(GSE31210.cli)=GSE31210.cli$Samples
table(GSE31210.cli$Status)
GSE31210.cli=GSE31210.cli[which(GSE31210.cli$Status!='NULL'),]
GSE31210.cli$OS=ifelse(GSE31210.cli$Status=='alive',0,1)
GSE31210.cli$OS.time
head(GSE31210.cli)

load('origin_datas/GEO/GSE31210_exp.RData')
GSE31210.exp=GSE31210_exp
range(GSE31210.exp)
GSE31210.exp=log2(GSE31210.exp+1)
range(GSE31210.exp)

############
dir.create('results/01.ADME.subtype')
##########
tcga.limma=mg_limma_DEG(exp = tcga.data[rownames(tcga.data)%in%mrna_genecode$SYMBOL,tcga_type$Samples],
                        group = tcga_type$type,ulab = "Tumor",dlab = "Normal")
tcga.limma$Summary
tcga.degs=tcga.limma$DEG[tcga.limma$DEG$adj.P.Val<0.05 & abs(tcga.limma$DEG$logFC)>1,]
write.csv(tcga.degs,'results/01.ADME.subtype/TCGA_DEGs.csv')
dim(tcga.degs)
#2601
my_volcano(dat = tcga.limma,p_cutoff = 0.05,fc_cutoff = 1,col=c("orange","blue","black"))
ggsave('results/01.ADME.subtype/TCGA_volcano.pdf',height = 5,width = 6)

#########
ADME.related.genes=read.xlsx('origin_datas/PDIM_38212774_ADMErelatedgenes.xlsx')
ADME.related.genes=ADME.related.genes$Gene.Symbol
length(ADME.related.genes)
#298

ADME.cox=cox_batch(dat = tcga.exp[ADME.related.genes,tcga.cli$Samples],time = tcga.cli$OS.time,event = tcga.cli$OS)
ADME.cox=na.omit(ADME.cox)
ADME.cox
table(ADME.cox$p.value<0.05)
write.csv(ADME.cox,'results/01.ADME.subtype/ADME_cox_genes.csv')

DE.ADME.genes=intersect(rownames(tcga.degs),rownames(ADME.cox[ADME.cox$p.value<0.05,]))
length(DE.ADME.genes)
#22
library(eulerr)
v=list(rownames(tcga.degs),rownames(ADME.cox[ADME.cox$p.value<0.05,]))
names(v)=c('TCGA DEGs','ADME.cox')
venn.plot=plot(venn(v),labels = list(col = "gray20", font = 2), 
               edges = list(col="gray60", lex=1),
               fills = list(fill = c("#FF707FFF", "#6B58EEFF"), alpha = 0.6),
               quantities = list(cex=.8, col='gray20'))
venn.plot
ggsave('results/01.ADME.subtype/venn.plot.pdf',venn.plot,height = 5,width = 9)

pdf('results/01.ADME.subtype/forest_plot.pdf',height = 5,width = 7)
bioForest(rt = ADME.cox[DE.ADME.genes,],col=c('blue','red'))
dev.off()

############
library(ConsensusClusterPlus)
clusterAlg_name=c('hc','pam','km','kmdist')[3]
distance_name=c('pearson','spearman','euclidean','binary','maximum','canberra','minkowski')[3]
consen_gene=DE.ADME.genes
length(consen_gene)
tcga_consen_data=as.matrix(tcga.exp[consen_gene,tcga.cli$Samples])
# tcga_consen_data=t(scale(t(tcga_consen_data),scale = T))   #11,2  21,2  33,23
tcga_consen_data=t(scale(t(tcga_consen_data),scale = F))   #11,3   21,23
#tcga_consen_data=sweep(tcga_consen_data,1,apply(tcga_consen_data, 1, mean))#
# tcga_consen_data=sweep(tcga_consen_data,1,apply(tcga_consen_data, 1, median))  #
#tcga_consen_data=as.dist(1-cor(tcga_consen_data,method = 'pearson'))
tcga_consen_data=as.matrix(tcga_consen_data)
dim(tcga_consen_data)
tcga_clust_subtype <- ConsensusClusterPlus(tcga_consen_data
                                           , maxK = 10, reps = 500, pItem = 0.8
                                           , pFeature = 1
                                           , title = "TCGA_subtype"
                                           , clusterAlg = clusterAlg_name
                                           , distance = distance_name
                                           , plot = "pdf"
                                           , writeTable = T
                                           , seed = 123456)

k=2
cluster.color=pal_nejm()(8)[c(3:6)]
tcga.subtype <- data.frame(Samples = names(tcga_clust_subtype[[k]]$consensusClass),
                           Cluster=tcga_clust_subtype[[k]]$consensusClass)
tcga.subtype$Cluster=paste0('C',tcga.subtype$Cluster)
table(tcga.subtype$Cluster)
#write.csv(tcga.subtype[order(tcga.subtype$Cluster),],'results/01.cluster/TCGA_subtype.csv',row.names = F)
tcga.subtype.cli=merge(tcga.subtype,tcga.cli,by='Samples')
rownames(tcga.subtype.cli)=tcga.subtype.cli$Samples


tcga.subtype.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Cluster,data = tcga.subtype.cli),
                           data=tcga.subtype.cli,
                           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,surv.median.line = 'hv',
                           title='TCGA-LUAD',ggtheme=custom_theme(),
                           linetype = c("solid", "dashed","strata")[1],
                           palette = cluster.color,
                           legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                           # legend = c(0.8,0.75), 
                           legend.title = "")
tcga.subtype.km=mg_merge_plot(tcga.subtype.km$plot,tcga.subtype.km$table,
                              ncol = 1,nrow = 2,heights = c(2.5,1),align = 'v')
tcga.subtype.km
ggsave('results/01.ADME.subtype/TCGA_KM.pdf',tcga.subtype.km,width = 5,height = 5)

library(ggbiplot)
tcga.subtype.pca <- prcomp(t(tcga.exp[consen_gene,tcga.subtype.cli$Samples]), scale=T)
ggbiplot(tcga.subtype.pca, scale=1, groups = tcga.subtype.cli$Cluster,
                       ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_manual(values = cluster.color) + 
  theme_bw() +
  theme(legend.direction = 'horizontal', legend.position = 'top', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab('PCA1') + ylab('PCA2')+ggtitle('TCGA')
ggsave('results/01.ADME.subtype/TCGA_PCA.pdf',width = 5,height = 5)


GSE31210_consen_data=as.matrix(GSE31210.exp[consen_gene,GSE31210.cli$Samples])
# GSE31210_consen_data=t(scale(t(GSE31210_consen_data),scale = T))   #11,2  21,2  33,23
GSE31210_consen_data=t(scale(t(GSE31210_consen_data),scale = F))   #11,3   21,23
#GSE31210_consen_data=sweep(GSE31210_consen_data,1,apply(GSE31210_consen_data, 1, mean))#
# GSE31210_consen_data=sweep(GSE31210_consen_data,1,apply(GSE31210_consen_data, 1, median))  #
#GSE31210_consen_data=as.dist(1-cor(GSE31210_consen_data,method = 'pearson'))
GSE31210_consen_data=as.matrix(GSE31210_consen_data)
dim(GSE31210_consen_data)
GSE31210_clust_subtype <- ConsensusClusterPlus(GSE31210_consen_data
                                               , maxK = 10, reps = 500, pItem = 0.8
                                               , pFeature = 1
                                               , title = "GSE31210_subtype"
                                               , clusterAlg = clusterAlg_name
                                               , distance = distance_name
                                               , plot = "pdf"
                                               , writeTable = T
                                               , seed = 123456)

GSE31210.subtype <- data.frame(Samples = names(GSE31210_clust_subtype[[k]]$consensusClass),
                               Cluster=GSE31210_clust_subtype[[k]]$consensusClass)
GSE31210.subtype$Cluster=paste0('C',GSE31210.subtype$Cluster)
GSE31210.subtype$Cluster=gsub('C1','S2',GSE31210.subtype$Cluster)
GSE31210.subtype$Cluster=gsub('C2','S1',GSE31210.subtype$Cluster)
GSE31210.subtype$Cluster=gsub('S','C',GSE31210.subtype$Cluster)
table(GSE31210.subtype$Cluster)
GSE31210.subtype.cli=merge(GSE31210.subtype,GSE31210.cli,by='Samples')
rownames(GSE31210.subtype.cli)=GSE31210.subtype.cli$Samples


GSE31210.subtype.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Cluster,data = GSE31210.subtype.cli),
                               data=GSE31210.subtype.cli,
                               conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,surv.median.line = 'hv',
                               title='GSE31210',ggtheme=custom_theme(),
                               linetype = c("solid", "dashed","strata")[1],
                               palette = cluster.color,
                               legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                               # legend = c(0.8,0.75), 
                               legend.title = "")
GSE31210.subtype.km=mg_merge_plot(GSE31210.subtype.km$plot,GSE31210.subtype.km$table,
                                  ncol = 1,nrow = 2,heights = c(2.5,1),align = 'v')
GSE31210.subtype.km
ggsave('results/01.ADME.subtype/GSE31210_KM.pdf',GSE31210.subtype.km,width = 5,height = 5)

library(ggbiplot)
GSE31210.subtype.pca <- prcomp(t(GSE31210.exp[consen_gene,GSE31210.subtype.cli$Samples]), scale=T)
ggbiplot(GSE31210.subtype.pca, scale=1, groups = GSE31210.subtype.cli$Cluster,
         ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_manual(values =cluster.color) + 
  theme_light() +
  theme(legend.direction = 'horizontal', legend.position = 'top', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab('PCA1') + ylab('PCA2')+ggtitle('GSE31210')
ggsave('results/01.ADME.subtype/GSE31210_PCA.pdf',width = 5,height = 5)

#######
head(tcga.subtype.cli)
cluster.anno=data.frame(Cluster=tcga.subtype.cli$Cluster[order(tcga.subtype.cli$Cluster)])
rownames(cluster.anno)=tcga.subtype.cli$Samples[order(tcga.subtype.cli$Cluster)]

cluster.color.use=cluster.color[1:2]
names(cluster.color.use)=c('C1','C2')

pdf('results/01.ADME.subtype/tcga_gene_expr.pdf',height = 5,width = 9,onefile = F)
pheatmap(tcga.exp[consen_gene,rownames(cluster.anno)], 
         scale='row',name = 'Expression',
         annotation_col = cluster.anno,
         annotation_colors = list(Cluster=cluster.color.use),
         color=colorRampPalette(c('#1A5592','snow',"orange"))(100),
         cluster_cols = F, cluster_rows = T,show_colnames = F,
         fontsize_row = 10,fontsize_col = 12)
dev.off()



head(GSE31210.subtype.cli)
cluster.anno=data.frame(Cluster=GSE31210.subtype.cli$Cluster[order(GSE31210.subtype.cli$Cluster)])
rownames(cluster.anno)=GSE31210.subtype.cli$Samples[order(GSE31210.subtype.cli$Cluster)]
pdf('results/01.ADME.subtype/GSE31210_gene_expr.pdf',height = 5,width =9,onefile = F)
pheatmap(GSE31210.exp[consen_gene,rownames(cluster.anno)], 
         scale='row',name = 'Expression',
         annotation_col = cluster.anno,
         annotation_colors = list(Cluster=cluster.color.use),
         color=colorRampPalette(c('#1A5592','snow',"orange"))(100),
         cluster_cols = F, cluster_rows = T,show_colnames = F,
         fontsize_row = 10,fontsize_col = 12)
dev.off()


#########
dir.create('results/02.ImmuneCell.infiltration')
########


#########
tcga.immune.ssgsea=immu_ssgsea(exp = tcga.exp)
tcga.immune.ssgsea[1:5,1:5]
tme.df2=tcga.immune.ssgsea[tcga.subtype.cli$Samples,]
tme.df2=as.data.frame(tme.df2)
tme.df2$Cluster=tcga.subtype.cli$Cluster
tme.df2=melt(tme.df2)
head(tme.df2)

fig2a=tme.df2 %>%
  ggplot(aes(x=variable,y=value,fill=Cluster))+
  geom_boxplot()+stat_compare_means(aes(group=Cluster), label = "p.signif", method = 'wilcox.test')+
  scale_fill_manual(values =cluster.color)+
  xlab('')+ylab('Score')+ggtitle('28 immune cell infiltration')+
  theme_bw()+theme(text = element_text(family = 'Times',size = 12),legend.position = 'top',
                   axis.text.x = element_text(color = "black", size = 12,angle = 45,hjust = 1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
fig2a

tcga.est=read.delim('results/02.ImmuneCell.infiltration/TCGA_ESTIMATE_score.txt',row.names = 1)
head(tcga.est)

tme.df1=tcga.est
tme.df1$Cluster=tcga.subtype.cli$Cluster
tme.df1=melt(tme.df1)
head(tme.df1)
table(tme.df1$variable)
fig2b=tme.df1 %>%subset(variable!='TumorPurity')%>%
  ggplot(aes(x=variable,y=value,fill=Cluster))+
  geom_boxplot()+stat_compare_means(aes(group=Cluster), label = "p.format", method = 'wilcox.test')+
  scale_fill_manual(values =cluster.color)+
  xlab('')+ylab('Score')+ggtitle('ESTIMATE')+
  theme_bw()+theme(text = element_text(family = 'Times',size = 12),legend.position = 'top',
                   # axis.text.x = element_text(color = "black", size = 12,angle = 30,hjust = 1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
fig2b

wb_beeswarm_plot <- function(dat = NULL,
                             show_compare = T,
                             xlab = 'Groups',
                             ylab = '',
                             method = c('t.test', 'wilcox.test')[1],
                             col = risktype.col,
                             leg.pos = c('top','left','right','bottom','none')[1],
                             title = NULL,
                             group = 'Cluster') {
  library(ggbeeswarm)
  colnames(dat) <- c('Cluster', 'Feature')
  
  
  p1 <- ggplot(dat, aes(Cluster, Feature, color = Cluster)) + geom_quasirandom(method = "frowney") +
    ggtitle(title) + scale_color_manual(values = col[1:length(unique(dat$Cluster))]) +
    xlab(xlab) + ylab(ylab) + guides(color=guide_legend(title = group)) + theme_classic() +
    theme(legend.position=leg.pos,text = element_text(family = 'Times',size = 14))
  
  
  if(show_compare){
    uni.group = as.character(unique(dat$Cluster))
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,
                                     method = method,
                                     label= "p.signif",
                                     step_increase = 0.0)
  }
  return(p1)
}

test.df=tme.df1 %>% subset(variable=='TumorPurity')
test.df=test.df[,-2]
head(test.df)
fig2c=wb_beeswarm_plot(dat =test.df,col = cluster.color,ylab = 'TumorPurity')
fig2c

##GSEA##########
tcga.geneList=getGeneFC(gene.exp=tcga.exp[,tcga.subtype.cli$Samples],group=tcga.subtype.cli$Cluster
                        ,ulab='C1',dlab = 'C2')

gmt2list <- function(annofile){
  
  if (!file.exists(annofile)) {
    stop("There is no such a gmt file!")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
  return(annoList)
}
kegmt=gmt2list
library(fgsea)
set.seed(111)
fgseaRes <- fgsea(pathways = kegmt,
                  stats =tcga.geneList ,
                  minSize=10,
                  maxSize=500,
                  nperm=1000)
head(fgseaRes)

nrow(fgseaRes[padj<0.05&NES > 0])
#11
nrow(fgseaRes[padj<0.05&NES < 0])
#3


cluster.gsea.res=fgseaRes[fgseaRes$padj<0.05,]
head(cluster.gsea.res)
write.xlsx(cluster.gsea.res,'results/02.ImmuneCell.infiltration/cluster_GSEA.xlsx',overwrite = T)
table(cluster.gsea.res$NES>0)
cluster.gsea.res$group=ifelse(cluster.gsea.res$NES>0,"Activated","Suppressed")
cluster.gsea.res$pathway=gsub('HALLMARK_','',cluster.gsea.res$pathway)
cluster.gsea.res$pathway=gsub('_',' ',cluster.gsea.res$pathway)
cluster.gsea.res$pathway <- factor(cluster.gsea.res$pathway,
                                       levels = cluster.gsea.res$pathway[order(cluster.gsea.res$NES)])


fig2d=ggplot(cluster.gsea.res, aes(x = pathway, y = NES, fill =group)) +
  geom_bar(stat = "identity") +
  coord_flip() +  
  scale_fill_manual(values = c("#FDB462", "#80B1D3")) +
  labs(x = "", y = "NES") +
  theme_classic()+ theme(text = element_text(family = 'Times',size = 14,face = 'bold'),
                         legend.position = 'top')
fig2d

fig2=mg_merge_plot(fig2a,mg_merge_plot(fig2b,fig2c,fig2d,ncol=3,widths = c(1.2,1,1.5),labels = LETTERS[2:4]),
                   nrow=2,labels = c('A',''))
ggsave('results/02.ImmuneCell.infiltration/Fig2.pdf',fig2,height = 12,width = 16)

#########
dir.create('results/03.Subtype.DEGs')
tcga.subtype.limma=mg_limma_DEG(exp=tcga.exp[,tcga.subtype.cli$Samples],group=tcga.subtype.cli$Cluster
                               ,ulab='C1',dlab = 'C2')
tcga.subtype.limma$Summary
tcga.subtype.degs=tcga.subtype.limma$DEG[tcga.subtype.limma$DEG$adj.P.Val<0.05 & abs(tcga.subtype.limma$DEG$logFC)>1,]
fig3a=my_volcano(dat = tcga.subtype.limma,p_cutoff = 0.05,fc_cutoff = 1,col=c("orange","blue","black"))+ggtitle('C1 vs C2')
fig3a
write.csv(tcga.subtype.degs,'results/03.Subtype.DEGs/tcga.subtype.degs.csv')

tcga.DEGs.filter=tcga.subtype.degs
tcga.DEGs.filter$type=ifelse(tcga.DEGs.filter$logFC>0,'up','down')
tcga.DEGs.filter <- tcga.DEGs.filter %>%
  mutate(gene = rownames(tcga.DEGs.filter))  %>%
  group_by(type) %>% slice_max(n =50, order_by = abs(logFC)) %>%  arrange(logFC)
head(tcga.DEGs.filter)

cli_anno=tcga_type[,'type',drop=F]
fig3b=pheatmap(tcga.data[tcga.DEGs.filter$gene,rownames(cli_anno)],
               scale = 'row',name = 'Expression',  main="TOP50 DEGs in TCGA", 
               color =  circlize::colorRamp2(c(-3, 0, 3), c('#009B9F', "white", '#C75DAA')),
               annotation_col = cli_anno,
               annotation_colors = list(type=c(Tumor='#FB8072',Normal='#80B1D3')),
               cluster_cols = F, 
               cluster_rows = T,
               show_rownames = F, 
               show_colnames = F)
library(ggplotify)
fig3b = as.ggplot(fig3b)
fig3b

up.DEGs.enrich=mg_clusterProfiler(genes =rownames(tcga.subtype.degs[tcga.subtype.degs$logFC>0,]))
enrichplot::dotplot(up.DEGs.enrich$KEGG)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))
fig3c=enrichplot::dotplot(up.DEGs.enrich$GO_BP,showCategory = 8)+ggtitle('Biological Process')+
  theme(text=element_text(family = 'Times'))
fig3c
write.xlsx(up.DEGs.enrich$GO_BP@result,'results/03.Subtype.DEGs/up.DEGs.enrich.xlsx')

dn.DEGs.enrich=mg_clusterProfiler(genes =rownames(tcga.subtype.degs[tcga.subtype.degs$logFC<0,]))
enrichplot::dotplot(dn.DEGs.enrich$KEGG)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))
fig3d=enrichplot::dotplot(dn.DEGs.enrich$GO_BP)+ggtitle('Biological Process')+
  theme(text=element_text(family = 'Times'))
fig3d
write.xlsx(dn.DEGs.enrich$GO_BP@result,'results/03.Subtype.DEGs/dn.DEGs.enrich.xlsx')

pdf('results/03.Subtype.DEGs/Fig3.pdf',height = 10,width = 12)
mg_merge_plot(fig3a,fig3b,fig3c,fig3d,labels = LETTERS[1:4],heights = c(1,1.2))
dev.off()

#############
dir.create('results/04.prognosis_model')
tcga.subtype.degs.cox=cox_batch(dat = tcga.exp[rownames(tcga.subtype.degs),tcga.cli$Samples],
                                time = tcga.cli$OS.time,event = tcga.cli$OS)
table(tcga.subtype.degs.cox$p.value<0.001)
tcga.subtype.degs.cox.fit=tcga.subtype.degs.cox[tcga.subtype.degs.cox$p.value<0.001,]

pre.genes=rownames(tcga.subtype.degs.cox.fit)
length(pre.genes)
#127

tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[pre.genes, tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))

#####
library(glmnet)
tcga.lasso=get_riskscore.lasso(dat = tcga_model_data[,-c(1:2)],
                               os = tcga_model_data$OS,
                               os.time = tcga_model_data$OS.time )
length(tcga.lasso$lasso.gene)
tcga.lasso$plot

# tcga.lasso$lasso.gene

######
fmla <- as.formula(paste0("Surv(OS.time, OS) ~"
                          ,paste0(tcga.lasso$lasso.gene,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
cox=step(cox)

lan <- coef(cox)
lan
length(lan)
paste0(round(lan, 3), '*', names(lan),collapse = '+')
# "0.143*ANLN+0.111*FAM83A+-0.105*IRX5+0.096*RHOV"

gene.coef=data.frame(gene=names(lan),coef=as.numeric(lan))
gene.coef
gene.coef$Type=ifelse(gene.coef$coef>0,'Risk','Protective')
gene.coef.fig=ggplot(gene.coef, aes(x = coef, y = reorder(gene,coef), fill =Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#FDB462", "#6154AB")) +
  labs(x = 'coefficient', y = "") +
  geom_text(aes(label = round(coef,2),hjust =2), data = subset(gene.coef, coef > 0))+ 
  geom_text(aes(label = round(coef,2), hjust = -1), data = subset(gene.coef, coef < 0))+  
  theme_bw()+ theme(text = element_text(family = 'Times',size = 14,face = 'bold'),legend.position = 'top')
gene.coef.fig


gene.forest=ggforest(cox, data = tcga_model_data, 
                     main = "Hazardratio", fontsize =1.0, 
                     noDigits = 2)
gene.forest




risktype.col=c("pink", "#009B9F")
###########
risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.cli,Riskscore=risk.tcga)
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>median(risk.tcga),'High','Low')


tcga.roc=ggplotTimeROC(time = tcga.risktype.cli$OS.time,
                       status = tcga.risktype.cli$OS,
                       score = tcga.risktype.cli$Riskscore,mks = c(1:5))
tcga.roc

tcga.km.OS=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                   data = tcga.risktype.cli),
                      data=tcga.risktype.cli,
                      conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                      surv.median.line = 'hv',title='TCGA',
                      linetype = c("solid", "dashed","strata")[1],
                      palette = risktype.col,ggtheme = custom_theme(),
                      legend = c(0.8,0.85), 
                      legend.title = "Risktype",legend.labs=c('High','Low'))
tcga.km.OS=mg_merge_plot(tcga.km.OS$plot,tcga.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
tcga.km.OS


tcga.km.PFI=ggsurvplot(fit=survfit( Surv(PFI.time/365, PFI) ~ Risktype,
                                   data = tcga.risktype.cli),
                      data=tcga.risktype.cli,
                      conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                      surv.median.line = 'hv',title='TCGA',
                      linetype = c("solid", "dashed","strata")[1],
                      palette = risktype.col,ggtheme = custom_theme(),
                      legend = c(0.8,0.85), 
                      legend.title = "Risktype",legend.labs=c('High','Low'))
tcga.km.PFI=mg_merge_plot(tcga.km.PFI$plot,tcga.km.PFI$table,nrow=2,heights = c(2.5,1),align = 'v')
tcga.km.PFI

data_ex = cbind.data.frame(t(tcga.exp[names(lan),tcga.risktype.cli$Samples]),Risktype=tcga.risktype.cli$Risktype)
data_ex=data_ex[order(data_ex$Risktype),]
data_ex
library(ggpubr)
p=list()
for(i in 1:length(names(lan))){
  dt<-data_ex[,c("Risktype",names(lan)[i])]
  p[[i]]<-ggboxplot(dt,x='Risktype',y=colnames(dt)[2],fill = "Risktype", 
                    palette =risktype.col,main=names(lan)[i],ylab="Expression",xlab="")+
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none",text = element_text(family = 'Times'))+
    stat_compare_means(method = "wilcox.test",label.x=1.5,label="p.signif")
  assign(paste("p", i, sep=""), p)
}

tcga.expr.fig=mg_merge_plot(p,ncol=4,nrow=1)
tcga.expr.fig



GSE31210_model_data <- cbind(GSE31210.cli[, c("OS.time", "OS")],
                             t(GSE31210.exp[pre.genes, GSE31210.cli$Samples]))
colnames(GSE31210_model_data) <- gsub('-', '_', colnames(GSE31210_model_data))

risk.GSE31210=as.numeric(lan%*%as.matrix(t(GSE31210_model_data[GSE31210.cli$Samples,names(lan)])))
GSE31210.risktype.cli=data.frame(GSE31210.cli,Riskscore=risk.GSE31210)
GSE31210.risktype.cli$Risktype=ifelse(GSE31210.risktype.cli$Riskscore>median(risk.GSE31210),'High','Low')

GSE31210.roc=ggplotTimeROC(GSE31210.risktype.cli$OS.time,
                           GSE31210.risktype.cli$OS,
                           GSE31210.risktype.cli$Riskscore,mks = c(1:5))
GSE31210.roc



GSE31210.km.OS=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                       data = GSE31210.risktype.cli),
                          data=GSE31210.risktype.cli,
                          conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                          surv.median.line = 'hv',title='GSE31210',
                          linetype = c("solid", "dashed","strata")[1],
                          palette = risktype.col,ggtheme = custom_theme(),
                          legend = 'top', 
                          legend.title = "Risktype",legend.labs=c('High','Low'))
GSE31210.km.OS=mg_merge_plot(GSE31210.km.OS$plot,GSE31210.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
GSE31210.km.OS

data_ex2 = cbind.data.frame(t(GSE31210.exp[names(lan),GSE31210.risktype.cli$Samples]),
                            Risktype=GSE31210.risktype.cli$Risktype)
data_ex2=data_ex2[order(data_ex2$Risktype),]
data_ex2
library(ggpubr)
p=list()
for(i in 1:length(names(lan))){
  dt<-data_ex2[,c("Risktype",names(lan)[i])]
  p[[i]]<-ggboxplot(dt,x='Risktype',y=colnames(dt)[2],fill = "Risktype", 
                    palette =risktype.col,main=names(lan)[i],ylab="Expression",xlab="")+
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none",text = element_text(family = 'Times'))+
    stat_compare_means(method = "wilcox.test",label.x=1.5,label="p.signif")
  assign(paste("p", i, sep=""), p)
}

GSE31210.expr.fig=mg_merge_plot(p,ncol=4,nrow=1)
GSE31210.expr.fig

fig4=mg_merge_plot(mg_merge_plot(tcga.lasso$plot,gene.coef.fig,gene.forest,ncol=3,widths = c(2,1,1.2),labels = c('','C','D')),
                   mg_merge_plot(tcga.roc,tcga.km.OS,GSE31210.roc,GSE31210.km.OS,ncol=4,labels = LETTERS[5:8]),
                   mg_merge_plot(tcga.expr.fig,GSE31210.expr.fig,labels = LETTERS[9:10]),
                   nrow=3)
ggsave('results/04.prognosis_model/Fig4.pdf',fig4,height = 15,width = 20)

###########
dir.create('results/05.nomogram')
################
head(tcga.risktype.cli)
tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)
table(tcga_cox_datas$Age1)

table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='I'|tcga_cox_datas$Stage=='II']<-'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='III'|tcga_cox_datas$Stage=='IV']<-'III+IV'

table(tcga_cox_datas$T.stage)
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T1'|tcga_cox_datas$T.stage=='T2']<-'T1+T2'
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T3'|tcga_cox_datas$T.stage=='T4']<-'T3+T4'

table(tcga_cox_datas$N.stage)
tcga_cox_datas$N.stage[tcga_cox_datas$N.stage=='N1'|tcga_cox_datas$N.stage=='N2'|tcga_cox_datas$N.stage=='N3']<-'N1+N2+N3'

table(tcga_cox_datas$M.stage)


#Age
tcga_cox_datas=crbind2DataFrame(tcga_cox_datas)
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age1,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat

##Sex
Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                             data=tcga_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                              HR = round(Gender_sig_cox[[7]][,2],3),
                              lower.95 = round(Gender_sig_cox[[8]][,3],3),
                              upper.95 = round(Gender_sig_cox[[8]][,4],3),
                              p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat


#stage
Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat

#t_stage
t_stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~T.stage,
                                 data=tcga_cox_datas))
t_stage_sig_cox_dat <- data.frame(Names=rownames(t_stage_sig_cox[[8]]),
                                  HR = round(t_stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(t_stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(t_stage_sig_cox[[8]][,4],3),
                                  p.value=round(t_stage_sig_cox[[7]][,5],3))
t_stage_sig_cox_dat

#n_stage
n_stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~N.stage,
                                 data=tcga_cox_datas))
n_stage_sig_cox_dat <- data.frame(Names=rownames(n_stage_sig_cox[[8]]),
                                  HR = round(n_stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(n_stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(n_stage_sig_cox[[8]][,4],3),
                                  p.value=round(n_stage_sig_cox[[7]][,5],3))
n_stage_sig_cox_dat

#m_stage
m_stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~M.stage,
                                 data=tcga_cox_datas))
m_stage_sig_cox_dat <- data.frame(Names=rownames(m_stage_sig_cox[[8]]),
                                  HR = round(m_stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(m_stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(m_stage_sig_cox[[8]][,4],3),
                                  p.value=round(m_stage_sig_cox[[7]][,5],3))
m_stage_sig_cox_dat


#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=tcga_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     Stage_sig_cox_dat,
                     t_stage_sig_cox_dat,
                     n_stage_sig_cox_dat,
                     m_stage_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Features=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
rownames(data.sig) <-c('Age','Gender','Stage','T.stage','N.stage','M.stage','Riskscore')
data.sig$Features=rownames(data.sig)
data.sig
data.sig$p.value=ifelse(data.sig$p.value<0.001,'<0.001',data.sig$p.value)
pdf('results/05.nomogram/Univariate.pdf',height = 4,width = 6,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap =5,lineheight = 9,zero = 1,
                 boxsize = .2,lwd.zero=1,lwd.ci=1.5,lwd.xaxis=1,
                 box_col='#4169E1',summary_col="black",lines_col='black',zero_col='grey',
                 xlab='Hazard Ratio',lty.ci = 6,graph.pos =4)
dev.off()



#########
#T.stage +  N.stage +M.stage
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage+T.stage +  N.stage +M.stage+Riskscore, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Features=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
data.muti
rownames(data.muti) <- c("Stage",'T.stage','N.stage','M.stage',"RiskScore")
data.muti$Features=rownames(data.muti)
data.muti
data.muti$p.value=ifelse(data.muti$p.value<0.001,'<0.001',data.muti$p.value)
pdf('results/05.nomogram/Multivariate.pdf',height = 4,width =6,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap =5,lineheight = 9,zero = 1,
                 boxsize = .2,lwd.zero=1,lwd.ci=1.5,lwd.xaxis=1,
                 box_col='#4169E1',summary_col="black",lines_col='black',zero_col='grey',
                 xlab='Hazard Ratio',lty.ci = 6,graph.pos =4)
dev.off()


#########

pdf('results/05.nomogram/nomogram.pdf', width = 12, height = 10)
nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$Riskscore,
                                T.stage=tcga_cox_datas$T.stage,
                                N.stage=tcga_cox_datas$N.stage),
                     os = tcga_cox_datas$OS.time,
                     status = tcga_cox_datas$OS,
                     mks = c(1,3,5)
)
dev.off()
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))


sum1=summary(coxph(formula=Surv(OS.time, OS)~T.stage+N.stage+Riskscore, data=tcga_cox_datas))
sum1$concordance[1]
c_index=as.numeric(sum1$concordance[1])

clinical.feature=c('Age','Gender','T.stage','N.stage','M.stage','Stage','Riskscore')
for (i in 1:length(clinical.feature)) {
  # i=2
  as.formula(paste0("Surv(OS.time, OS) ~",clinical.feature[i]))
  sum_res=summary(coxph(as.formula(paste0("Surv(OS.time, OS) ~",clinical.feature[i])), data=tcga_cox_datas))
  index=as.numeric(sum_res$concordance[1])
  c_index=append(c_index,index)
}
names(c_index)=c('Nomogram',clinical.feature)
c_index[order(c_index)]
pdf('results/05.nomogram/c_index_barplot.pdf',height = 5,width = 10,onefile = F)
barplot(c_index[order(c_index)],xlab="Clinical Features", ylab="C-index",ylim = c(0,0.9),main = 'Prediction effect',col = '#FFBFE5FF')
dev.off()


############
dir.create('results/06.TME')


tcga.immune.ssgsea[1:5,1:5]
tcga.immune.ssgsea=as.data.frame(tcga.immune.ssgsea)
fig6a=mg_PlotMutiBoxplot(data = tcga.immune.ssgsea[tcga.risktype.cli$Samples,],group = tcga.risktype.cli$Risktype,
                   group_cols = risktype.col,test_method = 'wilcox.test',legend.pos = 'top',add = 'boxplot')
fig6a

tcga_tide_res<-read.csv('results/06.TME/tcga.tide.res.csv',row.names = 1,stringsAsFactors = F)
head(tcga_tide_res)

tme.TIDE=data.frame(tcga_tide_res[tcga.risktype.cli$Samples,c('TIDE','Responder')],tcga.risktype.cli)
head(tme.TIDE)
fig6b=wb_beeswarm_plot(dat =tme.TIDE[,c('Risktype','TIDE')],col = risktype.col,
                       ylab = 'TIDE',leg.pos = 'none',xlab = '')
fig6b


library(IOBR)
IPS.res=IPS_calculation(project = 'LUAD',eset = tcga.exp,plot = F)
head(IPS.res)
dim(IPS.res)

tme.IPS=data.frame(IPS.res[tcga.risktype.cli$Samples,'IPS',drop=F],tcga.risktype.cli)
head(tme.IPS)
fig6c=tme.IPS %>%
  ggplot(aes(x=Risktype,y=IPS,fill=Risktype))+
  geom_violin()+  
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = NA,fill="white")+
  stat_compare_means(aes(group=Risktype), label = "p.format", method = 'wilcox.test')+
  scale_fill_manual(values =risktype.col)+theme_bw()+ylab('Immune Phenotype Score')+
  theme(text = element_text(family = 'Times',size = 15),legend.position = 'none',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
fig6c

#########
GSE135222.cli=read.delim('origin_datas/GEO/GSE135222_CLINICAL.txt')
head(GSE135222.cli)
colnames(GSE135222.cli)[1]='Samples'
dim(GSE135222.cli)

GSE135222.exp=read_tsv('origin_datas/GEO/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz')
GSE135222.exp[1:5,1:5]
GSE135222.exp=as.data.frame(GSE135222.exp)
GSE135222.exp$gene_id=str_split_fixed(GSE135222.exp$gene_id,'[.]',2)[,1]
rownames(GSE135222.exp)=GSE135222.exp$gene_id
GSE135222.exp=GSE135222.exp[,-1]
range(GSE135222.exp)
GSE135222.exp=log2(GSE135222.exp+1)

com.sample=intersect(colnames(GSE135222.exp),GSE135222.cli$Samples)
GSE135222.cli=GSE135222.cli[GSE135222.cli$Samples %in% com.sample,]
dim(GSE135222.cli)


library(clusterProfiler)
library(org.Hs.eg.db)
#ENSEMBL
keytypes(org.Hs.eg.db)
gene.df <- bitr(names(lan), fromType = "SYMBOL", #fromType
                toType = 'ENSEMBL', #toType
                OrgDb = org.Hs.eg.db)#Orgdb
head(gene.df)

intersect(GSE135222.exp$gene_id,gene.df$ENSEMBL)
GSE135222_model_data <- cbind(GSE135222.cli,
                              t(GSE135222.exp[gene.df$ENSEMBL, GSE135222.cli$Samples]))
rownames(GSE135222_model_data)=GSE135222_model_data$Samples
colnames(GSE135222_model_data)[5:8]=gene.df$SYMBOL


fmla1 <- as.formula(paste0("Surv(PFS, PFS_event) ~"
                           ,paste0(names(lan),collapse = '+')))
cox1 <- coxph(fmla1, data =as.data.frame(GSE135222_model_data))
lan1 <- coef(cox1)
lan1

risk.GSE135222=as.numeric(lan1%*%as.matrix(t(GSE135222_model_data[GSE135222.cli$Samples,names(lan1)])))
GSE135222.risktype.cli=data.frame(GSE135222.cli,Riskscore=risk.GSE135222)
GSE135222.risktype.cli$Risktype=ifelse(GSE135222.risktype.cli$Riskscore>median(risk.GSE135222),'High','Low')


table(GSE135222.risktype.cli$benefit,GSE135222.risktype.cli$Risktype)



fig6d=ggplotTimeROC(GSE135222.risktype.cli$PFS/12,
                    GSE135222.risktype.cli$PFS_event,
                    GSE135222.risktype.cli$Riskscore,mks = c(1,2,3))
fig6d


fig6e=plotMutiBar(table(GSE135222.risktype.cli$benefit,GSE135222.risktype.cli$Risktype),legTitle = 'benefit')
fig6e

fig6f=GSE135222.risktype.cli %>%
  ggplot(aes(x=benefit,y=Riskscore,fill=benefit))+
  geom_violin()+  
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = NA,fill="white")+
  stat_compare_means(aes(group=benefit), label = "p.format", method = 'wilcox.test')+
  scale_fill_manual(values =risktype.col)+theme_bw()+ylab('RiskScore')+ggtitle('GSE135222')+
  theme(text = element_text(family = 'Times',size = 15),legend.position = 'none',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
fig6f




fig6=mg_merge_plot(fig6a,
                   mg_merge_plot(fig6b,fig6c,labels = LETTERS[2:3]),
                   mg_merge_plot(fig6d,fig6e,fig6f,ncol=3,labels = LETTERS[4:6]),
                   nrow=3,labels = c('A','',''))
ggsave('results/06.TME/Fig6.pdf',fig6,height = 15,width = 16)


save.image(file = 'project.RData')
save(tcga.exp,file='results/tcga.exp.RData')
