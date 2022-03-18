data <- read.csv('~/Dropbox/IronMetab/2017_data/MetabMat2.csv')
names(data)
IS.ref <- c(14:23)
ID.pre <- c(24:35)
ID.post <- c(36:47)
data.is7mo = read.csv('~/Dropbox/IronMetab/2017_data/IS_7mo_metab.csv') 

IS.ID.rat <- 1/(rowMeans(data[,IS.ref])/rowMeans(data[,ID.pre]))
IS.ID.lograt <- rowMeans(log2(data[,IS.ref]))-rowMeans(log2(data[,ID.pre]))
IS.ID.pvals <- c()
for(i in 1:dim(data)[1]){
  if(i==356) IS.ID.pvals[i] <- NA
  if(i!=356) IS.ID.pvals[i] <- t.test(log2(data[i,IS.ref]),log2(data[i,ID.pre]))$p.value
}

IS.IDp.rat <- 1/(rowMeans(data[,IS.ref])/rowMeans(data[,ID.post]))
IS.IDp.lograt <- rowMeans(log2(data[,ID.post]))-rowMeans(log2(data[,IS.ref]))
IS.IDp.pvals <- c()
for(i in 1:dim(data)[1]){
  if(i==356) IS.IDp.pvals[i] <- NA
  if(i!=356) IS.IDp.pvals[i] <- t.test(log2(data[i,IS.ref]),log2(data[i,ID.post]))$p.value
}

ID.IDp.rat <- 1/(rowMeans(data[,ID.post])/rowMeans(data[,ID.pre]))
ID.IDp.lograt <- rowMeans(log2(data[,ID.pre]))-rowMeans(log2(data[,ID.post]))
ID.IDp.pvals <- c()
for(i in 1:dim(data)[1]){
  if(i==356) ID.IDp.pvals[i] <- NA
  if(i!=356) ID.IDp.pvals[i] <- t.test(log2(data[i,ID.post]),log2(data[i,ID.pre]), matched=TRUE)$p.value
}

IS.ISp.rat <- 1/(rowMeans(data[,IS.ref])/rowMeans(data.is7mo))
IS.ISp.lograt <- rowMeans(log2(data[,IS.ref]))-rowMeans(log2(data.is7mo))
IS.ISp.pvals <- c()
for(i in 1:dim(data)[1]){
  if(i==356) IS.ISp.pvals[i] <- NA
  if(i!=356) IS.ISp.pvals[i] <- t.test(log(data[i,IS.ref]),log(data.is7mo), matched=FALSE)$p.value
}

IDp.ISp.rat <- 1/(rowMeans(data[,ID.post])/rowMeans(data.is7mo))
IDp.ISp.lograt <- rowMeans(log2(data[,ID.post]))-rowMeans(log2(data.is7mo))
IDp.ISp.pvals <- c()
for(i in 1:dim(data)[1]){
  if(i==356) IDp.ISp.pvals[i] <- NA
  if(i!=356) IDp.ISp.pvals[i] <- t.test(log2(data[i,ID.post]),log2(data.is7mo), matched=FALSE)$p.value
}

ID.ISp.rat <- 1/(rowMeans(data[,ID.pre])/rowMeans(data.is7mo))
ID.ISp.lograt <- rowMeans(log2(data[,ID.pre]))-rowMeans(log2(data.is7mo))
ID.ISp.pvals <- c()
for(i in 1:dim(data)[1]){
  if(i==356) ID.ISp.pvals[i] <- NA
  if(i!=356) ID.ISp.pvals[i] <- t.test(log2(data[i,ID.pre]),log2(data.is7mo), matched=FALSE)$p.value
}



png(file='~/Dropbox/IronMetab/2017_data/pval_hists.png',width=600,height=500)
par(mfrow=c(3,2))
hist(IS.ID.pvals, main="IS 6mo vs. IS 6mo")
abline(h=65,lty=2)
hist(IS.IDp.pvals, main="IS 6mo vs. Post-IDA")
abline(h=65,lty=2)
hist(ID.IDp.pvals, main="ID 6mo vs. Post-IDA")
abline(h=65,lty=2)
hist(IDp.ISp.pvals, main="IS 7mo vs. Post-IDA")
abline(h=65,lty=2)
hist(IS.ISp.pvals, main="IS 6mo vs. IS 7mo")
abline(h=65,lty=2)
dev.off()

IS.ID.fdr=p.adjust(IS.ID.pvals,method='fdr')
IS.IDp.fdr=p.adjust(IS.IDp.pvals,method='fdr')
ID.IDp.fdr=p.adjust(ID.IDp.pvals,method='fdr')

Table <- cbind(data[,1:13],IS.ID.rat,IS.ID.lograt,IS.ID.pvals,IS.ID.fdr,IS.IDp.rat,IS.IDp.lograt,IS.IDp.pvals,IS.IDp.fdr,ID.IDp.rat,ID.IDp.lograt,ID.IDp.pvals,ID.IDp.fdr)

#write.csv(Table, file='~/Dropbox/IronMetab/2017_data/Metab_ComparisonTable.csv')
write.csv(Table, file='~/Dropbox/IronMetab/2017_data/Metab_ComparisonTable_with_FDR.csv')

pval.mat = cbind(IS.ISp.pvals,ID.ISp.pvals,ID.IDp.pvals,IS.IDp.pvals,IS.ID.pvals,IDp.ISp.pvals)
fdr.mat = pval.mat
for(i in 1:dim(pval.mat)[1]){
  fdr.mat[i,]=p.adjust(pval.mat[i,])
}
total.mat = cbind(pval.mat,fdr.mat) 

total.mat=round(total.mat,3)
total.mat[total.mat==0]='<0.001'
Table.fig <- cbind(data[,1:13],total.mat)
write.csv(Table.fig, file='~/Dropbox/IronMetab/2017_data/Metab_Comparison_For_Figs.csv')


data.mat <- as.matrix(log(data[,14:47]))



###CCA analysis
#load hematology
heme2017 <- read.csv('~/Dropbox/IronMetab/cbc_2017.csv',  header = TRUE) 
hememat2017 <- matrix(as.numeric(as.matrix(heme2017[,c(6,12,14,16,18,20,22)])), ncol=7,nrow=25)
y.pred=c()
for(i in 1:6){
  x=hememat2017[,1]
  y=hememat2017[,i+1]
  blah=lm(y~x)
  summary(blah)
  y.pred[i] =predict(blah)[18]
}
hememat2017.impute = hememat2017
hememat2017.impute[18,2:7]=y.pred



#hememat2017=hememat2017[-c(18),]#18 is missing
hememat2017.s = hememat2017.impute
for(i in 1:7){
  hememat2017.s[,i] = (hememat2017.impute[,i]-mean(hememat2017.impute[,i]))/sd(hememat2017.impute[,i])
}

id.heme = as.character(heme2017$ID)
metabnames = read.csv('~/Dropbox/IronMetab/2017_data/metabnames.csv',header=FALSE)
id.met = gsub(" .*","",as.character(metabnames[,1]))
group = c(rep('IS',10),rep('ID',12),rep('ID_post',12))
id.met.filt = id.met[group=='IS'|group=='ID']
data.mat.filt=data.mat[,group=='IS'|group=='ID']
#id.met.filt2=id.met.filt[-c(15)] #no BM02 in heme data
#data.mat.filt = data.mat.filt[,-c(15)]
group.filt=group[-c(23:34)]
datamat.filt.s=data.mat.filt
for(i in 1:654){
  datamat.filt.s[i,] = (data.mat.filt[i,]-mean(data.mat.filt[i,]))/sd(data.mat.filt[i,])
}
Match = match(id.met.filt,id.heme)
ID.match.pre = id.heme[Match]
hememat.s.match = hememat2017.s[Match,]
dim(hememat.s.match)

X=hememat.s.match
Y=t(datamat.filt.s)
Y[is.nan(Y)]=0
Y[is.na(Y)]=0

library('PMA')
#X is hematology data
#Y is metabolite data 
lambda.seq = seq(0.1,1,length.out=10) #metabolite penalties considered
scoresXcv <- matrix(nrow=n,ncol=10) #matrix of test scores for hematology
scoresYcv <- matrix(nrow=n,ncol=10) #matrix of test scores for metabolites
cor.vec = c() #test correlation under penalties considered
for(j in 1:10){ for(i in 1:n){ 
  res <- CCA(X[-i,],Y[-i,],penaltyx=1,penaltyz=lambda.seq[j], standardize=FALSE)
  scoresXcv[i,j] <- X[i,]%*%res$u
  scoresYcv[i,j] <- Y[i,]%*%res$v
}
  cor.vec[j]=cor(scoresXcv[,j],scoresYcv[,j])
  }
index = which.max(cor.vec)
#plot of scores under cross-validation:
plot(scoresXcv[,index],scoresYcv[,index])


points(scoresXcv[group.filt=='ID'],scoresYcv[group.filt=='ID'],col='red')
cor(scoresXcv,scoresYcv)
t.test(scoresYcv~group.filt,var.equal=TRUE)
t.test(scoresXcv~group.filt,var.equal=TRUE)
res=CCA(X,Y,penaltyx=1,penaltyz=0.1, standardize=FALSE)
plot(sort(res$v))
metab.post = data.mat[,group=='ID_post']
metab.post.s = metab.post
for(i in 1:654){
  metab.post.s[i,] = (metab.post[i,]-mean(data.mat.filt[i,]))/sd(data.mat.filt[i,])
}
#heme post
hememat2017.post <- matrix(as.numeric(as.matrix(heme2017[c(11:25),c(7,13,15,17,19,21,23)])), ncol=7,nrow=15)
#standardize:
hememat2017.post.s = hememat2017.post
for(i in 1:7){
  hememat2017.post.s[,i] = (hememat2017.post[,i]-mean(hememat2017.impute[,i]))/sd(hememat2017.impute[,i])
}
ID.met.post = id.met[group=='ID_post']
ID.hem.post=as.character(heme2017$ID)[11:25]
Match=match(ID.met.post,ID.hem.post)
X.post=hememat2017.post.s[Match,]
ID.match = ID.hem.post[Match]  #pickuphere
Y.post=t(metab.post.s)
Y.post[is.nan(Y.post)]=0
Y.post[is.na(Y.post)]=0
scoresX.post=X.post%*%res$u
scoresY.post=Y.post%*%res$v
##now IS treated data
#order: ZPP  HGB  HCT  MCV  MCH MCHC  RDW
heme.is7mo = read.csv('~/Dropbox/IronMetab/2017_data/Heme_7mo.csv',header=TRUE)
heme.is7mo.mat = as.matrix(heme.is7mo[,c(4,7,8,9,10,11,12)])
heme.is7mo.mat.s = heme.is7mo.mat
for(i in 1:7){
  heme.is7mo.mat.s[,i] = (heme.is7mo.mat[,i]-mean(hememat2017.impute[,i]))/sd(hememat2017.impute[,i])
}
X.ISpost=heme.is7mo.mat.s
scoresX.ISpost=X.ISpost%*%res$u
metab.is7mo = log(data.is7mo)
metab.is7mo.s = metab.is7mo
for(i in 1:654){
  metab.is7mo.s[i,] = (metab.is7mo[i,]-mean(data.mat.filt[i,]))/sd(data.mat.filt[i,])
}
Y.ISpost = t(metab.is7mo.s)
Y.ISpost[is.nan(Y.ISpost)]=0
Y.ISpost[is.na(Y.ISpost)]=0
Y.ISpost[is.infinite(Y.ISpost)]=0
scoresY.ISpost=Y.ISpost%*%res$v

t.test(scoresX.ISpost,scoresXcv[group.filt=='ID'],var.equal=TRUE)
t.test(scoresY.ISpost,scoresYcv[group.filt=='ID'],var.equal=TRUE)
       t.test(scoresXcv~group.filt))

plot(scoresXcv,scoresYcv)
points(scoresX.post,scoresY.post,pch='+')
points(scoresX.ISpost,scoresY.ISpost,pch='*')

###make plots

tiff('~/Dropbox/IronMetab/2017_data/cca.cv.2017_with_IS7mo_imputed.tiff',height=500,width=500,pointsize=20)
x.score.cv <- scoresXcv
y.score.cv <- scoresYcv
plot(x.score.cv,y.score.cv,col='blue',xlab='Hematology score',ylab='Metabolite score', main='CCA (cross-validation)',las=1)#,xlim=c(-1.5,1.5),ylim=c(-1,1))
points(x.score.cv[group.filt=='ID'],y.score.cv[group.filt=='ID'],col='red',pch='o')
points(x.score.cv[group.filt=='IS'],y.score.cv[group.filt=='IS'],col='blue', pch='o')
points(scoresX.post,scoresY.post,pch='+',col='red')
points(scoresX.ISpost,scoresY.ISpost,pch='*',col='blue')
legend(x=-3.3,y=3.2, legend = c("IS 6mo", "IS 7mo","IDA 6mo","Post-IDA"), pch = c('o','*', 'o','+'), col = c('blue','blue','red', 'red'))
dev.off()


allY = rbind(Y,Y.post,Y.ISpost)
allX = rbind(X,X.post,X.ISpost)
group2 = group[-c(15)]
IDS.match = c(ID.match.pre,ID.match)

blah <- prcomp(allX)
png(file='~/Dropbox/IronMetab/2017_data/heme_PCA_7mo.png',width=450,height=300)
plot(blah$x[,1],blah$x[,2],col='white',xlab='Component 1',ylab='Component 2',main = 'Hematology PCA')
points(blah$x[group2=='IS'&IDS.match!='BL97',1],blah$x[group2=='IS'&IDS.match!='BL97',2],col='red')
points(blah$x[group2=='ID'&IDS.match!='BL97',1],blah$x[group2=='ID'&IDS.match!='BL97',2],col='blue')
points(blah$x[group2=='ID_post'&IDS.match!='BL97',1],blah$x[group2=='ID_post'&IDS.match!='BL97',2],col='green')
points(blah$x[group2=='ID_post'&IDS.match!='BL97',1],blah$x[group2=='ID_post',2],col='blue')
points(blah$x[group2=='ID_post'&IDS.match=='BL97',1],blah$x[group2=='ID_post'&IDS.match=='BL97',2],col='green',pch='X')
points(blah$x[group2=='ID'&IDS.match=='BL97',1],blah$x[group2=='ID'&IDS.match=='BL97',2],col='blue',pch='X')
points(blah$x[34:38,1],blah$x[34:38,2],col='yellow')
legend(1.5,4.2,legend = c('IS 6mo','IS 7mo','IDA 6mo','Post-IDA','BL97 pre','BL97 post'),pch = c('o','o','o','o','X','X'),col=c('red','yellow','blue','green','blue','green'))
dev.off()

blah <- prcomp(allY)
png(file='~/Dropbox/IronMetab/2017_data/metab_PCA_7mo.png',width=450,height=300)
plot(blah$x[,1],blah$x[,2],col='white',xlab='Component 1',ylab='Component 2',main = 'Metabolite PCA')
points(blah$x[group2=='IS'&IDS.match!='BL97',1],blah$x[group2=='IS'&IDS.match!='BL97',2],col='red')
points(blah$x[group2=='ID'&IDS.match!='BL97',1],blah$x[group2=='ID'&IDS.match!='BL97',2],col='blue')
points(blah$x[group2=='ID_post'&IDS.match!='BL97',1],blah$x[group2=='ID_post'&IDS.match!='BL97',2],col='green')
points(blah$x[group2=='ID_post'&IDS.match!='BL97',1],blah$x[group2=='ID_post',2],col='blue')
points(blah$x[group2=='ID_post'&IDS.match=='BL97',1],blah$x[group2=='ID_post'&IDS.match=='BL97',2],col='green',pch='X')
points(blah$x[group2=='ID'&IDS.match=='BL97',1],blah$x[group2=='ID'&IDS.match=='BL97',2],col='blue',pch='X')
points(blah$x[34:38,1],blah$x[34:38,2],col='yellow')
legend(-0.5,36,legend = c('IS 6mo','IS 7mo','IDA 6mo','Post-IDA','BL97 pre','BL97 post'),pch = c('o','o','o','o','X','X'),col=c('red','yellow','blue','green','blue','green'))
dev.off()



###Pairwise t-test for hematology data
#datasets
X.ISpost
X.IDpost = X.post
X.ID = X[group.filt=='ID',]
X.IS = X[group.filt=='IS',]
#1=ZPP, 2= HGB, 3=HCT, 4=MCV, 5=MCH, 6=MCHC, 7=RDW
vars = c(2,4,1,7)
##for HGB
pvals = matrix(nrow=6,ncol=4)
for(i in 1:4){
  ii=vars[i]
  #IS 7mno vs IS 6mo
  pvals[1,i]=t.test(X.ISpost[,ii],X.IS[,ii],var.equal=TRUE)$p.value
  pvals[2,i]=t.test(X.ID[,ii],X.ISpost[,ii],var.equal=TRUE)$p.value
  pvals[3,i]=t.test(X.IDpost[,ii],X.ID[,ii],var.equal=TRUE,paired=TRUE)$p.value
  pvals[4,i]=t.test(X.IDpost[,ii],X.IS[,ii],var.equal=TRUE)$p.value
  pvals[5,i]=t.test(X.ID[,ii],X.IS[,ii],var.equal=TRUE)$p.value
  pvals[5,i]=t.test(X.ID[,ii],X.IS[,ii],var.equal=TRUE)$p.value
  #IS 7mo vs post-IDA missing
  pvals[6,i]=t.test(X.IDpost[,ii],X.ISpost[,ii],var.equal=TRUE)$p.value
}
pvals
#corrected
qvals = matrix(nrow=6,ncol=4)
for(i in 1:4) qvals[,i] = p.adjust(pvals[,i],method='fdr')
qvals

Table = matrix(nrow=6,ncol=8)
Table[,1]=pvals[,1]
Table[,2]=qvals[,1]
Table[,3]=pvals[,2]
Table[,4]=qvals[,2]
Table[,5]=pvals[,3]
Table[,6]=qvals[,3]
Table[,7]=pvals[,4]
Table[,8]=qvals[,4]
Table=round(Table,3)
Table[Table==0]='<0.001'
write.table(Table, file="~/Dropbox/IronMetab/2017_data/heme_comparison_table.txt", sep='\t',row.names=F)
ID7.heme = hememat2017.post
X.post


