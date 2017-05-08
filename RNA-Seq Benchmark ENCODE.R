setwd("/Volumes/My Book/ShaiM/Alex/")
library("preprocessCore")

star = read.table("countmatrix_star.tsv", stringsAsFactors=F, sep="\t", header=T)
kallisto = read.table("countmatrix_kallisto.tsv", stringsAsFactors=F, sep="\t")

inter = intersect(rownames(star), rownames(kallisto))
cint = intersect(colnames(star), colnames(kallisto))

star = as.matrix(star[inter, cint])
kallisto = as.matrix(kallisto[inter, cint])

star[is.na(star)] = 0
kallisto[is.na(kallisto)] = 0

nstar = log2(star+1)
nstar = normalize.quantiles(nstar)
dimnames(nstar) = dimnames(star)

nkallisto = log2(kallisto+1)
nkallisto = normalize.quantiles(nkallisto)
dimnames(nkallisto) = dimnames(kallisto)


## this is simple correlation between the counts. For the first 20 genes. Each point in the graph is a sample
pdf("count_comparison_genes.pdf")
for(i in 1:20){
  co = cor(star[i,], kallisto[i,])
  plot(star[i,], kallisto[i,], pch='.', cex=2, main=paste(rownames(star)[i], format(co, digits=3)), 
       xlab="STAR", ylab="Kallisto")
  abline(0,1, lwd=2, col="red")
}
dev.off()

# saving a vector which includes all the correlations of the rows
coo = c()
for(i in 1:nrow(star)){
  co = cor(star[i,], kallisto[i,])
  coo = c(coo, co)
}

# plot the density of the correlations
pdf("allgene_cor.pdf")
plot(density(coo, na.rm=T), lwd=3, xlab="gene expression correlation", main="STAR vs Kallisto")
dev.off()


gmt2 = readLines("ENCODE_TF_ChIP-seq_2015.txt")

## get all the TF from encode
ge = list()
for(ll in gmt2){
  sp = unlist(strsplit(ll, "\t"))
  gene = sp[1]
  t = c()
  for(i in 3:length(sp)){
    sp1 = unlist(strsplit( sp[i], ","))
    t = c(t, sp1[1])
  }
  ge[[length(ge)+1]] = intersect(rownames(star), t)
  names(ge)[length(ge)] = gsub("_$", "", gene)
}


## the most frequentely 8000 genes
ue = names(rev(sort(table(unlist(ge))))[1:8000])

## read the experiments file
rr = readLines("TF-SRA-signatures-for-benchmarking.txt")



################################################ 
################### t-test #####################
################################################ 

## function that do t-test
ctest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(0) else return(obj$p.value)
}

sigstar_ttest = list()
sigkallisto_ttest = list()

## getting list with the z-scores of the t-test for each experiment. each term in the list will have the 
## z-scores values of the experiment which is specified in the names of the list.
## Then, plotting kallisto and star for each experiment

pdf("ttest_kallisto_star.pdf")
for(line in rr){
  sp = unlist(strsplit(line, ";"))
  tf = sp[3]
  ctr = unlist(strsplit(gsub("ctrl:","",sp[6]), ","))
  pert = unlist(strsplit(gsub("pert:","",sp[7]), ","))
  
  ic = intersect(ctr, colnames(star))
  ip = intersect(pert, colnames(star))
  if(length(ic) == length(ctr) && length(ip) == length(pert)){
    print(tf)
    tstar <- apply(nstar[,c(ctr,pert)], 1, 
                   function(x){ctest(x[1:length(ctr)],x[(length(ctr)+1):(length(ctr)+length(pert))])})
    tstar <- sort(tstar)
    tkallisto <- apply(nkallisto[,c(ctr,pert)],1,
                       function(x){ctest(x[1:length(ctr)],x[(length(ctr)+1):(length(ctr)+length(pert))])})
    tkallisto <- sort(tkallisto)
    sigstar_ttest[[length(sigstar_ttest)+1]] = tstar
    names(sigstar_ttest)[length(sigstar_ttest)] = paste(sp[2], tf)
    sigkallisto_ttest[[length(sigkallisto_ttest)+1]] = tkallisto
    names(sigkallisto_ttest)[length(sigkallisto_ttest)] = paste(sp[2], tf)
    
    plot((tstar), (tkallisto), pch=".", cex=2, xlab="STAR", ylab="Kallisto", main=paste(sp[2], tf))
    abline(0,1, col="red", lwd=2)
  }
}
dev.off()


pdf("ttest_star_roc_encode.pdf")
aucs_ttest_star_encode = list()
out_star <- c()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigstar_ttest)){
  print(names(sigstar_ttest)[j])
  sp = strsplit(names(sigstar_ttest)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigstar_ttest[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  for(g in names(ge)){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if(!is.nan(aa) && (aa < 0.4 | aa > 0.6)){
        plot(cu, type="l", main=paste(strsplit(g, "_")[[1]][1], " sig:", sp))
        abline(0,length(genelist)/length(ue), col="red")
      }
      auc = c(auc, aa)
    }
    else out_star <- c(g,out_star)
  }
  aucs_ttest_star_encode[[length(aucs_ttest_star_encode)+1]] = auc
  names(aucs_ttest_star_encode)[length(aucs_ttest_star_encode)] = sp
}
dev.off()



pdf("ttest_kallisto_roc_encode.pdf")
aucs_ttest_kallisto_encode = list()
out_kallisto <- c()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigkallisto_ttest)){
  print(names(sigkallisto_ttest)[j])
  sp = strsplit(names(sigkallisto_ttest)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigkallisto_ttest[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  for(g in names(ge)){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if(!is.nan(aa) && (aa < 0.4 | aa > 0.6)){
        plot(cu, type="l", main=paste(strsplit(g, "_")[[1]][1], " sig:", sp))
        abline(0,length(genelist)/length(ue), col="red")
      }
      auc = c(auc, aa)
    }
    else out_kallisto <- c(g,out_kallisto)
  }
  aucs_ttest_kallisto_encode[[length(aucs_ttest_kallisto_encode)+1]] = auc
  names(aucs_ttest_kallisto_encode)[length(aucs_ttest_kallisto_encode)] = sp
}
dev.off()




aum_ttest_star_encode = do.call(cbind, aucs_ttest_star_encode)
write.table(aum_ttest_star_encode, file="aum_ttest_star_encode.txt", sep="\t")
aum_ttest_kallisto_encode = do.call(cbind, aucs_ttest_kallisto_encode)
write.table(aum_ttest_kallisto_encode, file="aum_ttest_kallisto_encode.txt", sep="\t")

dd_ttest_star_encode = c()
for(i in 1:length(tg)){
  # index of column to take (match between the TF gene and the chip seq gene)
  ww = which(colnames(aum_ttest_star_encode) == tg[i])
  
  # save the auc value of this chip-seq and the TF gene
  if(length(ww) > 0){
    dd_ttest_star_encode = c(dd_ttest_star_encode, aum_ttest_star_encode[i,ww])
  }
}

dd_ttest_kallisto_encode = c()
for(i in 1:length(tg)){
  # index of column to take (match between the TF gene and the chip seq gene)
  ww = which(colnames(aum_ttest_kallisto_encode) == tg[i])
  
  # save the auc value of this chip-seq and the TF gene
  if(length(ww) > 0){
    dd_ttest_kallisto_encode = c(dd_ttest_kallisto_encode, aum_ttest_kallisto_encode[i,ww])
  }
}


pdf("ttest_density_star_encode.pdf")
plot(density(aum_ttest_star_encode), lwd=3)
lines(density(dd_ttest_star_encode), lwd=3, col="red")
dev.off()


pdf("ttest_density_kallisto_encode.pdf")
plot(density(aum_ttest_kallisto_encode), lwd=3)
lines(density(dd_ttest_kallisto_encode), lwd=3, col="red")
dev.off()

out_star <- unique(out_star)
out_kallisto <- unique(out_kallisto)
rownames(aum_ttest_star_encode) <- names(ge)
rownames(aum_ttest_kallisto_encode) <- names(ge)



##################################################### 
################### fold change #####################
##################################################### 

## function that do the fold change calculations

fold_change_function <- function(matrix_control, matrix_condition){
  res <- abs(log2(sum(matrix_control)/length(matrix_control)) - log2(sum(matrix_condition)/length(matrix_condition)))
  return(res)
}


sigstar_fch <- list()
sigkallisto_fch <- list()



pdf("fch_kallisto_star.pdf")
for(line in rr){
  sp = unlist(strsplit(line, ";"))
  tf = sp[3]
  ctr = unlist(strsplit(gsub("ctrl:","",sp[6]), ","))
  pert = unlist(strsplit(gsub("pert:","",sp[7]), ","))
  
  ic = intersect(ctr, colnames(star))
  ip = intersect(pert, colnames(star))
  if(length(ic) == length(ctr) && length(ip) == length(pert)){
    print(tf)
    fch_star <- apply(nstar[,c(ctr,pert)], 1, function(x) fold_change_function(x[ctr],x[pert]))
    fch_kallisto <- apply(nkallisto[,c(ctr,pert)], 1, function(x) fold_change_function(x[ctr],x[pert]))
    sigstar_fch[[length(sigstar_fch)+1]] = fch_star
    names(sigstar_fch)[length(sigstar_fch)] = paste(sp[2], tf)
    sigkallisto_fch[[length(sigkallisto_fch)+1]] = fch_kallisto
    names(sigkallisto_fch)[length(sigkallisto_fch)] = paste(sp[2], tf)
    
    plot((fch_star), (fch_kallisto), pch=".", cex=2, xlab="STAR", ylab="Kallisto", main=paste(sp[2], tf))
    abline(0,1, col="red", lwd=2)
  }
}
dev.off()



pdf("fch_star_roc_encode.pdf")
aucs_fch_star_encode = list()
out_star <- c()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigstar_fch)){
  print(names(sigstar_fch)[j])
  sp = strsplit(names(sigstar_fch)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigstar_fch[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  for(g in names(ge)){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if(!is.nan(aa) && (aa < 0.4 | aa > 0.6)){
        plot(cu, type="l", main=paste(strsplit(g, "_")[[1]][1], " sig:", sp))
        abline(0,length(genelist)/length(ue), col="red")
      }
      auc = c(auc, aa)
    }
    else out_star <- c(g,out_star)
  }
  aucs_fch_star_encode[[length(aucs_fch_star_encode)+1]] = auc
  names(aucs_fch_star_encode)[length(aucs_fch_star_encode)] = sp
}
dev.off()



pdf("fch_kallisto_roc_encode.pdf")
aucs_fch_kallisto_encode = list()
out_kallisto <- c()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigkallisto_fch)){
  print(names(sigkallisto_fch)[j])
  sp = strsplit(names(sigkallisto_fch)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigkallisto_fch[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  for(g in names(ge)){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if(!is.nan(aa) && (aa < 0.4 | aa > 0.6)){
        plot(cu, type="l", main=paste(strsplit(g, "_")[[1]][1], " sig:", sp))
        abline(0,length(genelist)/length(ue), col="red")
      }
      auc = c(auc, aa)
    }
    else out_kallisto <- c(g,out_kallisto)
  }
  aucs_fch_kallisto_encode[[length(aucs_fch_kallisto_encode)+1]] = auc
  names(aucs_fch_kallisto_encode)[length(aucs_fch_kallisto_encode)] = sp
}
dev.off()


aum_fch_star_encode = do.call(cbind, aucs_fch_star_encode)
write.table(aum_fch_star_encode, file="aum_fch_star_encode.txt", sep="\t")
aum_fch_kallisto_encode = do.call(cbind, aucs_fch_kallisto_encode)
write.table(aum_fch_kallisto_encode, file="aum_fch_kallisto_encode.txt", sep="\t")

dd_fch_star_encode = c()
for(i in 1:length(tg)){
  # index of column to take (match between the TF gene and the chip seq gene)
  ww = which(colnames(aum_fch_star_encode) == tg[i])
  
  # save the auc value of this chip-seq and the TF gene
  if(length(ww) > 0){
    dd_fch_star_encode = c(dd_fch_star_encode, aum_fch_star_encode[i,ww])
  }
}

dd_fch_kallisto_encode = c()
for(i in 1:length(tg)){
  # index of column to take (match between the TF gene and the chip seq gene)
  ww = which(colnames(aum_fch_kallisto_encode) == tg[i])
  
  # save the auc value of this chip-seq and the TF gene
  if(length(ww) > 0){
    dd_fch_kallisto_encode = c(dd_fch_kallisto_encode, aum_fch_kallisto_encode[i,ww])
  }
}

pdf("fch_density_star_encode.pdf")
plot(density(aum_fch_star_encode), lwd=3)
lines(density(dd_fch_star_encode), lwd=3, col="red")
dev.off()


pdf("fch_density_kallisto_encode.pdf")
plot(density(aum_fch_kallisto_encode), lwd=3)
lines(density(dd_fch_kallisto_encode), lwd=3, col="red")
dev.off()


out_star <- unique(out_star)
out_kallisto <- unique(out_kallisto)
rownames(aum_fch_star_encode) <- names(ge)
rownames(aum_fch_kallisto_encode) <- names(ge)

################################################ 
################### DESeq2 #####################
################################################


DESeq2_function <- function(countData,g1,g2){
  colData <- as.data.frame(c(rep(c("Control"),length(g1)),rep(c("Condition"),length(g2))))
  rownames(colData) <- c(g1,g2)
  colnames(colData) <- c("group")
  colData$group <- relevel(colData$group, ref="Control")
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design=~(group))
  dds$group <- relevel(dds$group, ref="Control")
  dds <- DESeq(dds)
  res <- results(dds)
  res <- as.data.frame(res)
  res <- res[order(res$pvalue),, drop=F]
  res <- res[,'pvalue',drop=F]
  return(res)
}


sigstar_deseq2 <- list()
sigkallisto_deseq2 <- list()


pdf("deseq2_kallisto_star.pdf")
for(line in rr){
  sp = unlist(strsplit(line, ";"))
  tf = sp[3]
  ctr = unlist(strsplit(gsub("ctrl:","",sp[6]), ","))
  pert = unlist(strsplit(gsub("pert:","",sp[7]), ","))
  
  ic = intersect(ctr, colnames(star))
  ip = intersect(pert, colnames(star))
  if(length(ic) == length(ctr) && length(ip) == length(pert)){
    print(tf)
    dss <- DESeq2_function(star[,c(ctr,pert)],ctr,pert)
    deseq2_star <- apply(dss, 1, function(x) x)
    dsk <- DESeq2_function(round(kallisto[,c(ctr,pert)]),ctr,pert)
    deseq2_kallisto <- apply(dsk, 1, function(x) x)
    sigstar_deseq2[[length(sigstar_deseq2)+1]] = deseq2_star
    names(sigstar_deseq2)[length(sigstar_deseq2)] = paste(sp[2], tf)
    sigkallisto_deseq2[[length(sigkallisto_deseq2)+1]] = deseq2_kallisto
    names(sigkallisto_deseq2)[length(sigkallisto_deseq2)] = paste(sp[2], tf)
    
    plot((deseq2_star), (deseq2_kallisto), pch=".", cex=2, xlab="STAR", ylab="Kallisto", main=paste(sp[2], tf))
    abline(0,1, col="red", lwd=2)
  }
}
dev.off()




pdf("deseq2_star_roc_encode.pdf")
aucs_deseq2_star_encode = list()
out_star <- c()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigstar_deseq2)){
  print(names(sigstar_deseq2)[j])
  sp = strsplit(names(sigstar_deseq2)[j], " ")[[1]][2]
  sig = rev(sort(sigstar_deseq2[[j]][ue]))
  auc = c()
  rn = c()
  tg = c()
  
  for(g in names(ge)){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if(!is.nan(aa) && (aa < 0.4 | aa > 0.6)){
        plot(cu, type="l", main=paste(strsplit(g, "_")[[1]][1], " sig:", sp))
        abline(0,length(genelist)/length(ue), col="red")
      }
      auc = c(auc, aa)
    }
    else out_star <- c(g,out_star)
  }
  aucs_deseq2_star_encode[[length(aucs_deseq2_star_encode)+1]] = auc
  names(aucs_deseq2_star_encode)[length(aucs_deseq2_star_encode)] = sp
}
dev.off()



pdf("deseq2_kallisto_roc_encode.pdf")
aucs_deseq2_kallisto_encode = list()
out_kallisto <- c()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigkallisto_deseq2)){
  print(names(sigkallisto_deseq2)[j])
  sp = strsplit(names(sigkallisto_deseq2)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigkallisto_deseq2[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  for(g in names(ge)){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if(!is.nan(aa) && (aa < 0.4 | aa > 0.6)){
        plot(cu, type="l", main=paste(strsplit(g, "_")[[1]][1], " sig:", sp))
        abline(0,length(genelist)/length(ue), col="red")
      }
      auc = c(auc, aa)
    }
    else out_kallisto <- c(out_kallisto)
  }
  aucs_deseq2_kallisto_encode[[length(aucs_deseq2_kallisto_encode)+1]] = auc
  names(aucs_deseq2_kallisto_encode)[length(aucs_deseq2_kallisto_encode)] = sp
}
dev.off()



aum_deseq2_star_encode = do.call(cbind, aucs_deseq2_star_encode)
write.table(aum_deseq2_star_encode, file="aum_deseq2_star_encode.txt", sep="\t")
aum_deseq2_kallisto_encode = do.call(cbind, aucs_deseq2_kallisto_encode)
write.table(aum_deseq2_kallisto_encode, file="aum_deseq2_kallisto_encode.txt", sep="\t")

dd_deseq2_star_encode = c()
for(i in 1:length(tg)){
  # index of column to take (match between the TF gene and the chip seq gene)
  ww = which(colnames(aum_deseq2_star_encode) == tg[i])
  
  # save the auc value of this chip-seq and the TF gene
  if(length(ww) > 0){
    dd_deseq2_star_encode = c(dd_deseq2_star_encode, aum_deseq2_star_encode[i,ww])
  }
}

dd_deseq2_kallisto_encode = c()
for(i in 1:length(tg)){
  # index of column to take (match between the TF gene and the chip seq gene)
  ww = which(colnames(aum_deseq2_kallisto_encode) == tg[i])
  
  # save the auc value of this chip-seq and the TF gene
  if(length(ww) > 0){
    dd_deseq2_kallisto_encode = c(dd_deseq2_kallisto_encode, aum_deseq2_kallisto_encode[i,ww])
  }
}

pdf("deseq2_density_star_encode.pdf")
plot(density(aum_deseq2_star_encode), lwd=3)
lines(density(dd_deseq2_star_encode), lwd=3, col="red")
dev.off()




pdf("deseq2_density_kallisto_encode.pdf")
plot(density(aum_deseq2_kallisto_encode), lwd=3)
lines(density(dd_deseq2_kallisto_encode), lwd=3, col="red")
dev.off()


out_star <- unique(out_star)
out_kallisto <- unique(out_kallisto)
rownames(aum_deseq2_star_encode) <- names(ge)
rownames(aum_deseq2_kallisto_encode) <- names(ge)



############################################### 
################### edgeR #####################
###############################################


edgeR_function <- function(countData,g1,g2){
  colData <- as.data.frame(c(rep(c("Control"),length(g1)),rep(c("Condition"),length(g2))))
  rownames(colData) <- c(g1,g2)
  colnames(colData) <- c("group")
  colData$group <- relevel(colData$group, ref="Control")
  y <- DGEList(counts=countData, group=colData$group)
  y <- calcNormFactors(y)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  et <- exactTest(y)
  res <- topTags(et, n=Inf)
  res <- as.data.frame(res)
  res <- res[order(res$PValue),]
  res <- res[,'PValue',drop=F]
  return(res)
}


sigstar_edgeR <- list()
sigkallisto_edgeR <- list()



pdf("edgeR_kallisto_star.pdf")
for(line in rr){
  sp = unlist(strsplit(line, ";"))
  tf = sp[3]
  ctr = unlist(strsplit(gsub("ctrl:","",sp[6]), ","))
  pert = unlist(strsplit(gsub("pert:","",sp[7]), ","))
  
  ic = intersect(ctr, colnames(star))
  ip = intersect(pert, colnames(star))
  if(length(ic) == length(ctr) && length(ip) == length(pert)){
    print(tf)
    dss <- edgeR_function(star[,c(ctr,pert)],ctr,pert)
    edgeR_star <- apply(dss, 1, function(x) x)
    dsk <- edgeR_function(round(kallisto[,c(ctr,pert)]),ctr,pert)
    edgeR_kallisto <- apply(dsk, 1, function(x) x)
    sigstar_edgeR[[length(sigstar_edgeR)+1]] = edgeR_star
    names(sigstar_edgeR)[length(sigstar_edgeR)] = paste(sp[2], tf)
    sigkallisto_edgeR[[length(sigkallisto_edgeR)+1]] = edgeR_kallisto
    names(sigkallisto_edgeR)[length(sigkallisto_edgeR)] = paste(sp[2], tf)
    
    plot((edgeR_star), (edgeR_kallisto), pch=".", cex=2, xlab="STAR", ylab="Kallisto", main=paste(sp[2], tf))
    abline(0,1, col="red", lwd=2)
  }
}
dev.off()


pdf("edgeR_star_roc_encode.pdf")
aucs_edgeR_star_encode = list()
out_star <- c()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigstar_edgeR)){
  print(names(sigstar_edgeR)[j])
  sp = strsplit(names(sigstar_edgeR)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigstar_edgeR[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  for(g in names(ge)){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if(!is.nan(aa) && (aa < 0.4 | aa > 0.6)){
        plot(cu, type="l", main=paste(strsplit(g, "_")[[1]][1], " sig:", sp))
        abline(0,length(genelist)/length(ue), col="red")
      }
      auc = c(auc, aa)
    }
    else out_star <- c(g, out_star)
  }
  aucs_edgeR_star_encode[[length(aucs_edgeR_star_encode)+1]] = auc
  names(aucs_edgeR_star_encode)[length(aucs_edgeR_star_encode)] = sp
}
dev.off()



pdf("edgeR_kallisto_roc_encode.pdf")
aucs_edgeR_kallisto_encode = list()
out_kallisto <- c()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigkallisto_edgeR)){
  print(names(sigkallisto_edgeR)[j])
  sp = strsplit(names(sigkallisto_edgeR)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigkallisto_edgeR[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  for(g in names(ge)){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if(!is.nan(aa) && (aa < 0.4 | aa > 0.6)){
        plot(cu, type="l", main=paste(strsplit(g, "_")[[1]][1], " sig:", sp))
        abline(0,length(genelist)/length(ue), col="red")
      }
      auc = c(auc, aa)
    }
    else out_kallisto <- c(g, out_kallisto)
  }
  aucs_edgeR_kallisto_encode[[length(aucs_edgeR_kallisto_encode)+1]] = auc
  names(aucs_edgeR_kallisto_encode)[length(aucs_edgeR_kallisto_encode)] = sp
}
dev.off()



aum_edgeR_star_encode = do.call(cbind, aucs_edgeR_star_encode)
write.table(aum_edgeR_star_encode, file="aum_edgeR_star_encode.txt", sep="\t")
aum_edgeR_kallisto_encode = do.call(cbind, aucs_edgeR_kallisto_encode)
write.table(aum_edgeR_kallisto_encode, file="aum_edgeR_kallisto_encode.txt", sep="\t")

dd_edgeR_star_encode = c()
for(i in 1:length(tg)){
  # index of column to take (match between the TF gene and the chip seq gene)
  ww = which(colnames(aum_edgeR_star_encode) == tg[i])
  
  # save the auc value of this chip-seq and the TF gene
  if(length(ww) > 0){
    dd_edgeR_star_encode = c(dd_edgeR_star_encode, aum_edgeR_star_encode[i,ww])
  }
}

dd_edgeR_kallisto_encode = c()
for(i in 1:length(tg)){
  # index of column to take (match between the TF gene and the chip seq gene)
  ww = which(colnames(aum_edgeR_kallisto_encode) == tg[i])
  
  # save the auc value of this chip-seq and the TF gene
  if(length(ww) > 0){
    dd_edgeR_kallisto_encode = c(dd_edgeR_kallisto_encode, aum_edgeR_kallisto_encode[i,ww])
  }
}

pdf("edgeR_density_star_encode.pdf")
plot(density(aum_edgeR_star_encode), lwd=3)
lines(density(dd_edgeR_star_encode), lwd=3, col="red")
dev.off()



pdf("edgeR_density_kallisto_encode.pdf")
plot(density(aum_edgeR_kallisto_encode), lwd=3)
lines(density(dd_edgeR_kallisto_encode), lwd=3, col="red")
dev.off()

out_star <- unique(out_star)
out_kallisto <- unique(out_kallisto)
rownames(aum_edgeR_star_encode) <- names(ge)
rownames(aum_edgeR_kallisto_encode) <- names(ge)


############################################### 
################### limma #####################
###############################################


limma_function <- function(countData,g1,g2){
  colData <- as.data.frame(c(rep(c("Control"),length(g1)),rep(c("Condition"),length(g2))))
  rownames(colData) <- c(g1,g2)
  colnames(colData) <- c("group")
  colData$group <- relevel(colData$group, ref="Control")
  design <- c()
  for(k in colnames(countData)){
    if(colData[which(rownames(colData)==k),'group'] == "Control"){
      design <- c(design,1)
    }
    if(colData[which(rownames(colData)==k),'group'] == "Condition"){
      design <- c(design,-1)
    }
  }
  dm = model.matrix(~as.factor(design))
  dge <- DGEList(counts=countData)
  if (any(dge$samples$lib.size==0)){
    take_out <- which(dge$samples$lib.size==0)
    dge$samples <- dge$samples[-take_out,]
    dge$counts <- dge$counts[,-take_out]
  }
  dm = model.matrix(~as.factor(design))
  dge <- calcNormFactors(dge)
  v <- voom(dge, dm, plot=F)
  fit <- lmFit(v, dm)
  fit <- eBayes(fit)
  res <- toptable(fit,n=Inf)
  res <- as.data.frame(res)
  res <- res[order(res$P.Value),drop=F,]
  res <- res[,'P.Value', drop=F]
  return(res)
}


sigstar_limma <- list()
sigkallisto_limma <- list()


pdf("limma_kallisto_star.pdf")
for(line in rr){
  sp = unlist(strsplit(line, ";"))
  tf = sp[3]
  ctr = unlist(strsplit(gsub("ctrl:","",sp[6]), ","))
  pert = unlist(strsplit(gsub("pert:","",sp[7]), ","))
  
  ic = intersect(ctr, colnames(star))
  ip = intersect(pert, colnames(star))
  if(length(ic) == length(ctr) && length(ip) == length(pert)){
    print(tf)
    dss <- limma_function(star[,c(ctr,pert)],ctr,pert)
    limma_star <- apply(dss, 1, function(x) x)
    dsk <- limma_function(round(kallisto[,c(ctr,pert)]),ctr,pert)
    limma_kallisto <- apply(dsk, 1, function(x) x)
    sigstar_limma[[length(sigstar_limma)+1]] = limma_star
    names(sigstar_limma)[length(sigstar_limma)] = paste(sp[2], tf)
    sigkallisto_limma[[length(sigkallisto_limma)+1]] = limma_kallisto
    names(sigkallisto_limma)[length(sigkallisto_limma)] = paste(sp[2], tf)
    
    plot((limma_star), (limma_kallisto), pch=".", cex=2, xlab="STAR", ylab="Kallisto", main=paste(sp[2], tf))
    abline(0,1, col="red", lwd=2)
  }
}
dev.off()


pdf("limma_star_roc_encode.pdf")
aucs_limma_star_encode = list()
out_star <- c()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigstar_limma)){
  print(names(sigstar_limma)[j])
  sp = strsplit(names(sigstar_limma)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigstar_limma[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  for(g in names(ge)){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if(!is.nan(aa) && (aa < 0.4 | aa > 0.6)){
        plot(cu, type="l", main=paste(strsplit(g, "_")[[1]][1], " sig:", sp))
        abline(0,length(genelist)/length(ue), col="red")
      }
      auc = c(auc, aa)
    }
    else out_star <- c(g, out_star)
  }
  aucs_limma_star_encode[[length(aucs_limma_star_encode)+1]] = auc
  names(aucs_limma_star_encode)[length(aucs_limma_star_encode)] = sp
}
dev.off()



pdf("limma_kallisto_roc_encode.pdf")
aucs_limma_kallisto_encode = list()
out_kallisto <- c()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigkallisto_limma)){
  print(names(sigkallisto_limma)[j])
  sp = strsplit(names(sigkallisto_limma)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigkallisto_limma[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  for(g in names(ge)){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if(!is.nan(aa) && (aa < 0.4 | aa > 0.6)){
        plot(cu, type="l", main=paste(strsplit(g, "_")[[1]][1], " sig:", sp))
        abline(0,length(genelist)/length(ue), col="red")
      }
      auc = c(auc, aa)
    }
    else out_kallisto <- c(g, out_kallisto)
  }
  aucs_limma_kallisto_encode[[length(aucs_limma_kallisto_encode)+1]] = auc
  names(aucs_limma_kallisto_encode)[length(aucs_limma_kallisto_encode)] = sp
}
dev.off()



aum_limma_star_encode = do.call(cbind, aucs_limma_star_encode)
write.table(aum_limma_star_encode, file="aum_limma_star_encode.txt", sep="\t")
aum_limma_kallisto_encode = do.call(cbind, aucs_limma_kallisto_encode)
write.table(aum_limma_kallisto_encode, file="aum_limma_kallisto_encode.txt", sep="\t")

dd_limma_star_encode = c()
for(i in 1:length(tg)){
  # index of column to take (match between the TF gene and the chip seq gene)
  ww = which(colnames(aum_limma_star_encode) == tg[i])
  
  # save the auc value of this chip-seq and the TF gene
  if(length(ww) > 0){
    dd_limma_star_encode = c(dd_limma_star_encode, aum_limma_star_encode[i,ww])
  }
}

dd_limma_kallisto_encode = c()
for(i in 1:length(tg)){
  # index of column to take (match between the TF gene and the chip seq gene)
  ww = which(colnames(aum_limma_kallisto_encode) == tg[i])
  
  # save the auc value of this chip-seq and the TF gene
  if(length(ww) > 0){
    dd_limma_kallisto_encode = c(dd_limma_kallisto_encode, aum_limma_kallisto_encode[i,ww])
  }
}

pdf("limma_density_star_encode.pdf")
plot(density(aum_limma_star_encode), lwd=3)
lines(density(dd_limma_star_encode), lwd=3, col="red")
dev.off()



pdf("limma_density_kallisto_encode.pdf")
plot(density(aum_limma_kallisto_encode), lwd=3)
lines(density(dd_limma_kallisto_encode), lwd=3, col="red")
dev.off()

out_star <- unique(out_star)
out_kallisto <- unique(out_kallisto)
rownames(aum_limma_star_encode) <- names(ge)
rownames(aum_limma_kallisto_encode) <- names(ge)



############################################ 
################### CD #####################
############################################


"chdirfull" <- function(data,ctrls,expms,npc,r=1)
{
  ctrl = data[,ctrls]
  expm = data[,expms]
  pp = prcomp(t(data))
  last = npc
  V = pp$rotation[,1:last]
  R = pp$x[1:last,1:last]
  pcvars = pp$sdev[1:last]
  meanvec <- rowMeans(expm) - rowMeans(ctrl)
  Dd <- diag(pcvars)
  sigma <- mean(diag(Dd))
  shrunkMats <- r*Dd + sigma*(1-r)*diag(ncol(R))
  b <- V%*%solve(shrunkMats)%*%t(V)%*%meanvec
  b <- b*as.vector(sqrt(1/t(b)%*%b))
  names(b) = rownames(data)
  b
}



cd_function <- function(countData,g1,g2){
  colData <- as.data.frame(c(rep(c("Control"),length(g1)),rep(c("Condition"),length(g2))))
  rownames(colData) <- c(g1,g2)
  colnames(colData) <- c("group")
  colData$group <- relevel(colData$group, ref="Control")
  a <- c()
  for (k in 1:nrow(colData)){
    if (colData$group[k]=='Control'){
      a <-c(a,"0")
    }   
    if(colData$group[k]=='Condition'){
      a <-c(a,"1")
    }
  }
  colnames(countData) <- as.list(a)
  header <- colnames(countData)
  mat <- as.matrix(countData[1:dim(countData)[1],1:dim(countData)[2]])
  ctrlMat <- mat[,header==0]
  expmMat<- cbind(mat[,header==1])
  count_matrix <- cbind(ctrlMat,expmMat)
  unitV <- chdirfull(count_matrix,1:ncol(ctrlMat),(ncol(ctrlMat)+1):(ncol(ctrlMat)+ncol(expmMat)),ncol(count_matrix))
  unitV <- abs(unitV)
  unitV <- cbind(unitV[order(unitV[,1], decreasing = T),])
  unitV <- as.data.frame(unitV)
  return(unitV)
}

sigstar_cd <- list()
sigkallisto_cd <- list()


pdf("cd_kallisto_star.pdf")
for(line in rr){
  sp = unlist(strsplit(line, ";"))
  tf = sp[3]
  ctr = unlist(strsplit(gsub("ctrl:","",sp[6]), ","))
  pert = unlist(strsplit(gsub("pert:","",sp[7]), ","))
  
  ic = intersect(ctr, colnames(star))
  ip = intersect(pert, colnames(star))
  if(length(ic) == length(ctr) && length(ip) == length(pert)){
    print(tf)
    dss <- cd_function(star[,c(ctr,pert)],ctr,pert)
    cd_star <- apply(dss, 1, function(x) x)
    dsk <- cd_function(round(kallisto[,c(ctr,pert)]),ctr,pert)
    cd_kallisto <- apply(dsk, 1, function(x) x)
    sigstar_cd[[length(sigstar_cd)+1]] = cd_star
    names(sigstar_cd)[length(sigstar_cd)] = paste(sp[2], tf)
    sigkallisto_cd[[length(sigkallisto_cd)+1]] = cd_kallisto
    names(sigkallisto_cd)[length(sigkallisto_cd)] = paste(sp[2], tf)
    
    plot((cd_star), (cd_kallisto), pch=".", cex=2, xlab="STAR", ylab="Kallisto", main=paste(sp[2], tf))
    abline(0,1, col="red", lwd=2)
  }
}
dev.off()


pdf("cd_star_roc_encode.pdf")
aucs_cd_star_encode = list()
out_star <- c()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigstar_cd)){
  print(names(sigstar_cd)[j])
  sp = strsplit(names(sigstar_cd)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigstar_cd[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  for(g in names(ge)){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if(!is.nan(aa) && (aa < 0.4 | aa > 0.6)){
        plot(cu, type="l", main=paste(strsplit(g, "_")[[1]][1], " sig:", sp))
        abline(0,length(genelist)/length(ue), col="red")
      }
      auc = c(auc, aa)
    }
    else out_star <- c(g, out_star)
  }
  aucs_cd_star_encode[[length(aucs_cd_star_encode)+1]] = auc
  names(aucs_cd_star_encode)[length(aucs_cd_star_encode)] = sp
}
dev.off()



pdf("cd_kallisto_roc_encode.pdf")
aucs_cd_kallisto_encode = list()
out_kallisto <- c()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigkallisto_cd)){
  print(names(sigkallisto_cd)[j])
  sp = strsplit(names(sigkallisto_cd)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigkallisto_cd[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  for(g in names(ge)){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if(!is.nan(aa) && (aa < 0.4 | aa > 0.6)){
        plot(cu, type="l", main=paste(strsplit(g, "_")[[1]][1], " sig:", sp))
        abline(0,length(genelist)/length(ue), col="red")
      }
      auc = c(auc, aa)
    }
    else out_kallisto <- c(g, out_kallisto)
  }
  aucs_cd_kallisto_encode[[length(aucs_cd_kallisto_encode)+1]] = auc
  names(aucs_cd_kallisto_encode)[length(aucs_cd_kallisto_encode)] = sp
}
dev.off()



aum_cd_star_encode = do.call(cbind, aucs_cd_star_encode)
write.table(aum_cd_star_encode, file="aum_cd_star_encode.txt", sep="\t")
aum_cd_kallisto_encode = do.call(cbind, aucs_cd_kallisto_encode)
write.table(aum_cd_kallisto_encode, file="aum_cd_kallisto_encode.txt", sep="\t")

dd_cd_star_encode = c()
for(i in 1:length(tg)){
  # index of column to take (match between the TF gene and the chip seq gene)
  ww = which(colnames(aum_cd_star_encode) == tg[i])
  
  # save the auc value of this chip-seq and the TF gene
  if(length(ww) > 0){
    dd_cd_star_encode = c(dd_cd_star_encode, aum_cd_star_encode[i,ww])
  }
}


dd_cd_kallisto_encode = c()
for(i in 1:length(tg)){
  # index of column to take (match between the TF gene and the chip seq gene)
  ww = which(colnames(aum_cd_kallisto_encode) == tg[i])
  
  # save the auc value of this chip-seq and the TF gene
  if(length(ww) > 0){
    dd_cd_kallisto_encode = c(dd_cd_kallisto_encode, aum_cd_kallisto_encode[i,ww])
  }
}

pdf("cd_density_star_encode.pdf")
plot(density(aum_cd_star_encode), lwd=3)
lines(density(dd_cd_star_encode), lwd=3, col="red")
dev.off()

pdf("cd_density_kallisto_encode.pdf")
plot(density(aum_cd_kallisto_encode), lwd=3)
lines(density(dd_cd_kallisto_encode), lwd=3, col="red")
dev.off()

out_star <- unique(out_star)
out_kallisto <- unique(out_kallisto)
rownames(aum_cd_star_encode) <- names(ge)
rownames(aum_cd_kallisto_encode) <- names(ge)




####################################################################
##########  roc curves for matched (auc>0.6 | auc<0.4) ############# 
####################################################################

mean(dd_cd_star_encode[which(dd_cd_star_encode>=0.6)])
mean(dd_cd_star_encode[which(dd_cd_star_encode<=0.4)])

# t test
aucs_ttest_star_encode_test = list()
aucs_ttest_star_encode_test1 <- list()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigstar_ttest)){
  print(names(sigstar_ttest)[j])
  sp = strsplit(names(sigstar_ttest)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigstar_ttest[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  go_test <- names(ge)[which(grepl(sp,names(ge)))]
  for(g in go_test){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if((aa < 0.4)){
        aucs_ttest_star_encode_test[[length(aucs_ttest_star_encode_test)+1]] = cu
        names(aucs_ttest_star_encode_test)[length(aucs_ttest_star_encode_test)] = sp
      }
      if(aa>0.6) {
        aucs_ttest_star_encode_test1[[length(aucs_ttest_star_encode_test1)+1]] = cu
        names(aucs_ttest_star_encode_test1)[length(aucs_ttest_star_encode_test1)] = sp
      }
    }
  }
}

if(length(aucs_ttest_star_encode_test)>0){
  test_ttest <- c()
  aucs_ttest_star_encode_test <- as.data.frame(aucs_ttest_star_encode_test)
  for(i in 1:nrow(aucs_ttest_star_encode_test)){
    test_ttest <- c(test_ttest,(sum(aucs_ttest_star_encode_test[i,])/ncol(aucs_ttest_star_encode_test)))
  }
  test_ttest <- rescale(test_ttest, to=c(0,1))
}

if(length(aucs_ttest_star_encode_test1)>0){
  test_ttest1 <- c()
  aucs_ttest_star_encode_test1 <- as.data.frame(aucs_ttest_star_encode_test1)
  for(i in 1:nrow(aucs_ttest_star_encode_test1)){
    test_ttest1 <- c(test_ttest1,(sum(aucs_ttest_star_encode_test1[i,])/ncol(aucs_ttest_star_encode_test1)))
  }
  test_ttest1 <- rescale(test_ttest1, to=c(0,1))
}


# fch
aucs_fch_star_encode_test = list()
aucs_fch_star_encode_test1 <- list()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigstar_fch)){
  print(names(sigstar_fch)[j])
  sp = strsplit(names(sigstar_fch)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigstar_fch[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  go_test <- names(ge)[which(grepl(sp,names(ge)))]
  for(g in go_test){
    genelist = intersect(ge[[g]], ue)
    
    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if((aa < 0.4)){
        aucs_fch_star_encode_test[[length(aucs_fch_star_encode_test)+1]] = cu
        names(aucs_fch_star_encode_test)[length(aucs_fch_star_encode_test)] = sp
      }
      if(aa>0.6) {
        aucs_fch_star_encode_test1[[length(aucs_fch_star_encode_test1)+1]] = cu
        names(aucs_fch_star_encode_test1)[length(aucs_fch_star_encode_test1)] = sp
      }
    }
  }
}

if(length(aucs_fch_star_encode_test)>0){
  test_fch <- c()
  aucs_fch_star_encode_test <- as.data.frame(aucs_fch_star_encode_test)
  for(i in 1:nrow(aucs_fch_star_encode_test)){
    test_fch <- c(test_fch,(sum(aucs_fch_star_encode_test[i,])/ncol(aucs_fch_star_encode_test)))
  }
  test_fch <- rescale(test_fch, to=c(0,1))
}

if(length(aucs_fch_star_encode_test1)>0){
  test_fch1 <- c()
  aucs_fch_star_encode_test1 <- as.data.frame(aucs_fch_star_encode_test1)
  for(i in 1:nrow(aucs_fch_star_encode_test1)){
    test_fch1 <- c(test_fch1,(sum(aucs_fch_star_encode_test1[i,])/ncol(aucs_fch_star_encode_test1)))
  }
  test_fch1 <- rescale(test_fch1, to=c(0,1))
}


# deseq2
aucs_deseq2_star_encode_test = list()
aucs_deseq2_star_encode_test1 <- list()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigstar_deseq2)){
  print(names(sigstar_deseq2)[j])
  sp = strsplit(names(sigstar_deseq2)[j], " ")[[1]][2]
  sigstar_deseq2[[j]][is.na(sigstar_deseq2[[j]])] <- 1
  sig = rev(sort(sigstar_deseq2[[j]][ue]))
  auc = c()
  rn = c()
  tg = c()
  
  go_test <- names(ge)[which(grepl(sp,names(ge)))]
  for(g in go_test){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if((aa < 0.4)){
        aucs_deseq2_star_encode_test[[length(aucs_deseq2_star_encode_test)+1]] = cu
        names(aucs_deseq2_star_encode_test)[length(aucs_deseq2_star_encode_test)] = sp
      }
      if(aa>0.6) {
        aucs_deseq2_star_encode_test1[[length(aucs_deseq2_star_encode_test1)+1]] = cu
        names(aucs_deseq2_star_encode_test1)[length(aucs_deseq2_star_encode_test1)] = sp
      }
    }
  }
}

if(length(aucs_deseq2_star_encode_test)>0){
  test_deseq2 <- c()
  aucs_deseq2_star_encode_test <- as.data.frame(aucs_deseq2_star_encode_test)
  for(i in 1:nrow(aucs_deseq2_star_encode_test)){
    test_deseq2 <- c(test_deseq2,(sum(aucs_deseq2_star_encode_test[i,])/ncol(aucs_deseq2_star_encode_test)))
  }
  test_deseq2 <- rescale(test_deseq2, to=c(0,1))
}

if(length(aucs_deseq2_star_encode_test1)>0){
  test_deseq21 <- c()
  aucs_deseq2_star_encode_test1 <- as.data.frame(aucs_deseq2_star_encode_test1)
  for(i in 1:nrow(aucs_deseq2_star_encode_test1)){
    test_deseq21 <- c(test_deseq21,(sum(aucs_deseq2_star_encode_test1[i,])/ncol(aucs_deseq2_star_encode_test1)))
  }
  test_deseq21 <- rescale(test_deseq21, to=c(0,1))
}


# edgeR
aucs_edgeR_star_encode_test = list()
aucs_edgeR_star_encode_test1 <- list()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigstar_edgeR)){
  print(names(sigstar_edgeR)[j])
  sp = strsplit(names(sigstar_edgeR)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigstar_edgeR[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  go_test <- names(ge)[which(grepl(sp,names(ge)))]
  for(g in go_test){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if((aa < 0.4)){
        aucs_edgeR_star_encode_test[[length(aucs_edgeR_star_encode_test)+1]] = cu
        names(aucs_edgeR_star_encode_test)[length(aucs_edgeR_star_encode_test)] = sp
      }
      if(aa>0.6) {
        aucs_edgeR_star_encode_test1[[length(aucs_edgeR_star_encode_test1)+1]] = cu
        names(aucs_edgeR_star_encode_test1)[length(aucs_edgeR_star_encode_test1)] = sp
      }
    }
  }
}

if(length(aucs_edgeR_star_encode_test)>0){
  test_edgeR <- c()
  aucs_edgeR_star_encode_test <- as.data.frame(aucs_edgeR_star_encode_test)
  for(i in 1:nrow(aucs_edgeR_star_encode_test)){
    test_edgeR <- c(test_edgeR,(sum(aucs_edgeR_star_encode_test[i,])/ncol(aucs_edgeR_star_encode_test)))
  }
  test_edgeR <- rescale(test_edgeR, to=c(0,1))
}

if(length(aucs_edgeR_star_encode_test1)>0){
  test_edgeR1 <- c()
  aucs_edgeR_star_encode_test1 <- as.data.frame(aucs_edgeR_star_encode_test1)
  for(i in 1:nrow(aucs_edgeR_star_encode_test1)){
    test_edgeR1 <- c(test_edgeR1,(sum(aucs_edgeR_star_encode_test1[i,])/ncol(aucs_edgeR_star_encode_test1)))
  }
  test_edgeR1 <- rescale(test_edgeR1, to=c(0,1))
}



# limma
aucs_limma_star_encode_test = list()
aucs_limma_star_encode_test1 <- list()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigstar_limma)){
  print(names(sigstar_limma)[j])
  sp = strsplit(names(sigstar_limma)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigstar_limma[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  go_test <- names(ge)[which(grepl(sp,names(ge)))]
  for(g in go_test){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if((aa < 0.4)){
        aucs_limma_star_encode_test[[length(aucs_limma_star_encode_test)+1]] = cu
        names(aucs_limma_star_encode_test)[length(aucs_limma_star_encode_test)] = sp
      }
      if(aa>0.6) {
        aucs_limma_star_encode_test1[[length(aucs_limma_star_encode_test1)+1]] = cu
        names(aucs_limma_star_encode_test1)[length(aucs_limma_star_encode_test1)] = sp
      }
    }
  }
}

if(length(aucs_limma_star_encode_test)>0){
  test_limma <- c()
  aucs_limma_star_encode_test <- as.data.frame(aucs_limma_star_encode_test)
  for(i in 1:nrow(aucs_limma_star_encode_test)){
    test_limma <- c(test_limma,(sum(aucs_limma_star_encode_test[i,])/ncol(aucs_limma_star_encode_test)))
  }
  test_limma <- rescale(test_limma, to=c(0,1))
}

if(length(aucs_limma_star_encode_test1)>0){
  test_limma1 <- c()
  aucs_limma_star_encode_test1 <- as.data.frame(aucs_limma_star_encode_test1)
  for(i in 1:nrow(aucs_limma_star_encode_test1)){
    test_limma1 <- c(test_limma1,(sum(aucs_limma_star_encode_test1[i,])/ncol(aucs_limma_star_encode_test1)))
  }
  test_limma1 <- rescale(test_limma1, to=c(0,1))
}


# CD
aucs_cd_star_encode_test = list()
aucs_cd_star_encode_test1 <- list()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigstar_cd)){
  print(names(sigstar_cd)[j])
  sp = strsplit(names(sigstar_cd)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigstar_cd[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  go_test <- names(ge)[which(grepl(sp,names(ge)))]
  for(g in go_test){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if((aa < 0.4)){
        aucs_cd_star_encode_test[[length(aucs_cd_star_encode_test)+1]] = cu
        names(aucs_cd_star_encode_test)[length(aucs_cd_star_encode_test)] = sp
      }
      if(aa>0.6) {
        aucs_cd_star_encode_test1[[length(aucs_cd_star_encode_test1)+1]] = cu
        names(aucs_cd_star_encode_test1)[length(aucs_cd_star_encode_test1)] = sp
      }
    }
  }
}

if(length(aucs_cd_star_encode_test)>0){
  test_cd <- c()
  aucs_cd_star_encode_test <- as.data.frame(aucs_cd_star_encode_test)
  for(i in 1:nrow(aucs_cd_star_encode_test)){
    test_cd <- c(test_cd,(sum(aucs_cd_star_encode_test[i,])/ncol(aucs_cd_star_encode_test)))
  }
  test_cd <- rescale(test_cd, to=c(0,1))
}

if(length(aucs_cd_star_encode_test1)>0){
  test_cd1 <- c()
  aucs_cd_star_encode_test1 <- as.data.frame(aucs_cd_star_encode_test1)
  for(i in 1:nrow(aucs_cd_star_encode_test1)){
    test_cd1 <- c(test_cd1,(sum(aucs_cd_star_encode_test1[i,])/ncol(aucs_cd_star_encode_test1)))
  }
  test_cd1 <- rescale(test_cd1, to=c(0,1))
}




x_values <- c(1:8000)
x_values <- rescale(x_values, to=c(0,1))

test11 <- cbind(x_values, test_cd, test_edgeR, test_fch, test_limma, test_ttest)
head(test11)
test11 <- as.data.frame(test11)
names(test11) <- c("FPR","1","2","3","4","5")
newData <- melt(test11, id.vars="FPR")
newData$method <- rep(c("CD","edgeR","Fold Change","limma","t-test"),
                      times = c(length(test_cd),length(test_edgeR),length(test_fch), 
                                length(test_limma), length(test_ttest)))
colnames(newData) <- c("FPR","variable","TPR","method")


pdf("Benchmarking STAR ENCODE AUC<0.4.pdf")
title <- c("Benchmarking STAR ENCODE, AUC<0.4")
cols <- c("blue", "green","cyan","brown","magenta")
g <- ggplot(data=newData, aes(x=FPR, y=TPR, group=variable, colour=method)) + geom_line() + geom_blank()+ 
  theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(title) +
  geom_abline(linetype='dashed',slope=1,intercept=0)+  
  theme(axis.text=element_text(size=20), axis.title=element_text(size=26,face="bold"))+
  theme(plot.title = element_text(size=21))
print(g)
dev.off()


test12 <- cbind(x_values, test_cd1, test_deseq21, test_edgeR1, test_fch1, test_limma1, test_ttest1)
head(test12)
test12 <- as.data.frame(test12)
names(test12) <- c("FPR","1","2","3","4","5","6")
newData <- melt(test12, id.vars="FPR")
newData$method <- rep(c("CD","DESeq2","edgeR","Fold Change","limma","t-test"),
                      times = c(length(test_cd1),length(test_deseq21),length(test_edgeR1),length(test_fch1), 
                                length(test_limma1), length(test_ttest1)))
colnames(newData) <- c("FPR","variable","TPR","method")


pdf("Benchmarking STAR ENCODE AUC>0.6.pdf")
title <- c("Benchmarking STAR ENCODE, AUC>0.6")
cols <- c("blue", "darkred", "forestgreen","black","darkornage","magenta")
g <- ggplot(data=newData, aes(x=FPR, y=TPR, group=variable, colour=method)) + geom_line() + geom_blank()+ 
  theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(title) +
  geom_abline(linetype='dashed',slope=1,intercept=0)+  
  theme(axis.text=element_text(size=20), axis.title=element_text(size=26,face="bold"))+
  theme(plot.title = element_text(size=21))
print(g)
dev.off()



#############################################################################
##########  roc curves for matched (auc>0.6 | auc<0.4) Kallisto ############# 
#############################################################################

mean(dd_cd_kallisto_encode[which(dd_cd_kallisto_encode>=0.6)])
mean(dd_cd_kallisto_encode[which(dd_cd_kallisto_encode<=0.4)])

# t test
aucs_ttest_kallisto_encode_test = list()
aucs_ttest_kallisto_encode_test1 <- list()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigkallisto_ttest)){
  print(names(sigkallisto_ttest)[j])
  sp = strsplit(names(sigkallisto_ttest)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigkallisto_ttest[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  go_test <- names(ge)[which(grepl(sp,names(ge)))]
  for(g in go_test){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if((aa < 0.4)){
        aucs_ttest_kallisto_encode_test[[length(aucs_ttest_kallisto_encode_test)+1]] = cu
        names(aucs_ttest_kallisto_encode_test)[length(aucs_ttest_kallisto_encode_test)] = sp
      }
      if(aa>0.6) {
        aucs_ttest_kallisto_encode_test1[[length(aucs_ttest_kallisto_encode_test1)+1]] = cu
        names(aucs_ttest_kallisto_encode_test1)[length(aucs_ttest_kallisto_encode_test1)] = sp
      }
    }
  }
}

if(length(aucs_ttest_kallisto_encode_test)>0){
  test_ttest <- c()
  aucs_ttest_kallisto_encode_test <- as.data.frame(aucs_ttest_kallisto_encode_test)
  for(i in 1:nrow(aucs_ttest_kallisto_encode_test)){
    test_ttest <- c(test_ttest,(sum(aucs_ttest_kallisto_encode_test[i,])/ncol(aucs_ttest_kallisto_encode_test)))
  }
  test_ttest <- rescale(test_ttest, to=c(0,1))
}

if(length(aucs_ttest_kallisto_encode_test1)>0){
  test_ttest1 <- c()
  aucs_ttest_kallisto_encode_test1 <- as.data.frame(aucs_ttest_kallisto_encode_test1)
  for(i in 1:nrow(aucs_ttest_kallisto_encode_test1)){
    test_ttest1 <- c(test_ttest1,(sum(aucs_ttest_kallisto_encode_test1[i,])/ncol(aucs_ttest_kallisto_encode_test1)))
  }
  test_ttest1 <- rescale(test_ttest1, to=c(0,1))
}


# fch
aucs_fch_kallisto_encode_test = list()
aucs_fch_kallisto_encode_test1 <- list()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigkallisto_fch)){
  print(names(sigkallisto_fch)[j])
  sp = strsplit(names(sigkallisto_fch)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigkallisto_fch[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  go_test <- names(ge)[which(grepl(sp,names(ge)))]
  for(g in go_test){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if((aa < 0.4)){
        aucs_fch_kallisto_encode_test[[length(aucs_fch_kallisto_encode_test)+1]] = cu
        names(aucs_fch_kallisto_encode_test)[length(aucs_fch_kallisto_encode_test)] = sp
      }
      if(aa>0.6) {
        aucs_fch_kallisto_encode_test1[[length(aucs_fch_kallisto_encode_test1)+1]] = cu
        names(aucs_fch_kallisto_encode_test1)[length(aucs_fch_kallisto_encode_test1)] = sp
      }
    }
  }
}

if(length(aucs_fch_kallisto_encode_test)>0){
  test_fch <- c()
  aucs_fch_kallisto_encode_test <- as.data.frame(aucs_fch_kallisto_encode_test)
  for(i in 1:nrow(aucs_fch_kallisto_encode_test)){
    test_fch <- c(test_fch,(sum(aucs_fch_kallisto_encode_test[i,])/ncol(aucs_fch_kallisto_encode_test)))
  }
  test_fch <- rescale(test_fch, to=c(0,1))
}

if(length(aucs_fch_kallisto_encode_test1)>0){
  test_fch1 <- c()
  aucs_fch_kallisto_encode_test1 <- as.data.frame(aucs_fch_kallisto_encode_test1)
  for(i in 1:nrow(aucs_fch_kallisto_encode_test1)){
    test_fch1 <- c(test_fch1,(sum(aucs_fch_kallisto_encode_test1[i,])/ncol(aucs_fch_kallisto_encode_test1)))
  }
  test_fch1 <- rescale(test_fch1, to=c(0,1))
}


# deseq2
aucs_deseq2_kallisto_encode_test = list()
aucs_deseq2_kallisto_encode_test1 <- list()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigkallisto_deseq2)){
  print(names(sigkallisto_deseq2)[j])
  sp = strsplit(names(sigkallisto_deseq2)[j], " ")[[1]][2]
  sigstar_deseq2[[j]][is.na(sigkallisto_deseq2[[j]])] <- 1
  sig = rev(sort(sigkallisto_deseq2[[j]][ue]))
  auc = c()
  rn = c()
  tg = c()
  
  go_test <- names(ge)[which(grepl(sp,names(ge)))]
  for(g in go_test){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if((aa < 0.4)){
        aucs_deseq2_kallisto_encode_test[[length(aucs_deseq2_kallisto_encode_test)+1]] = cu
        names(aucs_deseq2_kallisto_encode_test)[length(aucs_deseq2_kallisto_encode_test)] = sp
      }
      if(aa>0.6) {
        aucs_deseq2_kallisto_encode_test1[[length(aucs_deseq2_kallisto_encode_test1)+1]] = cu
        names(aucs_deseq2_kallisto_encode_test1)[length(aucs_deseq2_kallisto_encode_test1)] = sp
      }
    }
  }
}

if(length(aucs_deseq2_kallisto_encode_test)>0){
  test_deseq2 <- c()
  aucs_deseq2_kallisto_encode_test <- as.data.frame(aucs_deseq2_kallisto_encode_test)
  for(i in 1:nrow(aucs_deseq2_kallisto_encode_test)){
    test_deseq2 <- c(test_deseq2,(sum(aucs_deseq2_kallisto_encode_test[i,])/ncol(aucs_deseq2_kallisto_encode_test)))
  }
  test_deseq2 <- rescale(test_deseq2, to=c(0,1))
}

if(length(aucs_deseq2_kallisto_encode_test1)>0){
  test_deseq21 <- c()
  aucs_deseq2_kallisto_encode_test1 <- as.data.frame(aucs_deseq2_kallisto_encode_test1)
  for(i in 1:nrow(aucs_deseq2_kallisto_encode_test1)){
    test_deseq21 <- c(test_deseq21,(sum(aucs_deseq2_kallisto_encode_test1[i,])/ncol(aucs_deseq2_kallisto_encode_test1)))
  }
  test_deseq21 <- rescale(test_deseq21, to=c(0,1))
}


# edgeR
aucs_edgeR_kallisto_encode_test = list()
aucs_edgeR_kallisto_encode_test1 <- list()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigkallisto_edgeR)){
  print(names(sigkallisto_edgeR)[j])
  sp = strsplit(names(sigkallisto_edgeR)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigkallisto_edgeR[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  go_test <- names(ge)[which(grepl(sp,names(ge)))]
  for(g in go_test){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if((aa < 0.4)){
        aucs_edgeR_kallisto_encode_test[[length(aucs_edgeR_kallisto_encode_test)+1]] = cu
        names(aucs_edgeR_kallisto_encode_test)[length(aucs_edgeR_kallisto_encode_test)] = sp
      }
      if(aa>0.6) {
        aucs_edgeR_kallisto_encode_test1[[length(aucs_edgeR_kallisto_encode_test1)+1]] = cu
        names(aucs_edgeR_kallisto_encode_test1)[length(aucs_edgeR_kallisto_encode_test1)] = sp
      }
    }
  }
}

if(length(aucs_edgeR_kallisto_encode_test)>0){
  test_edgeR <- c()
  aucs_edgeR_kallisto_encode_test <- as.data.frame(aucs_edgeR_kallisto_encode_test)
  for(i in 1:nrow(aucs_edgeR_kallisto_encode_test)){
    test_edgeR <- c(test_edgeR,(sum(aucs_edgeR_kallisto_encode_test[i,])/ncol(aucs_edgeR_kallisto_encode_test)))
  }
  test_edgeR <- rescale(test_edgeR, to=c(0,1))
}

if(length(aucs_edgeR_kallisto_encode_test1)>0){
  test_edgeR1 <- c()
  aucs_edgeR_kallisto_encode_test1 <- as.data.frame(aucs_edgeR_kallisto_encode_test1)
  for(i in 1:nrow(aucs_edgeR_kallisto_encode_test1)){
    test_edgeR1 <- c(test_edgeR1,(sum(aucs_edgeR_kallisto_encode_test1[i,])/ncol(aucs_edgeR_kallisto_encode_test1)))
  }
  test_edgeR1 <- rescale(test_edgeR1, to=c(0,1))
}



# limma
aucs_limma_kallisto_encode_test = list()
aucs_limma_kallisto_encode_test1 <- list()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigstar_limma)){
  print(names(sigkallisto_limma)[j])
  sp = strsplit(names(sigkallisto_limma)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigkallisto_limma[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  go_test <- names(ge)[which(grepl(sp,names(ge)))]
  for(g in go_test){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if((aa < 0.4)){
        aucs_limma_kallisto_encode_test[[length(aucs_limma_kallisto_encode_test)+1]] = cu
        names(aucs_limma_kallisto_encode_test)[length(aucs_limma_kallisto_encode_test)] = sp
      }
      if(aa>0.6) {
        aucs_limma_kallisto_encode_test1[[length(aucs_limma_kallisto_encode_test1)+1]] = cu
        names(aucs_limma_kallisto_encode_test1)[length(aucs_limma_kallisto_encode_test1)] = sp
      }
    }
  }
}

if(length(aucs_limma_kallisto_encode_test)>0){
  test_limma <- c()
  aucs_limma_kallisto_encode_test <- as.data.frame(aucs_limma_kallisto_encode_test)
  for(i in 1:nrow(aucs_limma_kallisto_encode_test)){
    test_limma <- c(test_limma,(sum(aucs_limma_kallisto_encode_test[i,])/ncol(aucs_limma_kallisto_encode_test)))
  }
  test_limma <- rescale(test_limma, to=c(0,1))
}

if(length(aucs_limma_kallisto_encode_test1)>0){
  test_limma1 <- c()
  aucs_limma_kallisto_encode_test1 <- as.data.frame(aucs_limma_kallisto_encode_test1)
  for(i in 1:nrow(aucs_limma_kallisto_encode_test1)){
    test_limma1 <- c(test_limma1,(sum(aucs_limma_kallisto_encode_test1[i,])/ncol(aucs_limma_kallisto_encode_test1)))
  }
  test_limma1 <- rescale(test_limma1, to=c(0,1))
}


# CD
aucs_cd_kallisto_encode_test = list()
aucs_cd_kallisto_encode_test1 <- list()
#building a mtrix with auc values. rows are chip sea data (644), cols are TF genes (44)
for(j in 1:length(sigkallisto_cd)){
  print(names(sigkallisto_cd)[j])
  sp = strsplit(names(sigkallisto_cd)[j], " ")[[1]][2]
  sig = rev(sort(abs(sigkallisto_cd[[j]][ue])))
  auc = c()
  rn = c()
  tg = c()
  
  go_test <- names(ge)[which(grepl(sp,names(ge)))]
  for(g in go_test){
    genelist = intersect(ge[[g]], ue)

    if(length(genelist) > 0){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      if((aa < 0.4)){
        aucs_cd_kallisto_encode_test[[length(aucs_cd_kallisto_encode_test)+1]] = cu
        names(aucs_cd_kallisto_encode_test)[length(aucs_cd_kallisto_encode_test)] = sp
      }
      if(aa>0.6) {
        aucs_cd_kallisto_encode_test1[[length(aucs_cd_kallisto_encode_test1)+1]] = cu
        names(aucs_cd_kallisto_encode_test1)[length(aucs_cd_kallisto_encode_test1)] = sp
      }
    }
  }
}

if(length(aucs_cd_kallisto_encode_test)>0){
  test_cd <- c()
  aucs_cd_kallisto_encode_test <- as.data.frame(aucs_cd_kallisto_encode_test)
  for(i in 1:nrow(aucs_cd_kallisto_encode_test)){
    test_cd <- c(test_cd,(sum(aucs_cd_kallisto_encode_test[i,])/ncol(aucs_cd_kallisto_encode_test)))
  }
  test_cd <- rescale(test_cd, to=c(0,1))
}

if(length(aucs_cd_kallisto_encode_test1)>0){
  test_cd1 <- c()
  aucs_cd_kallisto_encode_test1 <- as.data.frame(aucs_cd_kallisto_encode_test1)
  for(i in 1:nrow(aucs_cd_kallisto_encode_test1)){
    test_cd1 <- c(test_cd1,(sum(aucs_cd_kallisto_encode_test1[i,])/ncol(aucs_cd_kallisto_encode_test1)))
  }
  test_cd1 <- rescale(test_cd1, to=c(0,1))
}




x_values <- c(1:8000)
x_values <- rescale(x_values, to=c(0,1))

test11 <- cbind(x_values, test_cd, test_deseq2, test_edgeR, test_fch, test_limma)
head(test11)
test11 <- as.data.frame(test11)
names(test11) <- c("FPR","1","2","3","4","5")
newData <- melt(test11, id.vars="FPR")
newData$method <- rep(c("CD","DESeq2","edgeR","Fold Change","limma"),
                      times = c(length(test_cd),length(test_deseq2),length(test_edgeR),length(test_fch), 
                                length(test_limma)))
colnames(newData) <- c("FPR","variable","TPR","method")


pdf("Benchmarking Kallisto ENCODE AUC<0.4.pdf")
title <- c("Benchmarking Kallisto ENCODE, AUC<0.4")
cols <- c("blue", "red", "green","cyan","brown","magenta")
g <- ggplot(data=newData, aes(x=FPR, y=TPR, group=variable, colour=method)) + geom_line() + geom_blank()+ 
  theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(title) +
  geom_abline(linetype='dashed',slope=1,intercept=0)+  
  theme(axis.text=element_text(size=20), axis.title=element_text(size=26,face="bold"))+
  theme(plot.title = element_text(size=21))
print(g)
dev.off()


test12 <- cbind(x_values, test_cd1, test_deseq21, test_edgeR1, test_fch1, test_limma1, test_ttest1)
head(test12)
test12 <- as.data.frame(test12)
names(test12) <- c("FPR","1","2","3","4","5","6")
newData <- melt(test12, id.vars="FPR")
newData$method <- rep(c("CD","DESeq2","edgeR","Fold Change","limma","t-test"),
                      times = c(length(test_cd1),length(test_deseq21),length(test_edgeR1),length(test_fch1), 
                                length(test_limma1), length(test_ttest1)))
colnames(newData) <- c("FPR","variable","TPR","method")


pdf("Benchmarking Kallisto ENCODE AUC>0.6.pdf")
title <- c("Benchmarking Kallisto ENCODE, AUC>0.6")
cols <- c("blue", "darkred", "forestgreen","black","darkornage","magenta")
g <- ggplot(data=newData, aes(x=FPR, y=TPR, group=variable, colour=method)) + geom_line() + geom_blank()+ 
  theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(title) +
  geom_abline(linetype='dashed',slope=1,intercept=0)+  
  theme(axis.text=element_text(size=20), axis.title=element_text(size=26,face="bold"))+
  theme(plot.title = element_text(size=21))
print(g)
dev.off()









### construction an histogram graph

plot_ranking_TF <- function(mat, aligner, met){
  x_axis <- rep(0,length(names(ge)))
  
  for(i in 1:ncol(mat)){
    temp <- mat[,i]
    temp <- as.data.frame(temp)
    temp <- temp[order(temp$temp, decreasing = T),,drop=F]
    temp1 <- rownames(temp)[which(temp$temp>0.6 |temp$temp<0.4)]
    matching <- names(ge)[grepl(colnames(mat)[i], names(ge))]
    temp2 <- intersect(temp1,matching)
    if(length(temp2)>0){
      where <- match(temp2, rownames(temp))
      x_axis[where] <- x_axis[where] + 1
    }
  }
  
  
  y_axis <- as.vector(tapply(x_axis, (seq_along(x_axis)-1) %/% 25, sum))
  y_axis <- as.data.frame(y_axis)
  y_axis[,2] <- as.character(c("1-25","26-50","51-75","76-100",
                               "101-125","126-150","151-175","176-200",
                               "201-225","226-250","251-275","276-300",
                               "301-325","326-350","351-375","376-400",
                               "401-425","426-450","451-475","476-500",
                               "501-525","526-550","551-575","576-600",
                               "601-625","626-650","651-675","676-700",
                               "701-725","726-750","751-775","776-800",
                               "801-816"))
  
  colnames(y_axis) <- c("y","rank")
  pdf(sprintf("test-%s-%s.pdf",met,aligner))
  title <- sprintf("Distribution of the rank of the ChIP-Seqs by %s (%s)", met,aligner)
  p <- ggplot(data=y_axis, aes(x=rank , y=y)) +geom_bar(stat = "identity") + 
    theme(axis.text = element_text(size = 3)) +  ggtitle(title) +
    scale_x_discrete(limits=y_axis[,2]) + geom_blank() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x="Rank",y="") + scale_y_continuous(breaks= c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))
  
  print(p)
  dev.off()
}

plot_ranking_TF(aum_ttest_star_encode, aligner="STAR", met="t-test")
plot_ranking_TF(aum_fch_star_encode,  aligner="STAR", met="fold change")
plot_ranking_TF(aum_deseq2_star_encode,  aligner="STAR", met="DESeq2")
plot_ranking_TF(aum_limma_star_encode,  aligner="STAR", met="limma")
plot_ranking_TF(aum_edgeR_star_encode,  aligner="STAR", met="edgeR")
plot_ranking_TF(aum_cd_star_encode,  aligner="STAR", met="CD")
plot_ranking_TF(aum_ttest_kallisto_encode, aligner="Kallisto", met="t-test")
plot_ranking_TF(aum_fch_kallisto_encode, aligner="Kallisto", met="fold change")
plot_ranking_TF(aum_deseq2_kallisto_encode, aligner="Kallisto", met="DESeq2")
plot_ranking_TF(aum_limma_kallisto_encode, aligner="Kallisto", met="limma")
plot_ranking_TF(aum_edgeR_kallisto_encode, aligner="Kallisto", met="edgeR")
plot_ranking_TF(aum_cd_kallisto_encode, aligner="Kallisto", met="CD")

