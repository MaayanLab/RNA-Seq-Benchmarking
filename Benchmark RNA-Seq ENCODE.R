###########################################################
############## RNA-Seq Benchmarking - ENCODE ############## 
###########################################################

setwd("/Volumes/My Book/ShaiM/Alex/")
library("preprocessCore")

# Reading the counts files and keep only the data that aligned successfuly by both STAR and Kallisto
star = read.table("star_full.tsv", stringsAsFactors=F, sep="\t", header=T)
kallisto = read.table("countmatrix_kallisto.tsv", stringsAsFactors=F, sep="\t")


inter = intersect(rownames(star), rownames(kallisto))
cint = intersect(colnames(star), colnames(kallisto))

star = as.matrix(star[inter, cint])
kallisto = as.matrix(kallisto[inter, cint])

star[is.na(star)] = 0
kallisto[is.na(kallisto)] = 0

# Applying log transform on the data (will be used while using the t-test and fold change methods)
# In addition, applying the normalize quantiles function so the data could be compatable 
nstar = log2(star+1)
nstar = normalize.quantiles(nstar)
dimnames(nstar) = dimnames(star)

nkallisto = log2(kallisto+1)
nkallisto = normalize.quantiles(nkallisto)
dimnames(nkallisto) = dimnames(kallisto)


# Correlationo of the raw counts
pdf("count_comparison_genes.pdf") 
for(i in 1:20){
  co = cor(star[i,], kallisto[i,])
  plot(star[i,], kallisto[i,], pch='.', cex=2, main=paste(rownames(star)[i], format(co, digits=3)), 
       xlab="STAR", ylab="Kallisto")
  abline(0,1, lwd=2, col="red")
}
dev.off()

# Coo is a vector which hold the correlations of the raw counts. This vector will be used toplot the density of the correlations
coo = c()
for(i in 1:nrow(star)){
  co = cor(star[i,], kallisto[i,])
  coo = c(coo, co)
}


pdf("allgene_cor.pdf")
plot(density(coo, na.rm=T), lwd=3, xlab="gene expression correlation", main="STAR vs Kallisto")
dev.off()


gmt2 = readLines("ENCODE_TF_ChIP-seq_2015.txt")

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

uu = names(rev(sort(table(unlist(ge)))))

setwd("/Volumes/My Book/ShaiM/Alex/bench_encode/")

rr = readLines("TF studies.txt")

useful1 <- function(vector){
  t <- c()
  for(i in vector){
    ps <- gsub("_.*$", "",i)
    t <- c(t,ps)
  }
  return(t)
}

useful2 <- function(vector, num){
  t <- c()
  for(i in vector){
    ps <- unlist(strsplit(i,"_"))
    t <- c(t,ps[num])
  }
  return(t)
}

useful <- function(vector,num,dd){
  t <- c()
  for(i in vector){
    ps <- unlist(strsplit(i," "))
    t <- c(t,ps[num]) 
  }
  return(t)
}

useful4 <- function(vector,num){
  t <- c()
  for(i in vector){
    ps <- unlist(strsplit(i,"_"))
    t <- c(t,ps[num]) 
  }
  return(t)
}

plot_with_cutoffs <- function(aum,dd){
  vec <- c()
  for(i in 1:ncol(aum)){
    vec <- c(vec, aum[,i])
  }
  eee <- mean(vec)
  eee1 <- sd(vec)
  upper_cutoff <- eee + 2*eee1
  lower_cutoff <- eee - 2*eee1
  print(upper_cutoff)
  print(lower_cutoff)
  plot(density(aum), lwd=3)
  lines(density(dd), lwd=3, col="red")
  abline(v=lower_cutoff, lwd=1, col="blue")
  abline(v=upper_cutoff, lwd=1, col="blue")
  print(length(aum[aum<lower_cutoff]))
  print(length(aum[aum>upper_cutoff]))
  print(length(dd[dd<lower_cutoff]))
  print(length(dd[dd>upper_cutoff]))
}



################################################ 
################### t-test #####################
################################################ 

ctest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(0) else return(obj$p.value)
}

sigstar_ttest_encode = list()
sigkallisto_ttest_encode = list()

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
    tstar <- apply(nstar[,c(ctr,pert)], 1, function(x){ctest(x[1:length(ctr)],x[(length(ctr)+1):(length(ctr)+length(pert))])})
    tstar[is.na(tstar)] <- 1
    tkallisto <- apply(nkallisto[,c(ctr,pert)],1, function(x){ctest(x[1:length(ctr)],x[(length(ctr)+1):(length(ctr)+length(pert))])})
    tkallisto[is.na(tkallisto)] <- 1
    sigstar_ttest_encode[[length(sigstar_ttest_encode)+1]] = tstar
    names(sigstar_ttest_encode)[length(sigstar_ttest_encode)] = paste(sp[2], tf, sp[4],sp[5],sp[8], sp[9])
    sigkallisto_ttest_encode[[length(sigkallisto_ttest_encode)+1]] = tkallisto
    names(sigkallisto_ttest_encode)[length(sigkallisto_ttest_encode)] = paste(sp[2], tf, sp[4],sp[5],sp[8], sp[9])
    
    plot((tstar), (tkallisto), pch=".", cex=2, xlab="STAR", ylab="Kallisto", main=paste(sp[2], tf))
    abline(0,1, col="red", lwd=2)
  }
}
dev.off()

aucs_ttest_star_encode = list()
out_star <- c()
for(j in 1:length(sigstar_ttest_encode)){
  print(names(sigstar_ttest_encode)[j])
  sp <- names(sigstar_ttest_encode)[j]
  sig = sigstar_ttest_encode[[j]][uu]
  sig <- sig[order(sig)]
  auc = c()
  rn = c()
  tg = c()
  for(g in names(ge)){
    genelist = intersect(ge[[g]],names(sig))
    if(length(genelist) > 100){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      auc = c(auc, aa)
    }
    else out_star <- c(g,out_star)
  }
  aucs_ttest_star_encode[[length(aucs_ttest_star_encode)+1]] = auc
  names(aucs_ttest_star_encode)[length(aucs_ttest_star_encode)] = sp
}

aucs_ttest_kallisto_encode = list()
out_kallisto <- c()
for(j in 1:length(sigkallisto_ttest_encode)){
  print(names(sigkallisto_ttest_encode)[j])
  sp <- names(sigkallisto_ttest_encode)[j]
  sig = sigkallisto_ttest_encode[[j]][uu]
  sig <- sig[order(sig)]
  auc = c()
  rn = c()
  tg = c()
  for(g in names(ge)){
    genelist = intersect(ge[[g]],names(sig))
    if(length(genelist) > 100){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      auc = c(auc, aa)
    }
    else out_kallisto <- c(g,out_kallisto)
  }
  aucs_ttest_kallisto_encode[[length(aucs_ttest_kallisto_encode)+1]] = auc
  names(aucs_ttest_kallisto_encode)[length(aucs_ttest_kallisto_encode)] = sp
}

aum_ttest_star_encode = do.call(cbind, aucs_ttest_star_encode)
aum_ttest_kallisto_encode = do.call(cbind, aucs_ttest_kallisto_encode)

dd_ttest_star_encode = c()
for(i in 1:length(tg)){
  temp <- c()
  for(j in 1:length(colnames(aum_ttest_star_encode))){
    temp <- c(temp,strsplit(colnames(aum_ttest_star_encode)[j], " ")[[1]][2])
  }
  ww <- which(temp==tg[i])
  if(length(ww) > 0){
    dd_ttest_star_encode = c(dd_ttest_star_encode, aum_ttest_star_encode[i,ww])
  }
}

pdf("ttest_density_star_encode.pdf")
plot_with_cutoffs(aum_ttest_star_encode,dd_ttest_star_encode)
dev.off()

dd_ttest_kallisto_encode = c()
for(i in 1:length(tg)){
  temp <- c()
  for(j in 1:length(colnames(aum_ttest_kallisto_encode))){
    temp <- c(temp,strsplit(colnames(aum_ttest_kallisto_encode)[j], " ")[[1]][2])
  }
  ww <- which(temp==tg[i])
  if(length(ww) > 0){
    dd_ttest_kallisto_encode = c(dd_ttest_kallisto_encode, aum_ttest_kallisto_encode[i,ww])
  }
}

pdf("ttest_density_kallisto_encode.pdf")
plot_with_cutoffs(aum_ttest_kallisto_encode,dd_ttest_kallisto_encode)
dev.off()

out_star <- unique(out_star)
out_kallisto <- unique(out_kallisto)
rownames(aum_ttest_star_encode) <- names(ge)[-match(out_star,names(ge))]
rownames(aum_ttest_kallisto_encode) <- names(ge)[-match(out_kallisto,names(ge))]
write.table(aum_ttest_star_encode, file="aum_ttest_star_encode.txt", sep="\t")
write.table(aum_ttest_kallisto_encode, file="aum_ttest_kallisto_encode.txt", sep="\t")



##################################################### 
################### fold change #####################
##################################################### 

fold_change_function <- function(matrix_control, matrix_condition){
  res <- (sum(matrix_condition)/length(matrix_condition)) - (sum(matrix_control)/length(matrix_control))
  return(res)
}

sigstar_fch_encode <- list()
sigkallisto_fch_encode <- list()

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
    fch_star <- apply(nstar[,c(ctr,pert)], 1, function(x) fold_change_function(x[1:length(ctr)],x[(length(ctr)+1):(length(ctr)+length(pert))]))
    fch_star[(is.na(fch_star))] <- 0
    fch_kallisto <- apply(nkallisto[,c(ctr,pert)], 1, function(x) fold_change_function(x[1:length(ctr)],x[(length(ctr)+1):(length(ctr)+length(pert))]))
    fch_kallisto[(is.na(fch_kallisto))] <- 0
    sigstar_fch_encode[[length(sigstar_fch_encode)+1]] = fch_star
    names(sigstar_fch_encode)[length(sigstar_fch_encode)] =  paste(sp[2], tf, sp[4],sp[5],sp[8],sp[9])
    sigkallisto_fch_encode[[length(sigkallisto_fch_encode)+1]] = fch_kallisto
    names(sigkallisto_fch_encode)[length(sigkallisto_fch_encode)] =  paste(sp[2], tf, sp[4],sp[5],sp[8],sp[9])
    plot((fch_star), (fch_kallisto), pch=".", cex=2, xlab="STAR", ylab="Kallisto", main=paste(sp[2], tf))
    abline(0,1, col="red", lwd=2)
  }
}
dev.off()

aucs_fch_star_encode = list()
out_star <- c()
for(j in 1:length(sigstar_fch_encode)){
  print(names(sigstar_fch_encode)[j])
  sp <- names(sigstar_fch_encode)[j]
  sig = rev(sort(abs(sigstar_fch_encode[[j]][uu])))
  auc = c()
  rn = c()
  tg = c()
  for(g in names(ge)){
    genelist = intersect(ge[[g]], names(sig))
    if(length(genelist) > 500){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      auc = c(auc, aa)
    }
    else out_star <- c(g,out_star)
  }
  aucs_fch_star_encode[[length(aucs_fch_star_encode)+1]] = auc
  names(aucs_fch_star_encode)[length(aucs_fch_star_encode)] = sp
}

aucs_fch_kallisto_encode = list()
out_kallisto <- c()
for(j in 1:length(sigkallisto_fch_encode)){
  print(names(sigkallisto_fch_encode)[j])
  sp <- names(sigkallisto_fch_encode)[j]
  sig = rev(sort(abs(sigkallisto_fch_encode[[j]][uu])))
  auc = c()
  rn = c()
  tg = c()
  for(g in names(ge)){
    genelist = intersect(ge[[g]], names(sig))
    if(length(genelist) > 500){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      auc = c(auc, aa)
    }
    else out_kallisto <- c(g,out_kallisto)
  }
  aucs_fch_kallisto_encode[[length(aucs_fch_kallisto_encode)+1]] = auc
  names(aucs_fch_kallisto_encode)[length(aucs_fch_kallisto_encode)] = sp
}

aum_fch_star_encode = do.call(cbind, aucs_fch_star_encode)
aum_fch_kallisto_encode = do.call(cbind, aucs_fch_kallisto_encode)

dd_fch_star_encode = c()
for(i in 1:length(tg)){
  temp <- c()
  for(j in 1:length(colnames(aum_fch_star_encode))){
    temp <- c(temp,strsplit(colnames(aum_fch_star_encode)[j], " ")[[1]][2])
  }
  ww <- which(temp==tg[i])
  if(length(ww) > 0){
    dd_fch_star_encode = c(dd_fch_star_encode, aum_fch_star_encode[i,ww])
  }
}

pdf("fch_density_star_encode.pdf")
plot_with_cutoffs(aum_fch_star_encode,dd_fch_star_encode)
dev.off()

dd_fch_kallisto_encode = c()
for(i in 1:length(tg)){
  temp <- c()
  for(j in 1:length(colnames(aum_fch_kallisto_encode))){
    temp <- c(temp,strsplit(colnames(aum_fch_kallisto_encode)[j], " ")[[1]][2])
  }
  ww <- which(temp==tg[i])
  if(length(ww) > 0){
    dd_fch_kallisto_encode = c(dd_fch_kallisto_encode, aum_fch_kallisto_encode[i,ww])
  }
}

pdf("fch_density_kallisto_encode.pdf")
plot_with_cutoffs(aum_fch_kallisto_encode,dd_fch_kallisto_encode)
dev.off()

out_star <- unique(out_star)
out_kallisto <- unique(out_kallisto)
rownames(aum_fch_star_encode) <- names(ge)[-match(out_star,names(ge))]
rownames(aum_fch_kallisto_encode) <- names(ge)[-match(out_kallisto,names(ge))]
write.table(aum_fch_star_encode, file="aum_fch_star_encode.txt", sep="\t")
write.table(aum_fch_kallisto_encode, file="aum_fch_kallisto_encode.txt", sep="\t")



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
  res <- results(dds, independentFiltering=FALSE)
  res[which(is.na(res$padj)),] <- 1
  res <- as.data.frame(res)
  res <- res[,'padj',drop=F]
  return(res)
}

sigstar_deseq2_encode <- list()
sigkallisto_deseq2_encode <- list()

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
    deseq2_star[is.na(deseq2_star)] <- 1
    dsk <- DESeq2_function(round(kallisto[,c(ctr,pert)]),ctr,pert)
    deseq2_kallisto <- apply(dsk, 1, function(x) x)
    deseq2_kallisto[is.na(deseq2_kallisto)] <- 1
    sigstar_deseq2_encode[[length(sigstar_deseq2_encode)+1]] = deseq2_star
    names(sigstar_deseq2_encode)[length(sigstar_deseq2_encode)] = paste(sp[2], tf, sp[4],sp[5],sp[8],sp[9])
    sigkallisto_deseq2_encode[[length(sigkallisto_deseq2_encode)+1]] = deseq2_kallisto
    names(sigkallisto_deseq2_encode)[length(sigkallisto_deseq2_encode)] = paste(sp[2], tf, sp[4],sp[5],sp[8],sp[9])
    
    plot((deseq2_star), (deseq2_kallisto), pch=".", cex=2, xlab="STAR", ylab="Kallisto", main=paste(sp[2], tf))
    abline(0,1, col="red", lwd=2)
  }
}
dev.off()

aucs_deseq2_star_encode = list()
out_star <- c()
for(j in 1:length(sigstar_deseq2_encode)){
  print(names(sigstar_deseq2_encode)[j])
  sp <- names(sigstar_deseq2_encode)[j]
  sig = sort(sigstar_deseq2_encode[[j]][uu])
  auc = c()
  rn = c()
  tg = c()
  for(g in names(ge)){
    genelist = intersect(ge[[g]], names(sig))
    if(length(genelist) > 100){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      auc = c(auc, aa)
    }
    else out_star <- c(g,out_star)
  }
  aucs_deseq2_star_encode[[length(aucs_deseq2_star_encode)+1]] = auc
  names(aucs_deseq2_star_encode)[length(aucs_deseq2_star_encode)] = sp
}

aucs_deseq2_kallisto_encode = list()
out_kallisto <- c()
for(j in 1:length(sigkallisto_deseq2_encode)){
  print(names(sigkallisto_deseq2_encode)[j])
  sp <- names(sigkallisto_deseq2_encode)[j]
  sig = sort(sigkallisto_deseq2_encode[[j]][uu])
  auc = c()
  rn = c()
  tg = c()
  for(g in names(ge)){
    genelist = intersect(ge[[g]], names(sig))
    if(length(genelist) > 100){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      auc = c(auc, aa)
    }
    else out_kallisto <- c(g, out_kallisto)
  }
  aucs_deseq2_kallisto_encode[[length(aucs_deseq2_kallisto_encode)+1]] = auc
  names(aucs_deseq2_kallisto_encode)[length(aucs_deseq2_kallisto_encode)] = sp
}

aum_deseq2_star_encode = do.call(cbind, aucs_deseq2_star_encode)
aum_deseq2_kallisto_encode = do.call(cbind, aucs_deseq2_kallisto_encode)

dd_deseq2_star_encode = c()
for(i in 1:length(tg)){
  temp <- c()
  for(j in 1:length(colnames(aum_deseq2_star_encode))){
    temp <- c(temp,strsplit(colnames(aum_deseq2_star_encode)[j], " ")[[1]][2])
  }
  ww <- which(temp==tg[i])
  if(length(ww) > 0){
    dd_deseq2_star_encode = c(dd_deseq2_star_encode, aum_deseq2_star_encode[i,ww])
  }
}

pdf("deseq2_density_star_encode.pdf")
plot_with_cutoffs(aum_deseq2_star_encode,dd_deseq2_star_encode)
dev.off()

dd_deseq2_kallisto_encode = c()
for(i in 1:length(tg)){
  temp <- c()
  for(j in 1:length(colnames(aum_deseq2_kallisto_encode))){
    temp <- c(temp,strsplit(colnames(aum_deseq2_kallisto_encode)[j], " ")[[1]][2])
  }
  ww <- which(temp==tg[i])
  if(length(ww) > 0){
    dd_deseq2_kallisto_encode = c(dd_deseq2_kallisto_encode, aum_deseq2_kallisto_encode[i,ww])
  }
}

pdf("deseq2_density_kallisto_encode.pdf")
plot_with_cutoffs(aum_deseq2_kallisto_encode,dd_deseq2_kallisto_encode)
dev.off()

out_star <- unique(out_star)
out_kallisto <- unique(out_kallisto)
rownames(aum_deseq2_star_encode) <- names(ge)[-match(out_star,names(ge))]
rownames(aum_deseq2_kallisto_encode) <- names(ge)[-match(out_kallisto,names(ge))]
write.table(aum_deseq2_star_encode, file="aum_deseq2_star_encode.txt", sep="\t")
write.table(aum_deseq2_kallisto_encode, file="aum_deseq2_kallisto_encode.txt", sep="\t")



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
  res <- res[,'FDR',drop=F]
  res <- res[match(names_order, rownames(res)),drop=F,]
  return(res)
}

sigstar_edgeR_encode <- list()
sigkallisto_edgeR_encode <- list()

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
    sigstar_edgeR_encode[[length(sigstar_edgeR_encode)+1]] = edgeR_star
    names(sigstar_edgeR_encode)[length(sigstar_edgeR_encode)] = paste(sp[2], tf, sp[4],sp[5],sp[8],sp[9])
    sigkallisto_edgeR_encode[[length(sigkallisto_edgeR_encode)+1]] = edgeR_kallisto
    names(sigkallisto_edgeR_encode)[length(sigkallisto_edgeR_encode)] = paste(sp[2], tf, sp[4],sp[5],sp[8],sp[9])
    plot((edgeR_star), (edgeR_kallisto), pch=".", cex=2, xlab="STAR", ylab="Kallisto", main=paste(sp[2], tf))
    abline(0,1, col="red", lwd=2)
  }
}
dev.off()

aucs_edgeR_star_encode = list()
out_star <- c()
for(j in 1:length(sigstar_edgeR_encode)){
  print(names(sigstar_edgeR_encode)[j])
  sp <- names(sigstar_edgeR_encode)[j]
  sig = sort(sigstar_edgeR_encode[[j]][uu])
  auc = c()
  rn = c()
  tg = c()
  for(g in names(ge)){
    genelist = intersect(ge[[g]], names(sig))
    if(length(genelist) > 100){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      auc = c(auc, aa)
    }
    else out_star <- c(g, out_star)
  }
  aucs_edgeR_star_encode[[length(aucs_edgeR_star_encode)+1]] = auc
  names(aucs_edgeR_star_encode)[length(aucs_edgeR_star_encode)] = sp
}

aucs_edgeR_kallisto_encode = list()
out_kallisto <- c()
for(j in 1:length(sigkallisto_edgeR_encode)){
  print(names(sigkallisto_edgeR_encode)[j])
  sp <- names(sigkallisto_edgeR_encode)[j]
  sig = sort(sigkallisto_edgeR_encode[[j]][uu])
  auc = c()
  rn = c()
  tg = c()
  for(g in names(ge)){
    genelist = intersect(ge[[g]], names(sig))
    if(length(genelist) > 100){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      auc = c(auc, aa)
    }
    else out_kallisto <- c(g, out_kallisto)
  }
  aucs_edgeR_kallisto_encode[[length(aucs_edgeR_kallisto_encode)+1]] = auc
  names(aucs_edgeR_kallisto_encode)[length(aucs_edgeR_kallisto_encode)] = sp
}

aum_edgeR_star_encode = do.call(cbind, aucs_edgeR_star_encode)
aum_edgeR_kallisto_encode = do.call(cbind, aucs_edgeR_kallisto_encode)

dd_edgeR_star_encode = c()
for(i in 1:length(tg)){
  temp <- c()
  for(j in 1:length(colnames(aum_edgeR_star_encode))){
    temp <- c(temp,strsplit(colnames(aum_edgeR_star_encode)[j], " ")[[1]][2])
  }
  ww <- which(temp==tg[i])
  if(length(ww) > 0){
    dd_edgeR_star_encode = c(dd_edgeR_star_encode, aum_edgeR_star_encode[i,ww])
  }
}

pdf("edgeR_density_star_encode.pdf")
plot_with_cutoffs(aum_edgeR_star_encode,dd_edgeR_star_encode)
dev.off()

dd_edgeR_kallisto_encode = c()
for(i in 1:length(tg)){
  temp <- c()
  for(j in 1:length(colnames(aum_edgeR_kallisto_encode))){
    temp <- c(temp,strsplit(colnames(aum_edgeR_kallisto_encode)[j], " ")[[1]][2])
  }
  ww <- which(temp==tg[i])
  if(length(ww) > 0){
    dd_edgeR_kallisto_encode = c(dd_edgeR_kallisto_encode, aum_edgeR_kallisto_encode[i,ww])
  }
}

pdf("edgeR_density_kallisto_encode.pdf")
plot_with_cutoffs(aum_edgeR_kallisto_encode,dd_edgeR_kallisto_encode)
dev.off()

out_star <- unique(out_star)
out_kallisto <- unique(out_kallisto)
rownames(aum_edgeR_star_encode) <- names(ge)[-match(out_star,names(ge))]
rownames(aum_edgeR_kallisto_encode) <- names(ge)[-match(out_kallisto,names(ge))]
write.table(aum_edgeR_star_encode, file="aum_edgeR_star_encode.txt", sep="\t")
write.table(aum_edgeR_kallisto_encode, file="aum_edgeR_kallisto_encode.txt", sep="\t")



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
  #res <- res[order(res$adj.P.Val),drop=F,]
  res <- res[,'adj.P.Val', drop=F]
  return(res)
}

sigstar_limma_encode <- list()
sigkallisto_limma_encode <- list()

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
    sigstar_limma_encode[[length(sigstar_limma_encode)+1]] = limma_star
    names(sigstar_limma_encode)[length(sigstar_limma_encode)] = paste(sp[2], tf, sp[4],sp[5],sp[8],sp[9])
    sigkallisto_limma_encode[[length(sigkallisto_limma_encode)+1]] = limma_kallisto
    names(sigkallisto_limma_encode)[length(sigkallisto_limma_encode)] = paste(sp[2], tf, sp[4],sp[5],sp[8],sp[9])
    
    plot((limma_star), (limma_kallisto), pch=".", cex=2, xlab="STAR", ylab="Kallisto", main=paste(sp[2], tf))
    abline(0,1, col="red", lwd=2)
  }
}
dev.off()

aucs_limma_star_encode = list()
out_star <- c()
for(j in 1:length(sigstar_limma_encode)){
  print(names(sigstar_limma_encode)[j])
  sp <- names(sigstar_limma_encode)[j]
  sig = sort(sigstar_limma_encode[[j]][uu])
  auc = c()
  rn = c()
  tg = c()
  for(g in names(ge)){
    genelist = intersect(ge[[g]], names(sig))
    if(length(genelist) > 100){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      auc = c(auc, aa)
    }
    else out_star <- c(g, out_star)
  }
  aucs_limma_star_encode[[length(aucs_limma_star_encode)+1]] = auc
  names(aucs_limma_star_encode)[length(aucs_limma_star_encode)] = sp
}

aucs_limma_kallisto_encode = list()
out_kallisto <- c()
for(j in 1:length(sigkallisto_limma_encode)){
  print(names(sigkallisto_limma_encode)[j])
  sp <- names(sigkallisto_limma_encode)[j]
  sig = sort(sigkallisto_limma_encode[[j]][uu])
  auc = c()
  rn = c()
  tg = c()
  for(g in names(ge)){
    genelist = intersect(ge[[g]], names(sig))
    if(length(genelist) > 100){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      auc = c(auc, aa)
    }
    else out_kallisto <- c(g, out_kallisto)
  }
  aucs_limma_kallisto_encode[[length(aucs_limma_kallisto_encode)+1]] = auc
  names(aucs_limma_kallisto_encode)[length(aucs_limma_kallisto_encode)] = sp
}

aum_limma_star_encode = do.call(cbind, aucs_limma_star_encode)
aum_limma_kallisto_encode = do.call(cbind, aucs_limma_kallisto_encode)

dd_limma_star_encode = c()
for(i in 1:length(tg)){
  temp <- c()
  for(j in 1:length(colnames(aum_limma_star_encode))){
    temp <- c(temp,strsplit(colnames(aum_limma_star_encode)[j], " ")[[1]][2])
  }
  ww <- which(temp==tg[i])
  if(length(ww) > 0){
    dd_limma_star_encode = c(dd_limma_star_encode, aum_limma_star_encode[i,ww])
  }
}

pdf("limma_density_star_encode.pdf")
plot_with_cutoffs(aum_limma_star_encode,dd_limma_star_encode)
dev.off()

dd_limma_kallisto_encode = c()
for(i in 1:length(tg)){
  temp <- c()
  for(j in 1:length(colnames(aum_limma_kallisto_encode))){
    temp <- c(temp,strsplit(colnames(aum_limma_kallisto_encode)[j], " ")[[1]][2])
  }
  ww <- which(temp==tg[i])
  if(length(ww) > 0){
    dd_limma_kallisto_encode = c(dd_limma_kallisto_encode, aum_limma_kallisto_encode[i,ww])
  }
}

pdf("limma_density_kallisto_encode.pdf")
plot_with_cutoffs(aum_limma_kallisto_encode,dd_limma_kallisto_encode)
dev.off()

out_star <- unique(out_star)
out_kallisto <- unique(out_kallisto)
rownames(aum_limma_star_encode) <- names(ge)[-match(out_star,names(ge))]
rownames(aum_limma_kallisto_encode) <- names(ge)[-match(out_kallisto,names(ge))]
write.table(aum_limma_star_encode, file="aum_limma_star_encode.txt", sep="\t")
write.table(aum_limma_kallisto_encode, file="aum_limma_kallisto_encode.txt", sep="\t")



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
  unitV <- as.data.frame(unitV)
  return(unitV)
}

sigstar_cd_encode <- list()
sigkallisto_cd_encode <- list()

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
    sigstar_cd_encode[[length(sigstar_cd_encode)+1]] = cd_star
    names(sigstar_cd_encode)[length(sigstar_cd_encode)] = paste(sp[2], tf, sp[4],sp[5],sp[8],sp[9])
    sigkallisto_cd_encode[[length(sigkallisto_cd_encode)+1]] = cd_kallisto
    names(sigkallisto_cd_encode)[length(sigkallisto_cd_encode)] = paste(sp[2], tf, sp[4],sp[5],sp[8],sp[9])
    
    plot((cd_star), (cd_kallisto), pch=".", cex=2, xlab="STAR", ylab="Kallisto", main=paste(sp[2], tf))
    abline(0,1, col="red", lwd=2)
  }
}
dev.off()

aucs_cd_star_encode = list()
out_star <- c()
for(j in 1:length(sigstar_cd_encode)){
  print(names(sigstar_cd_encode)[j])
  sp <- names(sigstar_cd_encode)[j]
  sig = rev(sort(abs(sigstar_cd_encode[[j]][uu])))
  auc = c()
  rn = c()
  tg = c()
  for(g in names(ge)){
    genelist = intersect(ge[[g]], names(sig))
    if(length(genelist) > 100){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      auc = c(auc, aa)
    }
    else out_star <- c(g, out_star)
  }
  aucs_cd_star_encode[[length(aucs_cd_star_encode)+1]] = auc
  names(aucs_cd_star_encode)[length(aucs_cd_star_encode)] = sp
}

aucs_cd_kallisto_encode = list()
out_kallisto <- c()
for(j in 1:length(sigkallisto_cd_encode)){
  print(names(sigkallisto_cd_encode)[j])
  sp <- names(sigkallisto_cd_encode)[j]
  sig = rev(sort(abs(sigkallisto_cd_encode[[j]][uu])))
  auc = c()
  rn = c()
  tg = c()
  for(g in names(ge)){
    genelist = intersect(ge[[g]], names(sig))
    if(length(genelist) > 100){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      auc = c(auc, aa)
    }
    else out_kallisto <- c(g, out_kallisto)
  }
  aucs_cd_kallisto_encode[[length(aucs_cd_kallisto_encode)+1]] = auc
  names(aucs_cd_kallisto_encode)[length(aucs_cd_kallisto_encode)] = sp
}

aum_cd_star_encode = do.call(cbind, aucs_cd_star_encode)
aum_cd_kallisto_encode = do.call(cbind, aucs_cd_kallisto_encode)

dd_cd_star_encode = c()
for(i in 1:length(tg)){
  temp <- c()
  for(j in 1:length(colnames(aum_cd_star_encode))){
    temp <- c(temp,strsplit(colnames(aum_cd_star_encode)[j], " ")[[1]][2])
  }
  ww <- which(temp==tg[i])
  if(length(ww) > 0){
    dd_cd_star_encode = c(dd_cd_star_encode, aum_cd_star_encode[i,ww])
  }
}

pdf("cd_density_star_encode.pdf")
plot_with_cutoffs(aum_cd_star_encode,dd_cd_star_encode)
dev.off()

dd_cd_kallisto_encode = c()
for(i in 1:length(tg)){
  temp <- c()
  for(j in 1:length(colnames(aum_cd_kallisto_encode))){
    temp <- c(temp,strsplit(colnames(aum_cd_kallisto_encode)[j], " ")[[1]][2])
  }
  ww <- which(temp==tg[i])
  if(length(ww) > 0){
    dd_cd_kallisto_encode = c(dd_cd_kallisto_encode, aum_cd_kallisto_encode[i,ww])
  }
}

pdf("cd_density_kallisto_encode.pdf")
plot_with_cutoffs(aum_cd_kallisto_encode,dd_cd_kallisto_encode)
dev.off()

out_star <- unique(out_star)
out_kallisto <- unique(out_kallisto)
rownames(aum_cd_star_encode) <- names(ge)[-match(out_star,names(ge))]
rownames(aum_cd_kallisto_encode) <- names(ge)[-match(out_kallisto,names(ge))]
write.table(aum_cd_star_encode, file="aum_cd_star_encode.txt", sep="\t")
write.table(aum_cd_kallisto_encode, file="aum_cd_kallisto_encode.txt", sep="\t")



############################################ 
################### CD log #################
############################################

cd_log_function <- function(countData,g1,g2){
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
  count_matrix <- log2(count_matrix+1)
  unitV <- chdirfull(count_matrix,1:ncol(ctrlMat),(ncol(ctrlMat)+1):(ncol(ctrlMat)+ncol(expmMat)),ncol(count_matrix))
  unitV <- as.data.frame(unitV)
  return(unitV)
}

sigstar_cd_log_encode <- list()
sigkallisto_cd_log_encode <- list()

pdf("cd_log_kallisto_star.pdf")
for(line in rr){
  sp = unlist(strsplit(line, ";"))
  tf = sp[3]
  ctr = unlist(strsplit(gsub("ctrl:","",sp[6]), ","))
  pert = unlist(strsplit(gsub("pert:","",sp[7]), ","))
 
  ic = intersect(ctr, colnames(star))
  ip = intersect(pert, colnames(star))
  if(length(ic) == length(ctr) && length(ip) == length(pert)){
    print(tf)
    dss <- cd_log_function(star[,c(ctr,pert)],ctr,pert)
    cd_log_star <- apply(dss, 1, function(x) x)
    dsk <- cd_log_function(round(kallisto[,c(ctr,pert)]),ctr,pert)
    cd_log_kallisto <- apply(dsk, 1, function(x) x)
    sigstar_cd_log_encode[[length(sigstar_cd_log_encode)+1]] = cd_log_star
    names(sigstar_cd_log_encode)[length(sigstar_cd_log_encode)] = paste(sp[2], tf, sp[4],sp[5],sp[8],sp[9])
    sigkallisto_cd_log_encode[[length(sigkallisto_cd_log_encode)+1]] = cd_log_kallisto
    names(sigkallisto_cd_log_encode)[length(sigkallisto_cd_log_encode)] = paste(sp[2], tf, sp[4],sp[5],sp[8],sp[9])
    
    plot((cd_log_star), (cd_log_kallisto), pch=".", cex=2, xlab="STAR", ylab="Kallisto", main=paste(sp[2], tf))
    abline(0,1, col="red", lwd=2)
  }
}
dev.off()

aucs_cd_log_star_encode = list()
out_star <- c()
for(j in 1:length(sigstar_cd_log_encode)){
  print(names(sigstar_cd_log_encode)[j])
  sp <- names(sigstar_cd_log_encode)[j]
  sig = rev(sort(abs(sigstar_cd_log_encode[[j]][uu])))
  auc = c()
  rn = c()
  tg = c()
  for(g in names(ge)){
    genelist = intersect(ge[[g]], names(sig))
    if(length(genelist) > 100){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      auc = c(auc, aa)
    }
    else out_star <- c(g, out_star)
  }
  aucs_cd_log_star_encode[[length(aucs_cd_log_star_encode)+1]] = auc
  names(aucs_cd_log_star_encode)[length(aucs_cd_log_star_encode)] = sp
}

aucs_cd_log_kallisto_encode = list()
out_kallisto <- c()
for(j in 1:length(sigkallisto_cd_log_encode)){
  print(names(sigkallisto_cd_log_encode)[j])
  sp <- names(sigkallisto_cd_log_encode)[j]
  sig = rev(sort(abs(sigkallisto_cd_log_encode[[j]][uu])))
  auc = c()
  rn = c()
  tg = c()
  for(g in names(ge)){
    genelist = intersect(ge[[g]], names(sig))
    if(length(genelist) > 100){
      tg = c(tg, strsplit(g, "_")[[1]][1])
      rn = c(rn, g)
      ww = names(sig) %in% genelist
      cu = cumsum(ww)
      aa = sum(cu)/(length(cu)*length(genelist))
      auc = c(auc, aa)
    }
    else out_kallisto <- c(g, out_kallisto)
  }
  aucs_cd_log_kallisto_encode[[length(aucs_cd_log_kallisto_encode)+1]] = auc
  names(aucs_cd_log_kallisto_encode)[length(aucs_cd_log_kallisto_encode)] = sp
}

aum_cd_log_star_encode = do.call(cbind, aucs_cd_log_star_encode)
aum_cd_log_kallisto_encode = do.call(cbind, aucs_cd_log_kallisto_encode)

dd_cd_log_star_encode = c()
for(i in 1:length(tg)){
  temp <- c()
  for(j in 1:length(colnames(aum_cd_log_star_encode))){
    temp <- c(temp,strsplit(colnames(aum_cd_log_star_encode)[j], " ")[[1]][2])
  }
  ww <- which(temp==tg[i])
  if(length(ww) > 0){
    dd_cd_log_star_encode = c(dd_cd_log_star_encode, aum_cd_log_star_encode[i,ww])
  }
}

pdf("cd_log_density_star_encode.pdf")
plot_with_cutoffs(aum_cd_log_star_encode,dd_cd_log_star_encode)
dev.off()

dd_cd_log_kallisto_encode = c()
for(i in 1:length(tg)){
  temp <- c()
  for(j in 1:length(colnames(aum_cd_log_kallisto_encode))){
    temp <- c(temp,strsplit(colnames(aum_cd_log_kallisto_encode)[j], " ")[[1]][2])
  }
  ww <- which(temp==tg[i])
  if(length(ww) > 0){
    dd_cd_log_kallisto_encode = c(dd_cd_log_kallisto_encode, aum_cd_log_kallisto_encode[i,ww])
  }
}

pdf("cd_log_density_kallisto_encode.pdf")
plot_with_cutoffs(aum_cd_log_kallisto_encode,dd_cd_log_kallisto_encode)
dev.off()

out_star <- unique(out_star)
out_kallisto <- unique(out_kallisto)
rownames(aum_cd_log_star_encode) <- names(ge)[-match(out_star,names(ge))]
rownames(aum_cd_log_kallisto_encode) <- names(ge)[-match(out_kallisto,names(ge))]
write.table(aum_cd_log_star_encode, file="aum_cd_log_star_encode.txt", sep="\t")
write.table(aum_cd_log_kallisto_encode, file="aum_cd_log_kallisto_encode.txt", sep="\t")



##############################################################
##################### Benchmark Aligners #####################
##############################################################


reformats <- function(p){
  as.numeric(format(p, digits=3, drop0trailing = T))
}
# plot star vs kallisto for each method while providing the mean AUC and the t.test stats

pdf("ttest_bench.pdf")
plot(density(dd_ttest_star_encode), lwd=3, col="red", main="ttest (star vs. kallisto) ENCODE")
lines(density(dd_ttest_kallisto_encode), lwd=3, col="blue")
abline(v=0.5, lwd=1, col="black")
Corner_text(text=sprintf("mean AUC STAR = %s\nmean AUC Kallisto = %s\np.val=%s",
                         round(mean(dd_ttest_star_encode),digits=3), round(mean(dd_ttest_kallisto_encode),digits=3), 
                         round(ctest(unlist(unname(aucs_ttest_star_encode)),unname(dd_ttest_star_encode)),digits=4)))
legend('right', c("STAR","Kallisto") , lty=1, col=c('red', 'blue'), bty='n', cex=.8)
dev.off()


pdf("fch_bench.pdf")
plot(density(dd_fch_star_encode), lwd=3, col="red", main="fold change (star vs. kallisto) ENCODE")
lines(density(dd_fch_kallisto_encode), lwd=3, col="blue")
abline(v=0.5, lwd=1, col="black")
Corner_text(text=sprintf("mean AUC STAR = %s\nmean AUC Kallisto = %s\np.val=%s",
                         round(mean(dd_fch_star_encode),digits=3), round(mean(dd_fch_kallisto_encode),digits=3), 
                         round(ctest(unlist(unname(aucs_fch_star_encode)),unname(dd_fch_star_encode)),digits=4)))
legend('right', c("STAR","Kallisto") , lty=1, col=c('red', 'blue'), bty='n', cex=.8)
dev.off()


pdf("deseq2_bench.pdf")
plot(density(dd_deseq2_star_encode), lwd=3, col="red", main="DESeq2 (star vs. kallisto) ENCODE")
lines(density(dd_deseq2_kallisto_encode), lwd=3, col="blue")
Corner_text(text=sprintf("mean AUC STAR = %s\nmean AUC Kallisto = %s\np.val=%s",
                         round(mean(dd_deseq2_star_encode),digits=3), round(mean(dd_deseq2_kallisto_encode),digits=3), 
                         round(ctest(unlist(unname(aucs_deseq2_star_encode)),unname(dd_deseq2_star_encode)),digits=4)))
abline(v=0.5, lwd=1, col="black")
legend('right', c("STAR","Kallisto") , lty=1, col=c('red', 'blue'), bty='n', cex=.8)
dev.off()


pdf("edgeR_bench.pdf")
plot(density(dd_edgeR_star_encode), lwd=3, col="red", main="edgeR (star vs. kallisto) ENCODE")
lines(density(dd_edgeR_kallisto_encode), lwd=3, col="blue")
Corner_text(text=sprintf("mean AUC STAR = %s\nmean AUC Kallisto = %s\np.val=%s",
                         round(mean(dd_edgeR_star_encode),digits=3), round(mean(dd_edgeR_kallisto_encode),digits=3), 
                         round(ctest(unlist(unname(aucs_edgeR_star_encode)),unname(dd_edgeR_star_encode)),digits=4)))
abline(v=0.5, lwd=1, col="black")
legend('right', c("STAR","Kallisto") , lty=1, col=c('red', 'blue'), bty='n', cex=.8)
dev.off()


pdf("limma_bench.pdf")
plot(density(dd_limma_star_encode), lwd=3, col="red", main="limma (star vs. kallisto) ENCODE")
lines(density(dd_limma_kallisto_encode), lwd=3, col="blue")
Corner_text(text=sprintf("mean AUC STAR = %s\nmean AUC Kallisto = %s\np.val=%s",
                         round(mean(dd_limma_star_encode),digits=3), round(mean(dd_limma_kallisto_encode),digits=3), 
                         round(ctest(unlist(unname(aucs_limma_star_encode)),unname(dd_limma_star_encode)),digits=4)))
abline(v=0.5, lwd=1, col="black")
legend('right', c("STAR","Kallisto") , lty=1, col=c('red', 'blue'), bty='n', cex=.8)
dev.off()


pdf("cd_bench.pdf")
plot(density(dd_cd_star_encode), lwd=3, col="red", main="CD (star vs. kallisto) ENCODE")
lines(density(dd_cd_kallisto_encode), lwd=3, col="blue")
Corner_text(text=sprintf("mean AUC STAR = %s\nmean AUC Kallisto = %s\np.val=%s",
                         round(mean(dd_cd_star_encode),digits=3), round(mean(dd_cd_kallisto_encode),digits=3), 
                         round(ctest(unlist(unname(aucs_cd_star_encode)),unname(dd_cd_star_encode)),digits=4)))
abline(v=0.5, lwd=1, col="black")
legend('right', c("STAR","Kallisto") , lty=1, col=c('red', 'blue'), bty='n', cex=.8)
dev.off()


pdf("cd_log_bench.pdf")
plot(density(dd_cd_log_star_encode), lwd=3, col="red", main="CD (star vs. kallisto) ENCODE")
lines(density(dd_cd_log_kallisto_encode), lwd=3, col="blue")
Corner_text(text=sprintf("mean AUC STAR = %s\nmean AUC Kallisto = %s\np.val=%s",
                         round(mean(dd_cd_log_star_encode),digits=3), round(mean(dd_cd_log_kallisto_encode),digits=3), 
                         round(ctest(unlist(unname(aucs_cd_log_star_encode)),unname(dd_cd_log_star_encode)),digits=4)))
abline(v=0.5, lwd=1, col="black")
legend('right', c("STAR","Kallisto") , lty=1, col=c('red', 'blue'), bty='n', cex=.8)
dev.off()


###############################################################
###################### Benchmark Methods ######################
###############################################################

# plot all methods in 1 graph 

pdf("star_bench.pdf")
plot(density(dd_ttest_star_encode), lwd=4, col="red", main="STAR (all methods)", xlim=c(0.2,0.9), ylim=c(0,20))
lines(density(dd_fch_star_encode), lwd=4, col="blue")
lines(density(dd_deseq2_star_encode), lwd=4, col="green")
lines(density(dd_edgeR_star_encode), lwd=4, col="purple")
lines(density(dd_limma_star_encode), lwd=4, col="cyan")
lines(density(dd_cd_star_encode), lwd=4, col="orange")
lines(density(dd_cd_log_star_encode), lwd=4, col="orange")
abline(v=0.5, lwd=1, col="black")
legend('right', c("t-test","fold change","DESeq2","edgeR","limma","CD","CD log") , lty=1, 
       col=c('red', 'blue','green','purple','cyan','orange','black'), bty='n', cex=.8)
dev.off()

pdf("kallisto_bench.pdf")
plot(density(dd_ttest_kallisto_encode), lwd=4, col="red", main="Kallisto (all methods)", xlim=c(0.2,0.9), ylim=c(0,20))
lines(density(dd_fch_kallisto_encode), lwd=4, col="blue")
lines(density(dd_deseq2_kallisto_encode), lwd=4, col="green")
lines(density(dd_edgeR_kallisto_encode), lwd=4, col="purple")
lines(density(dd_limma_kallisto_encode), lwd=4, col="cyan")
lines(density(dd_cd_kallisto_encode), lwd=4, col="orange")
abline(v=0.5, lwd=1, col="black")
legend('right', c("t-test","fold change","DESeq2","edgeR","limma","CD") , lty=1, 
       col=c('red', 'blue','green','purple','cyan','orange','black'), bty='n', cex=.8)
dev.off()


################################################################
############ Normalization of matched to background ############
################################################################

normalization_to_background <- function(aucs, dd){
  mean <- mean(unlist(aucs))
  sd <- sd(unlist(aucs))
  temp <- c()
  for(i in 1:length(dd)){
    temp <- c(temp,(dd[[i]]-mean)/sd)
  }
  return(temp)
}

temp1 <- normalization_to_background(aucs_ttest_star_encode, dd_ttest_star_encode)
temp2 <- normalization_to_background(aucs_fch_star_encode, dd_fch_star_encode)
temp3 <- normalization_to_background(aucs_deseq2_star_encode, dd_deseq2_star_encode)
temp4 <- normalization_to_background(aucs_edgeR_star_encode, dd_edgeR_star_encode)
temp5 <- normalization_to_background(aucs_limma_star_encode, dd_limma_star_encode)
temp6 <- normalization_to_background(aucs_cd_star_encode, dd_cd_star_encode)
temp7 <- normalization_to_background(aucs_cd_log_star_encode, dd_cd_log_star_encode)

temp8 <- normalization_to_background(aucs_ttest_kallisto_encode, dd_ttest_kallisto_encode)
temp9 <- normalization_to_background(aucs_fch_kallisto_encode, dd_fch_kallisto_encode)
temp10 <- normalization_to_background(aucs_deseq2_kallisto_encode, dd_deseq2_kallisto_encode)
temp11 <- normalization_to_background(aucs_edgeR_kallisto_encode, dd_edgeR_kallisto_encode)
temp12 <- normalization_to_background(aucs_limma_kallisto_encode, dd_limma_kallisto_encode)
temp13 <- normalization_to_background(aucs_cd_kallisto_encode, dd_cd_kallisto_encode)
temp14 <- normalization_to_background(aucs_cd_log_kallisto_encode, dd_cd_log_kallisto_encode)


# plot all the methods again in 1 graph (this time z-scored values of AUC)

pdf("star all encode.pdf")
plot(density(temp1), lwd=3, col="red", main="STAR (all methods) ENCODE", xlim=c(-5.0,5.0), ylim=c(0.0,1.0))
lines(density(temp2), lwd=3, col="blue")
lines(density(temp3), lwd=3, col="green")
lines(density(temp4), lwd=3, col="purple")
lines(density(temp5), lwd=3, col="cyan")
lines(density(temp6), lwd=3, col="orange")
lines(density(temp7), lwd=3, col="black")
abline(v=c(-1.0,1.0), lwd=1, col="black")
Corner_text(text=sprintf("-1< z-score (AUCs) < 1\nt-test = %s\nfold change = %s\nDESeq2 = %s\nedgeR = %s\nlimma = %s\nCD = %s\nCD log = %s\n",
                         round(mean(temp1[which(temp1<1 & temp1>(-1))]),digits=3),
                         round(mean(temp2[which(temp2<1 & temp2>(-1))]),digits=3),
                         round(mean(temp3[which(temp3<1 & temp3>(-1))]),digits=3),
                         round(mean(temp4[which(temp4<1 & temp4>(-1))]),digits=3),
                         round(mean(temp5[which(temp5<1 & temp5>(-1))]),digits=3),
                         round(mean(temp6[which(temp6<1 & temp6>(-1))]),digits=3),
                         round(mean(temp7[which(temp7<1 & temp7>(-1))]),digits=3)))
legend('right', c("t-test","fold change","DESeq2","edgeR","limma","CD","CD log") , lty=1, 
       col=c('red', 'blue','green','purple','cyan','orange','black'), bty='n', cex=.8)
dev.off()

pdf("kallisto all encode.pdf")
plot(density(temp8), lwd=3, col="red", main="Kallisto (all methods) ENCODE", xlim=c(-5.0,5.0), ylim=c(0.0,1.0))
lines(density(temp9), lwd=3, col="blue")
lines(density(temp10), lwd=3, col="green")
lines(density(temp11), lwd=3, col="purple")
lines(density(temp12), lwd=3, col="cyan")
lines(density(temp13), lwd=3, col="orange")
lines(density(temp14), lwd=3, col="black")
abline(v=0.0, lwd=1, col="black")
Corner_text(text=sprintf("-1< z-score (AUCs) < 1\nt-test = %s\nfold change = %s\nDESeq2 = %s\nedgeR = %s\nlimma = %s\nCD = %s\nCD log = %s\n",
                         round(mean(temp8[which(temp8<1 & temp8>(-1))]),digits=3),
                         round(mean(temp9[which(temp9<1 & temp9>(-1))]),digits=3),
                         round(mean(temp10[which(temp10<1 & temp10>(-1))]),digits=3),
                         round(mean(temp11[which(temp11<1 & temp11>(-1))]),digits=3),
                         round(mean(temp12[which(temp12<1 & temp12>(-1))]),digits=3),
                         round(mean(temp13[which(temp13<1 & temp13>(-1))]),digits=3),
                         round(mean(temp14[which(temp14<1 & temp14>(-1))]),digits=3)))
legend('right', c("t-test","fold change","DESeq2","edgeR","limma","CD","CD log") , lty=1, 
       col=c('red', 'blue','green','purple','cyan','orange','black'), bty='n', cex=.8)
dev.off()



##################################################################################################
############ Gaussian mixer model in attempts to identify subgroups in the population ############
##################################################################################################

pdf('STAR_ttest_ENCODE.pdf')
em <- normalmixEM(as.numeric(temp1), lambda = NULL, k=2)
pop1 <- which(temp3>em$mu[2]-em$sigma[2])
pop2 <- which(temp3<em$mu[1]+em$sigma[1])
plot(em, density=T, loglik = F)
dev.off()

pdf('STAR_fch_ENCODE.pdf')
em <- normalmixEM(as.numeric(temp2), lambda = NULL, k=2)
pop1 <- which(temp3>em$mu[2]-em$sigma[2])
pop2 <- which(temp3<em$mu[1]+em$sigma[1])
plot(em, density=T, loglik = F)
dev.off()

pdf('STAR_deseq2_ENCODE.pdf')
em <- normalmixEM(as.numeric(temp3), lambda = NULL, k=2)
pop1 <- which(temp3>em$mu[2]-em$sigma[2])
pop2 <- which(temp3<em$mu[1]+em$sigma[1])
plot(em, density=T, loglik = F)
dev.off()

pdf('STAR_edgeR_ENCODE.pdf')
em <- normalmixEM(as.numeric(temp4), lambda = NULL, k=2)
pop1 <- which(temp3>em$mu[2]-em$sigma[2])
pop2 <- which(temp3<em$mu[1]+em$sigma[1])
plot(em, density=T, loglik = F)
dev.off()

pdf('STAR_limma_ENCODE.pdf')
em <- normalmixEM(as.numeric(temp5), lambda = NULL, k=2)
pop1 <- which(temp3>em$mu[2]-em$sigma[2])
pop2 <- which(temp3<em$mu[1]+em$sigma[1])
plot(em, density=T, loglik = F)
dev.off()

pdf('STAR_cd_ENCODE.pdf')
em <- normalmixEM(as.numeric(temp6), lambda = NULL, k=2)
pop1 <- which(temp3>em$mu[2]-em$sigma[2])
pop2 <- which(temp3<em$mu[1]+em$sigma[1])
plot(em, density=T, loglik = F)
dev.off()

pdf('STAR_cd_log_ENCODE.pdf')
em <- normalmixEM(as.numeric(temp7), lambda = NULL, k=2)
pop1 <- which(temp3>em$mu[2]-em$sigma[2])
pop2 <- which(temp3<em$mu[1]+em$sigma[1])
plot(em, density=T, loglik = F)
dev.off()



pdf('Kallisto_ttest_ENCODE.pdf')
em <- normalmixEM(as.numeric(temp8), lambda = NULL, k=2)
pop1 <- which(temp3>em$mu[2]-em$sigma[2])
pop2 <- which(temp3<em$mu[1]+em$sigma[1])
plot(em, density=T, loglik = F)
dev.off()

pdf('Kallisto_fch_ENCODE.pdf')
em <- normalmixEM(as.numeric(temp9), lambda = NULL, k=2)
pop1 <- which(temp3>em$mu[2]-em$sigma[2])
pop2 <- which(temp3<em$mu[1]+em$sigma[1])
plot(em, density=T, loglik = F)
dev.off()

pdf('Kallisto_deseq2_ENCODE.pdf')
em <- normalmixEM(as.numeric(temp10), lambda = NULL, k=2)
pop1 <- which(temp3>em$mu[2]-em$sigma[2])
pop2 <- which(temp3<em$mu[1]+em$sigma[1])
plot(em, density=T, loglik = F)
dev.off()

pdf('Kallisto_edgeR_ENCODE.pdf')
em <- normalmixEM(as.numeric(temp11), lambda = NULL, k=2)
pop1 <- which(temp3>em$mu[2]-em$sigma[2])
pop2 <- which(temp3<em$mu[1]+em$sigma[1])
plot(em, density=T, loglik = F)
dev.off()

pdf('Kallisto_limma_ENCODE.pdf')
em <- normalmixEM(as.numeric(temp12), lambda = NULL, k=2)
pop1 <- which(temp3>em$mu[2]-em$sigma[2])
pop2 <- which(temp3<em$mu[1]+em$sigma[1])
plot(em, density=T, loglik = F)
dev.off()

pdf('Kallisto_cd_ENCODE.pdf')
em <- normalmixEM(as.numeric(temp13), lambda = NULL, k=2)
pop1 <- which(temp3>em$mu[2]-em$sigma[2])
pop2 <- which(temp3<em$mu[1]+em$sigma[1])
plot(em, density=T, loglik = F)
dev.off()

pdf('Kallisto_cd_log_ENCODE.pdf')
em <- normalmixEM(as.numeric(temp14), lambda = NULL, k=2)
pop1 <- which(temp3>em$mu[2]-em$sigma[2])
pop2 <- which(temp3<em$mu[1]+em$sigma[1])
plot(em, density=T, loglik = F)
dev.off()
