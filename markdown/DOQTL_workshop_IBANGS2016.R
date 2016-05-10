## ---- results='hide',message=FALSE,warning=FALSE-------------------------
options(na.action = na.pass)
library(DOQTL)
library(VariantAnnotation)
library(AnnotationHub)

## ----load_pheno,warning=FALSE--------------------------------------------
pheno = readRDS("data/ibangs_phenotypes.rds")

## ----pheno_hist,warning=FALSE--------------------------------------------
hist(pheno$Pct.Time.Corner.Slope, breaks = 20, main = "Pct.Time.Corner.Slope")

## ----load_haplo----------------------------------------------------------
probs = readRDS("data/ibangs_haploprobs.rds")

## ----probs_image,warning=FALSE-------------------------------------------
image(1:ncol(probs), 1:20, t(probs[20:1,,1]), axes = F, ann = F,
      breaks = c(-0.25, 0.25, 0.75, 1.25), col = c("white", "grey50", "black"))
box()
abline(v = 0:9+0.5, col = "grey80")
abline(h = 0:20+0.5, col = "grey80")
mtext(side = 3, line = 0.5, at = 1:8, text = LETTERS[1:8], cex = 1.5)
mtext(side = 2, line = 0.5, at = 20:1, text = rownames(probs)[1:20], las = 1)

## ----load_snps-----------------------------------------------------------
load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
snps = muga_snps[muga_snps[,1] %in% dimnames(probs)[[3]],]
rm(muga_snps)

## ----kinship,message=FALSE,warning=FALSE---------------------------------
K = kinship.probs(probs = probs, snps = snps, bychr = TRUE)

## ----covar,warning=FALSE-------------------------------------------------
covar = model.matrix(~Sex + Generation, data = pheno)[,-1]
colnames(covar)[1] = "sex"

## ----linkage,warning=FALSE-----------------------------------------------
qtl = scanone(pheno = pheno, pheno.col = "Pct.Time.Center.Slope", 
              probs = probs, K = K, addcovar = covar, snps = snps)

## ----linkage_plot,warning=FALSE------------------------------------------
plot(qtl, main = "Pct.Time.Center.Slope")

## ----linkage_perms,warning=FALSE-----------------------------------------
link_perms = readRDS("data/linkage_perms.rds")
thr = get.sig.thr(link_perms, alpha = c(0.05, 0.63))

## ----perms_hist,warning=FALSE--------------------------------------------
hist(link_perms[,1], breaks = 20)
abline(v = thr[1,1], col = "red", lwd = 2)

## ----linkage_plot2,warning=FALSE-----------------------------------------
plot(qtl, main = "Pct.Time.Center.Slope", sig.thr = thr,
     sig.col = c("red", "goldenrod"))

## ----coefplot,warning=FALSE----------------------------------------------
coefplot(qtl, chr = 4, main = "Pct.Time.Center.Slope")

## ----gwas,warning=FALSE--------------------------------------------------
gwas = scanone.assoc(pheno = pheno, pheno.col = "Pct.Time.Center.Slope",
       probs = probs, K = K, addcovar = covar, markers = snps, ncl = 1,
       sdp.file = "data/DO_Sanger_SDPs.txt.bgz")

## ----load_gwas_perms,warnings=FALSE--------------------------------------
assoc_perms = readRDS("data/assoc_perms.rds")
thr = get.sig.thr(-log10(assoc_perms))

## ----gwas_plot,warnings=FALSE--------------------------------------------
plot(gwas, bin.size = 1000, sig.thr = thr)

## ----assoc_map,warning=FALSE---------------------------------------------
assoc = assoc.map(pheno = pheno, pheno.col = "Pct.Time.Center.Slope", 
        probs = probs, K = K[[4]], addcovar = covar, snps = snps, chr = 4,
        start = 149, end = 152, output = "p-value")
high.snps = assoc.plot(assoc, thr = thr[1], show.sdps = TRUE)

## ----get_sanger,warnings=FALSE-------------------------------------------
snp.range = range(high.snps[,2]) * 1e6
snp.file = "ftp://ftp.jax.org/SNPtools/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
header = scanVcfHeader(file = snp.file)
samples = samples(header)[c(5,2,26,28,16,30,35)]
gr = GRanges(seqnames = 4, ranges = IRanges(start = snp.range[1],
             end = snp.range[2]))
param = ScanVcfParam(geno = "GT", samples = samples, which = gr)
vcf = readVcf(file = snp.file, genome = "mm10", param = param)

## ----intersect-snps,warnings=FALSE---------------------------------------
high.snps.gr = GRanges(seqnames = high.snps[,1], ranges = IRanges(start = high.snps[,2] * 1e6, width = 1))
wh = which(!start(high.snps.gr) %in% start(vcf))
start(high.snps.gr)[wh] = start(high.snps.gr)[wh] + 1
all(start(high.snps.gr) %in% start(vcf))
vcf = vcf[start(vcf) %in% start(high.snps.gr)]

## ----conseq,warning=FALSE------------------------------------------------
csq = info(vcf)$CSQ
nonsyn = lapply(csq, function(z) { grep("missense|nonsense|stop|splice", z) })
vcf = vcf[sapply(nonsyn, length) > 0]

## ----hub_query,warning=FALSE---------------------------------------------
hub = query(AnnotationHub(), c("ensembl", "gtf", "mus musculus"))
hub

## ----ensembl_gtf,message=FALSE,warning=FALSE-----------------------------
ensembl = hub[[names(hub)[grep("80", hub$title)]]]

## ----gene_symbols,warning=FALSE------------------------------------------
csq = as.list(info(vcf)$CSQ)
genes = rep(NA, length(csq))
for(i in 1:length(csq)) {
  tmp = csq[[i]][grep("missense|nonsense|stop|splice", csq[[i]])]
  tmp = strsplit(tmp, split = "\\|")
  genes[i] = unique(sapply(tmp, "[", 2))
} # for(i)

## ----query_ensembl,warning=FALSE-----------------------------------------
tmp = ensembl[ensembl$type == "gene" & seqnames(ensembl) == 4]
misense.genes = cbind(genes, tmp$gene_name[match(genes, tmp$gene_id)])
rm(tmp)
misense.genes[!duplicated(misense.genes[,1]),]

## ----load_expr_qtl,warning=FALSE-----------------------------------------
sig.qtl = readRDS("data/expr_sig_qtl.rds")

## ----get_cic_eqtl,warning=FALSE------------------------------------------
high.snps.range = c(min(start(high.snps.gr)), max(start(high.snps.gr))) * 1e-6
chr4.sig.qtl = sig.qtl[sig.qtl$chr == 4 & 
                      sig.qtl$qtl.chr == 4 & 
                      (abs(sig.qtl$start - high.snps.range[1]) < 1 |                                       abs(sig.qtl$start - high.snps.range[2]) < 1),]
chr4.sig.qtl = chr4.sig.qtl[chr4.sig.qtl$qtl.lod > 7,]

## ----load_expr,warning=FALSE---------------------------------------------
expr = readRDS("data/ibangs_expr.rds")

## ----new_cover,warning=FALSE---------------------------------------------
covar2 = covar[rownames(expr),]
expr.subset = expr[,rownames(chr4.sig.qtl)]

## ----mediation_loop,warning=FALSE,eval=FALSE-----------------------------
## pv.expr = rep(0, ncol(expr.subset))
## names(pv.expr) = colnames(expr.subset)
## for(i in 1:ncol(expr.subset)) {
## 
##   local.covar = cbind(covar2, expr.subset[,i])
##   assoc = assoc.map(pheno = pheno, pheno.col = "Pct.Time.Center.Slope",
##         probs = probs, K = K[[4]], addcovar = local.covar, snps = snps, chr = 4,
##         start = 149, end = 150, output = "p-value")
##   pv.expr[i] = -log10(min(assoc[,12]))
## 
## } # for(i)

## ----map_neutral_gene,warning=FALSE--------------------------------------
i = 1
local.covar = cbind(covar2, expr.subset[,i])
assoc1 = assoc.map(pheno = pheno, pheno.col = "Pct.Time.Center.Slope", 
        probs = probs, K = K[[4]], addcovar = local.covar, snps = snps, chr = 4,
        start = 149, end = 150, output = "p-value")
tmp = assoc.plot(assoc1)

## ----map_causal_gene,warning=FALSE---------------------------------------
i = 8
local.covar = cbind(covar2, expr.subset[,i])
assoc2 = assoc.map(pheno = pheno, pheno.col = "Pct.Time.Center.Slope", 
        probs = probs, K = K[[4]], addcovar = local.covar, snps = snps, chr = 4,
        start = 149, end = 150, output = "p-value")
tmp = assoc.plot(assoc2)

