## analyses for final revision

####### explore mixture models and combined p-values within the HD project
library(ggplot2)
library(mixtools)
library(Exact)
library(coin)
library(EmpiricalBrownsMethod)


# set working dir and read in the data
pvalsfile <- "HomoD_HemiD_4m_rawp_orig.txt"
pvals <- read.table(file = pvalsfile, header = T, sep = "\t", as.is = T, fileEncoding = "latin1")

known_FS <- grep(value = T, pattern = "FRA.*", perl = T, x = pvals$Gene)
known_TS <- c("CDKN2C", "FANCD2/VHL", "CDKN2A", "PTEN", "BRCA2", "RB1", "CYLD", "CDH1", "TP53", "MAP2K4", "NF1", "SMARCB1", "TET1", "FAT1", "BIRC2/BIRC3")

pvals$status <- ifelse(pvals$Gene %in% known_FS, "FS", ifelse(pvals$Gene %in% known_TS, "TS", "unknown"))

### check HD-HemiD counts and reselect tumour type to test
hdhemidfile <- "20170605_HD-HemiD_tumortype_table.txt"
hdhemidcounts <- read.table(file = hdhemidfile, header = T, sep = "\t", as.is = T)[, 1:5]
hdhemidcounts$peakregion <- factor(hdhemidcounts$peakregion, levels = pvals$peakregion)


get_type_pval <- function(df) {
  # max no deletions, break ties biased fashion
  total_idx <- order(df$HDs > 0, df$HDs + df$hemiDs, df$HDs, decreasing = T)[1]
  total_type <- df[total_idx, "tumor_type"]
  total_vect <- c(unlist(df[total_idx, c("HDs", "hemiDs")]), colSums(df[-total_idx, c("HDs", "hemiDs")]))
  total_mat <- matrix(data = total_vect, nrow = 2, byrow = F)
  total_pval <- tryCatch(exact.test(data = total_mat, alternative = "greater", interval = T, method = "Boschloo", model = "Multinomial", to.plot = F)$p.value,
                         error = function(err) {
                           print(paste("ERROR: ", err))
                           return(1)})
  data.frame(HomoD_tumor = total_vect[1], HemiD_tumor = total_vect[2], HomoD_rest = total_vect[3], HemiD_rest = total_vect[4], FBexact_pval = total_pval)
}

outdf <- do.call(rbind, by(data = hdhemidcounts, INDICES = hdhemidcounts$peakregion, FUN = get_type_pval))
pvals <- cbind(pvals[, !colnames(pvals) %in% c("HomoD_tumor", "HemiD_tumor", "HomoD_rest", "HemiD_rest")], outdf)

pvals$var1 <- log10((pvals$HDs) / (pvals$hemiDs))
pvals$var2 <- log10((pvals$nmajor2t) / (pvals$nmajor))
pvals$var3 <- log10((pvals$nHemilarge) / (pvals$nHemishort))

## redo reasoning:
# define variables in which the two groups differ: i.e. var1 - var2 - var3
# and test for this using a permutation test (e.g. oneway_test in the coin package,
# which performs a Fisher-Pitman permutation test) to show difference in the distribs.
# next use the resampled distrib of the FS group to compute p-values for all regions

pvals_sub <- pvals[pvals$status %in% c("FS", "TS"),]
var1 <- pvals_sub$var1
var2 <- pvals_sub$var2
var3 <- pvals_sub$var3
pvals_status <- factor(x = pvals_sub$status, levels = c("FS", "TS"))
oneway_test(var1pc ~ pvals_status, alternative = "less", distribution=exact())
oneway_test(var2pc ~ pvals_status, alternative = "less", distribution=exact())
oneway_test(var3pc ~ pvals_status, alternative = "less", distribution=exact())

# plot(lm(log10(HDs+.5) ~ log10(hemiDs+.5), data = pvals))

nsamples <- 1e6

sampleidxs <- c(sample(x = which(pvals$Chr != "X" & pvals$status == "FS"), size = nsamples, replace = T),
                sample(x = which(pvals$Chr != "X" & pvals$status == "TS"), size = nsamples, replace = T))

probmat <- pvals[sampleidxs, c("HDs", "hemiDs", "nHemilarge")]
probmat$none <- 2137 - rowSums(probmat)
resampledvars <- as.data.frame(t(apply(X = probmat, MARGIN = 1, FUN = function(x) rmultinom(n = 1, size = 2137, prob = x))))[,-4]
colnames(resampledvars) <- c("HDs", "hemiDs", "nHemilarge")
rownames(resampledvars) <- NULL
resampledvars$nHemishort <- resampledvars$hemiDs

resampledvars$nmajor2t <- rbinom(n = 2*nsamples, size = resampledvars$HDs, prob = (pvals[sampleidxs, c("nmajor2t")] + .5) / (rowSums(pvals[sampleidxs, c("nmajor2t", "nmajor")]) + 1) )
resampledvars$nmajor <- resampledvars$HDs - resampledvars$nmajor2t
resampledvars$status <- rep(c("FS", "TS"), c(nsamples, nsamples))



plot_simulated_distribution <- function(simuldata, dataset, model_variable, xlims = c(0.01,100), xbreaks = c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000), mmodel = NULL) {

  p1 <- ggplot(data = simuldata)
  p1 <- p1 + scale_x_continuous(breaks = log10(xbreaks), labels = as.character(xbreaks), limits = log10(xlims))
  p1 <- p1 + geom_density(mapping = aes(x = log10((simuldata[,1]+.5)/(simuldata[,2]+.5)), fill = status, colour = status), size = .05, alpha = .25, show.legend = FALSE)
  if (!is.null(mmodel))
    p1 <- p1 + stat_function(fun = function(x) .75*dnorm(x, mean = mmodel$mu[1], sd = mmodel$sigma[1])) + stat_function(fun = function(x) .75*dnorm(x, mean = mmodel$mu[2], sd = mmodel$sigma[2]))
  p1 <- p1 + geom_point(data = dataset[dataset$status != "unknown",], mapping = aes(x = dataset[dataset$status != "unknown",model_variable], y = 0, colour = status), shape = "|", size = 5, alpha = 0.6, position = position_jitter(w = 0.01, h = 0), show.legend = FALSE)
  p1 <- p1 + geom_point(data = dataset[dataset$status == "unknown",], mapping = aes(x = dataset[dataset$status == "unknown",model_variable], y = 0), colour = "grey40", shape = "|", size = 3, alpha = 0.6, position = position_jitter(w = 0.01, h = 0), show.legend = FALSE)
  p1 <- p1 + theme(axis.title = element_blank(), panel.grid = element_blank(), axis.text.y = element_blank(), axis.line.x = element_line(size = 0.25), axis.ticks.y = element_blank(), axis.ticks.x = element_line(size = 0.25), axis.ticks.length = unit(1, "mm"),
                   panel.background = element_rect(fill = "white"))
  return(p1)
}

p1 <- plot_simulated_distribution(simuldata = resampledvars[, c("HDs", "hemiDs", "status")], dataset = pvals[!(pvals$Chr == "X" | pvals$Gene %in% c("TCR", "IGH", "IGL")),], model_variable = "var1pc", xlims = c(0.01,100))
p1

p2 <- plot_simulated_distribution(simuldata = resampledvars[, c("nmajor2t", "nmajor", "status")], dataset = pvals[!(pvals$Chr == "X" | pvals$Gene %in% c("TCR", "IGH", "IGL")),], model_variable = "var2pc")
p2

p3 <- plot_simulated_distribution(simuldata = resampledvars[, c("nHemilarge", "nHemishort", "status")], dataset = pvals[!(pvals$Chr == "X" | pvals$Gene %in% c("TCR", "IGH", "IGL")),], model_variable = "var3pc", xlims = c(3,3000))
p3

ggsave("model1_binom_resample.pdf", plot = p1, width = 4, height = 4)
ggsave("model2_binom_resample.pdf", plot = p2, width = 4, height = 4)
ggsave("model3_binom_resample.pdf", plot = p3, width = 4, height = 4)




nsamples <- 1e7
sampleidxs <- sample(x = which(pvals$Chr != "X" & pvals$status == "FS"), size = nsamples, replace = T)

probmat <- pvals[sampleidxs, c("HDs", "hemiDs", "nHemilarge")]
probmat$none <- 2137 - rowSums(probmat)
resampledvars <- t(apply(X = probmat, MARGIN = 1, FUN = function(x) rmultinom(n = 1, size = 2137, prob = x)))[,-4]
dimnames(resampledvars) <- list(NULL, c("HDs", "hemiDs", "nHemilarge"))
nm2t <- rbinom(n = nsamples, size = resampledvars[,"HDs"], prob = (pvals[sampleidxs, c("nmajor2t")] + .5) / (rowSums(pvals[sampleidxs, c("nmajor2t", "nmajor")]) + 1) )
resampledvars <- cbind(resampledvars, nHemishort = resampledvars[,"hemiDs"], nmajor2t = nm2t, nmajor = resampledvars[,"HDs"] - nm2t)

sampledmat <- cbind(HH = (resampledvars[, "HDs"] + .5) / (resampledvars[, "hemiDs"] + .5),
                    len = (resampledvars[, "nmajor2t"] + .5) / (resampledvars[, "nmajor"] + .5),
                    Hemi = (resampledvars[, "nHemilarge"] + .5) / (resampledvars[, "nHemishort"] + .5))


binompvals_1_smpl <- apply(X = pvals[, c("HDs", "hemiDs")], MARGIN = 1, FUN = function(x, simvar) sum(simvar >= ((x[["HDs"]]) / (x[["hemiDs"]]))) / nsamples, simvar = sampledmat[, "HH"])
binompvals_2_smpl <- apply(X = pvals[, c("nmajor2t", "nmajor")], MARGIN = 1, FUN = function(x, simvar) sum(simvar >= ((x[["nmajor2t"]]) / (x[["nmajor"]]))) / nsamples, simvar = sampledmat[, "len"])
binompvals_3_smpl <- apply(X = pvals[, c("nHemilarge", "nHemishort")], MARGIN = 1, FUN = function(x, simvar) sum(simvar >= ((x[["nHemilarge"]]) / (x[["nHemishort"]]))) / nsamples, simvar = sampledmat[, "Hemi"])

pvals[, c("bpval_HH_smpl", "bpval_len_smpl", "bpval_Hemi_smpl")] <- cbind(binompvals_1_smpl, binompvals_2_smpl, binompvals_3_smpl)


empiricalBrownsWrapper <- function(datarow) {
  results <- empiricalBrownsMethod(data_matrix = t(pvals[,c("bpval_HH_smpl", "bpval_Hemi_smpl", "FBexact_pval")]), p_values = datarow, extra_info = TRUE)
  return(results)
}

pvals$bpval_Hemi_smpl <- ifelse(pvals$bpval_Hemi_smpl == 0, 1/nsamples, pvals$bpval_Hemi_smpl)
pvals$bpval_HH_smpl <- ifelse(pvals$bpval_HH_smpl == 0, 1/nsamples, pvals$bpval_HH_smpl)



brownspval <- apply(pvals[pvals$Chr != "X", c("bpval_HH_smpl", "bpval_Hemi_smpl", "FBexact_pval")], MARGIN = 1, FUN = empiricalBrownsWrapper)
pvals$Brown_pval <- NA
pvals[pvals$Chr != "X", "Brown_pval"] <- unlist(lapply(X = brownspval, FUN = function(x) x$P_test))
pvals$Brown_adj <- p.adjust(pvals$Brown_pval, method = "BH")


pairs(pvals[, c("comp.1_HH", "comp.1_len", "comp.2_Hemi", "Pvalue", "padjraw", "bpval_HH_smpl", "bpval_len_smpl", "bpval_Hemi_smpl", "FBexact_pval", "Brown_adj")])


FDR <- 0.05



### plotting p-val distrib

plot_empirical_pvalue_distribution <- function(dataset, model_variable, xlim = c(0,1), xbreaks = seq(0,1,0.2), bw = 0.05) {

  p1 <- ggplot(data = dataset, mapping = aes(dataset[,model_variable])) + theme_minimal()
  p1 <- p1 + scale_x_continuous(limits = xlim, breaks = xbreaks)
  # p1 <- p1 + scale_x_continuous(breaks = log10(xbreaks), labels = as.character(xbreaks), limits = log10(xlims))
  p1 <- p1 + geom_density(mapping = aes(y = ..density..), fill = "grey", size = 0.1, alpha = 0.2, bw = bw, trim = F)
  # p1 <- p1 + geom_density(mapping = aes(x = (dataset[,numerator]/dataset[,denominator]), y = ..count.., fill = is_TS), size = 0.1, alpha = 0.4, bw = 0.1)
  p1 <- p1 + geom_point(data = dataset[dataset$status != "unknown",], mapping = aes(x = ifelse(dataset[dataset$status != "unknown",model_variable] <= 0.05, dataset[dataset$status != "unknown",model_variable] + runif(nrow(dataset[dataset$status != "unknown",]), min = 0, max = 0.005),
                                                                                               ifelse(dataset[dataset$status != "unknown",model_variable] <= 0.95, dataset[dataset$status != "unknown",model_variable] + runif(nrow(dataset[dataset$status != "unknown",]), min = -0.012, max = 0.012),
                                                                                                      dataset[dataset$status != "unknown",model_variable] - runif(nrow(dataset[dataset$status != "unknown",]), min = 0, max = 0.025))), y = 0, colour = status), shape = "|", size = 5, alpha = 0.6, show.legend = FALSE)
  p1 <- p1 + geom_point(data = dataset[dataset$status == "unknown",], mapping = aes(x = ifelse(dataset[dataset$status == "unknown",model_variable] <= 0.05, dataset[dataset$status == "unknown",model_variable] + runif(nrow(dataset[dataset$status == "unknown",]), min = 0, max = 0.01),
                                                                                               ifelse(dataset[dataset$status == "unknown",model_variable] <= 0.95, dataset[dataset$status == "unknown",model_variable] + runif(nrow(dataset[dataset$status == "unknown",]), min = -0.012, max = 0.012),
                                                                                                      dataset[dataset$status == "unknown",model_variable] - runif(nrow(dataset[dataset$status == "unknown",]), min = 0, max = 0.025))), y = 0), colour = "grey40", shape = "|", size = 3, alpha = 0.6, show.legend = FALSE)
  p1 <- p1 + theme(axis.title = element_blank(), panel.grid = element_blank(), axis.text.y = element_blank(), axis.line.x = element_line(size = 0.25), axis.ticks.x = element_line(size = 0.25), axis.ticks.length = unit(1, "mm"))
  return(p1)
}

pl4 <- plot_empirical_pvalue_distribution(pvals[pvals$Chr != "X" & !pvals$Gene %in% c("TCR", "IGH", "IGL"),], "bpval_HH_smpl", bw = 0.08)
pl4
ggsave("pvalmodel1.pdf", plot = pl4, width = 4, height = 4)
pl5 <- plot_empirical_pvalue_distribution(pvals[pvals$Chr != "X" & !pvals$Gene %in% c("TCR", "IGH", "IGL"),], "bpval_len_smpl", bw = 0.08)
pl5
ggsave("pvalmodel2.pdf", plot = pl5, width = 4, height = 4)
pl6 <- plot_empirical_pvalue_distribution(pvals[pvals$Chr != "X" & !pvals$Gene %in% c("TCR", "IGH", "IGL"),], "bpval_Hemi_smpl", bw = 0.08)
pl6
ggsave("pvalmodel3.pdf", plot = pl6, width = 4, height = 4)
pl7 <- plot_empirical_pvalue_distribution(pvals[pvals$Chr != "X" & !pvals$Gene %in% c("TCR", "IGH", "IGL"),], "FBexact_pval", bw = 0.08)
pl7
ggsave("pvalmodel4.pdf", plot = pl7, width = 4, height = 4)


## final plot

sum(pvals[pvals$Chr != "X" & !pvals$Gene %in% c("TCR", "IGH", "IGL"), "Brown_adj"] < FDR)*FDR/nrow(pvals[pvals$Chr != "X" & !pvals$Gene %in% c("TCR", "IGH", "IGL"), ])
# pvals[pvals$Gene %in% c("MAP2K4", "CDH1"), "status"] <- "TS"
p1 <- ggplot(data = pvals[pvals$Chr != "X" & !pvals$Gene %in% c("TCR", "IGH", "IGL"), ], mapping = aes(Brown_pval)) + theme_minimal()
p1 <- p1 + scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2))
p1 <- p1 + geom_density(mapping = aes(y = ..density..), fill = "grey", size = 0.1, alpha = 0.2, bw = 0.08, trim = F)
p1 <- p1 + geom_vline(xintercept = sum(pvals[pvals$Chr != "X" & !pvals$Gene %in% c("TCR", "IGH", "IGL"), "Brown_adj"] < FDR)*FDR/nrow(pvals[pvals$Chr != "X" & !pvals$Gene %in% c("TCR", "IGH", "IGL"),]), linetype = "dashed")
# p1 <- p1 + geom_density(mapping = aes(x = unlist(P_Fisher), y = ..density..), fill = "blue", size = 0.1, alpha = 0.2, bw = 0.08, trim = F)
p1 <- p1 + geom_point(data = pvals[pvals$status == "unknown",], mapping = aes(x = ifelse(pvals[pvals$status == "unknown","Brown_pval"] <= 0.01,
                                                                                         pvals[pvals$status == "unknown","Brown_pval"] + runif(nrow(pvals[pvals$status == "unknown",]), min = 0, max = 0.01),
                                                                                         pvals[pvals$status == "unknown","Brown_pval"]),
                                                                              y = 0), colour = "grey40", shape = "|", size = 3, alpha = 0.6, show.legend = FALSE)
p1 <- p1 + geom_point(data = pvals[pvals$status != "unknown",], mapping = aes(x = ifelse(pvals[pvals$status != "unknown","Brown_pval"] <= 0.01,
                                                                                         pvals[pvals$status != "unknown","Brown_pval"] + runif(nrow(pvals[pvals$status != "unknown",]), min = 0, max = 0.01),
                                                                                         pvals[pvals$status != "unknown","Brown_pval"])
                                                                              , y = 0, colour = status), shape = "|", size = 5, alpha = 0.6, show.legend = FALSE)
p1 <- p1 + theme(axis.title = element_blank(), panel.grid = element_blank(), axis.text.y = element_blank(), axis.line.x = element_line(size = 0.25), axis.ticks.x = element_line(size = 0.25), axis.ticks.length = unit(1, "mm"))
p1

ggsave("pvalfinalmodel.pdf", plot = p1, width = 12.75, height = 2)
# pvals[pvals$Gene %in% c("MAP2K4", "CDH1"), "status"] <- "TS in FS"



## write out new tables
pvalsout <- pvals[,c("peakregion", "Chr", "Start", "Stop", "HDs", "Gene", "status", "Brown_pval", "Brown_adj")]
pvalsout[pvalsout$Chr == "X", "Brown_pval"] <- pvals[pvals$Chr == "X", "FBexact_pval"]
write.table(x = pvalsout, file = "20170609_allTS-FRAinfo.txt", quote = F, sep = "\t", row.names = F, col.names = T)
options(scipen = -3)
pvalsout$P_valFDR <- paste0(signif(pvalsout$Brown_pval, 3), " (", signif(pvalsout$Brown_adj, 3), ")")
pvalsout$Chr <- factor(pvalsout$Chr, levels = c(1:22, "X"))
pvalsout <- pvalsout[order(pvalsout$Chr),]

write.table(x = subset(pvalsout, status %in% c("TS"), c(peakregion, HDs, P_valFDR, Gene)), file = "20170609_STable5_a_knowTS.tsv", quote = F, sep = "\t", row.names = F, col.names = T)
immuneloci <- c("TCR", "IGH", "IGL")
write.table(x = subset(pvalsout, Gene %in% immuneloci, c(peakregion, HDs, P_valFDR, Gene)), file = "20170609_STable5_b_Immune.tsv", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = subset(pvalsout, Gene %in% c(known_FS, "CDH1", "MAP2K4"), c(peakregion, HDs, P_valFDR, Gene)), file = "20170609_STable5_c_knownFS.tsv", quote = F, sep = "\t", row.names = F, col.names = T)

telomeric <- c("chr2:0.92-0.94", "chr6:170.76-170.91", "chr7:158.91-159.13", "chr8:0.42-0.78",
               "chr9:0.76-0.88", "chr13:113.09-115.05", "chr17:80.94-81.01",
               "chr18:76.71-77.80", "chrX:0.10-1")
write.table(x = subset(pvalsout, (Brown_adj > FDR | is.na(Brown_adj)) & !peakregion %in% telomeric & !Gene %in% c(known_FS, known_TS, immuneloci), c(peakregion, HDs, P_valFDR, Gene)), file = "20170609_STable5_d_predFS.tsv", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = subset(pvalsout, peakregion %in% telomeric, c(peakregion, HDs, P_valFDR, Gene)), file = "20170609_STable5_e_telomeric.tsv", quote = F, sep = "\t", row.names = F, col.names = T)

write.table(x = subset(pvalsout, Brown_adj <= FDR & !Gene %in% known_TS & !Gene %in% immuneloci, c(peakregion, HDs, P_valFDR, Gene)), file = "20170609_STable5_f_newTS-unkn.tsv", quote = F, sep = "\t", row.names = F, col.names = T)

pvalsout2 <- pvalsout
pvalsout2$finstatus <- ifelse(pvalsout$status %in% c("FS", "TS"), pvalsout$status,
                              ifelse(pvals$Gene %in% immuneloci, "Immune",
                                     ifelse(pvalsout$peakregion %in% telomeric, "tel",
                                            ifelse(is.na(pvalsout$Brown_adj), "canFS",
                                                   ifelse(pvalsout$Brown_adj <= 0.05, "canTS", "canFS")))))
write.table(x = pvalsout2, file = "20170609_allTS-FRAinfo_FinalAnnot.txt", quote = F, sep = "\t", row.names = F, col.names = T)
options(scipen = 0)
