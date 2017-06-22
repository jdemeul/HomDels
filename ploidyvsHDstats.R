### Response to Reviewer 3 question 1: relation ploidy - HD
library(ggplot2)
library(Exact)
library(coin)
# library(plyr)

ploidytable <- read.delim("SuppTable2.txt")
ploidytable$ploidyclass <- ifelse(ploidytable$Ploidy..n. > 2.7, "tetraploid", "diploid")
hdtable <- read.delim("SuppTable3.txt")

p1 <- ggplot(data = ploidytable) + theme_minimal()
p1 <- p1 + geom_point(mapping = aes(x = Ploidy..n., y = X..HDs, colour = ploidyclass), size = 0.25, position = position_jitter(width = 0.1, height = 0.5), alpha = 0.45, show.legend = F)
# p1 <- p1 + geom_point(mapping = aes(x = Ploidy..n., y = log10(Total.length.of.HDs)), shape = ".", position = position_jitter(width = 0.025, height = 0.025), alpha = 0.75)
p1 <- p1 + geom_vline(xintercept = 2.7, linetype = "dashed")
p1 <- p1 + labs(x = "Tumor ploidy", y = "Total number of HDs", size = 0.5)
p1
ggsave("q1_ploidy_vs_noHDs.pdf", plot = p1, width = 4, height = 4)

p1 <- ggplot(data = ploidytable[ploidytable$Total.length.of.HDs > 1,]) + theme_minimal()
# p1 <- p1 + geom_point(mapping = aes(x = Ploidy..n., y = X..HDs), size = 0.25, position = position_jitter(width = 0.1, height = 0.5), alpha = 0.25)
p1 <- p1 + geom_point(mapping = aes(x = Ploidy..n., y = Total.length.of.HDs, colour = ploidyclass), size = 0.25, position = position_jitter(width = 0.0, height = 0), alpha = 0.70, show.legend = F)
p1 <- p1 + labs(x = "Tumor ploidy", y = "Total length of HDs (Mb)") + scale_y_log10(breaks = c(1e+3, 1e+4, 1e+5, 1e+6, 1e+7, 1e+8), labels = c("0.001", "0.01", "0.1", "1", "10", "100"))
p1 <- p1 + geom_vline(xintercept = 2.7, linetype = "dashed", size = 0.5)
p1
ggsave("q1_ploidy_vs_totHDlength.pdf", plot = p1, width = 4, height = 4)

# hdtable$ploidyclass <- ploidytable[match(hdtable$sample , ploidytable$Sample.name), "ploidyclass"]
hdtable_aug <- merge(hdtable, ploidytable[, c("Sample.name", "Ploidy..n.", "ploidyclass")], by.x = "sample", by.y = "Sample.name")
hdtable_aug$size <- hdtable_aug$end.position - hdtable_aug$start.position

p1 <- ggplot(data = hdtable_aug[hdtable_aug$size > 1, ]) + theme_minimal()
# p1 <- p1 + geom_point(mapping = aes(x = Ploidy..n., y = X..HDs), size = 0.25, position = position_jitter(width = 0.1, height = 0.5), alpha = 0.25)
p1 <- p1 + geom_point(mapping = aes(x = Ploidy..n., y = size, colour = ploidyclass), size = 0.25, position = position_jitter(width = 0.0, height = 0), alpha = 0.70, show.legend = F)
p1 <- p1 + labs(x = "Tumor ploidy", y = "Total length of HDs (Mb)") + scale_y_log10(breaks = c(1e+3, 1e+4, 1e+5, 1e+6, 1e+7, 1e+8), labels = c("0.001", "0.01", "0.1", "1", "10", "100"))
p1 <- p1 + geom_vline(xintercept = 2.7, linetype = "dashed", size = 0.5)
p1
ggsave("q1_ploidy_vs_totHDlength_redo.pdf", plot = p1, width = 4, height = 4)


p1 <- ggplot(data = ploidytable) + geom_bar(mapping = aes(x = X..HDs, y = ..prop.., fill = ploidyclass), position = "dodge", show.legend = F)
p1 <- p1 + theme_minimal() + coord_cartesian(xlim = c(0,20)) + labs(x = "Total number of HDs", y = "Proportion of samples")
ggsave("q2_noHDsbyPloidystate.pdf", plot = p1, width = 4, height = 4)


p1 <- ggplot(data = ploidytable[ploidytable$Total.length.of.HDs > 1,]) + geom_density(mapping = aes(x = Total.length.of.HDs, y = ..density.., fill = ploidyclass, show.legend = F), alpha = 0.5, show.legend = F)
p1 <- p1 + theme_minimal() + scale_x_log10(breaks = c(1e+3, 1e+4, 1e+5, 1e+6, 1e+7, 1e+8), labels = c("0.001", "0.01", "0.1", "1", "10", "100")) + labs(x = "Total length of HDs (Mb)", y = "Proportion of samples")
ggsave("q2_totalHDlengthbyPloidystate.pdf", plot = p1, width = 4, height = 4)

p1 <- ggplot(data = hdtable_aug[hdtable_aug$size > 1, ]) + geom_density(mapping = aes(x = size, y = ..density.., fill = ploidyclass, show.legend = F), alpha = 0.5, show.legend = F)
p1 <- p1 + theme_minimal() + scale_x_log10(breaks = c(1e+3, 1e+4, 1e+5, 1e+6, 1e+7, 1e+8), labels = c("0.001", "0.01", "0.1", "1", "10", "100")) + labs(x = "Length of HD (Mb)", y = "Proportion of HDs")
p1
ggsave("q2_totalHDlengthbyPloidystate_redo.pdf", plot = p1, width = 4, height = 4)


plclass <- factor(ploidytable[!is.na(ploidytable$ploidyclass), "ploidyclass"])
nohds <- ploidytable[!is.na(ploidytable$ploidyclass), "X..HDs"]
oneway_test(nohds ~ plclass, alternative = "greater", distribution=approximate(B = 1e6))

hdsize <- hdtable_aug[!is.na(hdtable_aug$ploidyclass), "size"]
plclass <- factor(hdtable_aug[!is.na(hdtable_aug$ploidyclass), "ploidyclass"])
wilcox_test(hdsize ~ plclass, alternative = "two.sided")




hdtable <- read.delim("SuppTable3.txt")
knownts <- read.delim("mart_export.txt", as.is = T)


library(GenomicRanges)
HDranges <- GRanges(seqnames = Rle(hdtable$chr), ranges = IRanges(start = hdtable$start.position,
                                                                  end = hdtable$end.position))
knowntsranges <- GRanges(seqnames = Rle(knownts$Chromosome.Name), ranges = IRanges(start = knownts$Gene.Start..bp.,
                                                                       end = knownts$Gene.End..bp.))
knowntsranges_ext <- knowntsranges
start(knowntsranges_ext) <- start(knowntsranges) - 1000000
end(knowntsranges_ext) <- end(knowntsranges) + 1000000

hits <- findOverlaps(HDranges, knowntsranges)

hdtable[queryHits(hits), "ts"] <- knownts[subjectHits(hits),"HGNC.symbol"]


### subquestion 2 - TS hit vs ploidy
hdtable_hitonly <- hdtable[!is.na(hdtable$ts),]
hdtable_hitonly_aug <- merge(hdtable_hitonly, ploidytable, by.x = "sample", by.y = "Sample.name")

tsfreq_ploidy <- as.data.frame(table(hdtable_hitonly_aug$ploidyclass, hdtable_hitonly_aug$ts))
by(data = tsfreq_ploidy$Freq, INDICES = tsfreq_ploidy$Var1, sum)
table(ploidytable$ploidyclass)
tsfreq_ploidy$Var2 <- as.vector(tsfreq_ploidy$Var2)
tsfreq_ploidy <- rbind(tsfreq_ploidy, data.frame(Var1 = "diploid", Var2 = "any", Freq = 143), data.frame(Var1 = "tetraploid", Var2 = "any", Freq = 59))
tsfreq_ploidy$Freqn <- ifelse(tsfreq_ploidy$Var1 == "diploid", tsfreq_ploidy$Freq/1501, tsfreq_ploidy$Freq/636)

pl1 <- ggplot(tsfreq_ploidy) + geom_bar(mapping = aes(x = Var2, y = Freqn*100, fill = Var1), stat = "identity", position = "dodge", show.legend = F)
pl1 <- pl1 + theme_minimal() + theme(axis.text.x = element_text(angle = 90))
pl1 <- pl1 + labs(x = "Tumor suppressor gene", y = "HD frequency per 100 cases")
pl1
ggsave("q3_TShitsbyPloidystate.pdf", plot = pl1, width = 4, height = 4)


get_TS_pval <- function(df, TS, totaldip, totaltet) {
  hddiptet <- df[df$Var2 == TS, "Freq"]
  mat <- matrix(data = c(hddiptet, c(totaldip, totaltet) - hddiptet), nrow = 2, ncol = 2, byrow = T,
                dimnames = list(c(TS, "other") , c("Diploid", "Tetraploid")))
  tryCatch(exact.test(data = t(mat), alternative = "two.sided", interval = T, method = "Boschloo", model = "Binomial", to.plot = F)$p.value,
                         error = function(err) {
                           print(paste("ERROR: ", err))
                           return(1)})
}


tsg_tests <- list()
for (tsg in unique(tsfreq_ploidy$Var2)) {
    tsg_tests[[tsg]] <- get_TS_pval(tsfreq_ploidy, tsg, 1501, 636)
}

sapply(tsg_tests, FUN = function(x) x$p.value)

table(hdtable_hitonly_aug[hdtable_hitonly_aug$ts == "RB1", c("ploidyclass", "cancer.type")])
table(ploidytable$Tumor.type, ploidytable$ploidyclass)
rb1mat <- matrix(data = c(16, 2, 63, 163), nrow = 2, byrow = F, dimnames = list(c("Tetraploid", "Diploid"), c("RB1-HD", "RB1-noHD")))
exact.test(data = rb1mat, alternative = "two.sided", interval = T, method = "Boschloo", model = "Binomial", to.plot = F)
