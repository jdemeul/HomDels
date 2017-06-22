## HDrecurrence

library(GenomicRanges)
library(doParallel)
library(biomaRt)

allhds <- read.delim(file = "SuppTable3.txt", as.is = T)

chrominfodf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/refdata/hg19.chrom.info", as.is = T, header = F)
chrominfo <- Seqinfo(seqnames = chrominfodf$V1, seqlengths = chrominfodf$V2, isCircular = rep(F, nrow(chrominfodf)), genome = "hg19")
seqlevels(chrominfo) <- c(1:22, "X")

hd_ranges <- GRanges(seqnames = allhds$chr, ranges = IRanges(start = allhds$start.position,
                                                                        end = allhds$end.position),
                        seqinfo = chrominfo,
                        mcols = allhds[ , c("sample", "cancer.type")])

hd_ranges <- sort(hd_ranges)


hd_ranges_disj <- disjoin(hd_ranges)
mcols(hd_ranges_disj) <- DataFrame(countOverlaps(hd_ranges_disj, hd_ranges))
observed_nohds <- mcols(hd_ranges_disj)[,1]


names(hd_ranges) <- NULL
hd_ranges_df <- as.data.frame(hd_ranges)


perform_iteration_opt2 <- function(hdranges_df, testranges, chrominfo, observed) {
  chrominfo_df <- as.data.frame(chrominfo)
  rand_chr <- sample(c(1:22,"X"), nrow(hdranges_df), replace = T)
  rand_start <- sapply(X = chrominfo_df[rand_chr, "seqlengths"] - hdranges_df$width, FUN = sample, size = 1)
  rand_hd_ranges <- GRanges(seqnames = rand_chr, IRanges(start = rand_start,
                                                         end = rand_start + hdranges_df$width - 1),
                            seqinfo = chrominfo)
  excess_observed <- countOverlaps(testranges, rand_hd_ranges) >= observed
  return(excess_observed)
}


NTHREADS <- 18
niter <- 1e4
clp <- makeCluster(NTHREADS)
registerDoParallel(clp)


for (i in 1:1000) {
print(paste0("iteration: ", i))
no_excesses <- rep(0, length(hd_ranges_disj))
no_excesses <- foreach(randrun=1:niter, .export = c("GRanges", "IRanges", "countOverlaps", "as.data.frame"), .combine = "+") %dopar% {
  perform_iteration_opt2(hdranges = hd_ranges_df,
                         testranges = hd_ranges_disj, chrominfo = chrominfo,
                         observed = observed_nohds)
}
write.table(x = no_excesses, file = paste0("no_excesses_batch", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = F)
}

p_obs <- no_excesses / niter
q_obs <- p.adjust(p_obs, method = "BH")

stopCluster(clp)



## Gene annotations
## Get all cancer gene census TSGs (== TSG | Rec)
ensus <- read.delim("/srv/data/vanloo/jdemeul/refdata/Census_allMon%20Feb%2013%2009-36-09%202017.tsv", as.is = T)
tsgs <- census[grepl("Rec", census$Molecular.Genetics) | grepl("TSG", census$Role.in.Cancer),]

## Get genes inside regions:
ensembl37 <- useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl37)

## significant ranges
mcols(subcl_ranges_disj) <- DataFrame(hdcount = observed_nohds, pval = p_obs, qval = q_obs)
subcl_ranges_disj_sign <- subcl_ranges_disj[q_obs <= 0.05]

subcl_ranges_disj_sign_merged <- reduce(subcl_ranges_disj_sign)
chromregions <- paste(seqnames(subcl_ranges_disj_sign_merged),
                      start(subcl_ranges_disj_sign_merged),
                      end(subcl_ranges_disj_sign_merged), sep = ":")

out <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name", "gene_biotype")
             , filters = c("chromosomal_region", "biotype"), values = list(chromregions, c("protein_coding")), mart = ensembl)
genes_ranges <- GRanges(seqnames = out$chromosome_name, IRanges(start = out$start_position,
                                                                end = out$end_position),
                        seqinfo = chrominfo, mcols = DataFrame(geneID = out$external_gene_name))
names(genes_ranges) <- out$external_gene_name


hittsgs <- out[out$external_gene_name %in% tsgs$Gene.Symbol, "external_gene_name"]
hittsgs_ranges <- genes_ranges[which(names(genes_ranges) %in% tsgs$Gene.Symbol)]

tsghits <- findOverlaps(hittsgs_ranges, subcl_ranges_disj)
overlappingranges <- by(data = subjectHits(tsghits), INDICES = queryHits(tsghits), FUN = c, simplify = T)
overlappingranges <- lapply(X = overlappingranges, FUN = function(x) mcols(subcl_ranges_disj)[x, ])
overlappingindices <- lapply(X = overlappingranges, FUN = function(x) which.max(x$hdcount))
maxoverlaps <- do.call(rbind, mapply(FUN = function(x, y) x[y, ], x = overlappingranges, y = overlappingindices, SIMPLIFY = T))
mcols(hittsgs_ranges) <- cbind(mcols(hittsgs_ranges), maxoverlaps)

## write files
hdrecurence_alldata <- data.frame(as.data.frame(subcl_ranges_disj))
colnames(hdrecurence_alldata)[c(1,6)] <- c("chr", "hdcount")
write.table(hdrecurence_alldata[ , -c(4,5)],
            file = "/srv/data/vanloo/jdemeul/ICGC/HDrecurrence/20170213_hdrecurrence.txt",
            quote = F, sep = "\t", row.names = F, col.names = T)

hdrecurrence_tsgs <- as.data.frame(hittsgs_ranges)[, c(1:3, 6:9)]
colnames(hdrecurrence_tsgs) <- c("chr", "start", "end", "gene", "hdcount_at_max", "pval", "qval")
write.table(hdrecurrence_tsgs,
            file = "/srv/data/vanloo/jdemeul/ICGC/HDrecurrence/20170213_hdrecurrence_tsgs.txt",
            quote = F, sep = "\t", row.names = F, col.names = T)
