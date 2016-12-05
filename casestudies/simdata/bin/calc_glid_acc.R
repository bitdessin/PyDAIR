#


read_files <- function(gene) {
    tmpl <- matrix(NA, ncol = 6, nrow = 6)
    colnames(tmpl) <- paste0("J", 1:6)
    rownames(tmpl) <- paste0("V", 1:6)
    acc <- fdr <- tmpl
    for (vi in 1:6) {
        for (ji in 1:6) {
            f <- paste0("estperformance.sim", vi - 1, "_", ji - 1, ".", gene, ".csv")
            x <- read.table(f, fill = TRUE, stringsAsFactors = FALSE, sep = "\t")
            s <- apply(x[, 2:6], 2, sum, na.rm = TRUE)
            acc[vi, ji] <- s[4] / s[1]
            fdr[vi, ji] <- s[5] / s[3]
        }
    }
    return (list(acc = acc, fdr = fdr))
}


v <- read_files('v')
j <- read_files('j')

write.table(v$acc, file = "estperformance.all.acc.v.csv", sep = "\t", quote = FALSE)
write.table(j$acc, file = "estperformance.all.acc.j.csv", sep = "\t", quote = FALSE)
write.table(v$fdr, file = "estperformance.all.fdr.v.csv", sep = "\t", quote = FALSE)
write.table(j$fdr, file = "estperformance.all.fdr.j.csv", sep = "\t", quote = FALSE)





