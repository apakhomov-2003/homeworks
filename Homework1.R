motifs2 <- matrix(c(
  "a", "C", "g", "G", "T", "A", "A", "t", "t", "C", "a", "G",
  "t", "G", "G", "G", "C", "A", "A", "T", "t", "C", "C", "a",
  "A", "C", "G", "t", "t", "A", "A", "t", "t", "C", "G", "G",
  "T", "G", "C", "G", "G", "G", "A", "t", "t", "C", "C", "C",
  "t", "C", "G", "a", "A", "A", "A", "t", "t", "C", "a", "G",
  "A", "C", "G", "G", "C", "G", "A", "a", "t", "T", "C", "C",
  "T", "C", "G", "t", "G", "A", "A", "t", "t", "a", "C", "G",
  "t", "C", "G", "G", "G", "A", "A", "t", "t", "C", "a", "C",
  "A", "G", "G", "G", "T", "A", "A", "t", "t", "C", "C", "G",
  "t", "C", "G", "G", "A", "A", "A", "a", "t", "C", "a", "C"
), nrow = 10, byrow = TRUE)

motifs2_upper <- apply(motifs2, 2, toupper)

count_motif2_matrix <- apply(motifs2_upper, 2, function(col) table
                             (factor(col, levels=c("A", "C", "G", "T"))))

profile_motif2_matrix <- apply(motifs2_upper, 2, function(x){
  counts <- table(factor(x, levels=c("A", "C", "G", "T")))
  counts/sum(counts)
})

scoreMotifs <- sum(apply(motifs2_upper, 2, function(x) length(x) - max(table(x))))
print(scoreMotifs)

getConsensus <- function(motifs) {
  consensus_vector <- apply(motifs, 2, function(x) {
    c("A", "C", "G", "T")[which.max(table(factor(x, levels=c("A", "C", "G", "T"))))]
    })
    paste(consensus_vector, collapse="")
}

barplot(prop.table(table(factor(motifs2_upper[,5], levels= c("A", "C", "G", "T")))),
        col = "skyblue",
        main = "Частота нуклеотидов в 5-ом столбце",
        xlab = "Нуклеотиды",
        ylab = "Доля",
        ylim = c(0,1))