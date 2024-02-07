


s <- as.character(stats.tf$tf) %in% tfs_to_remove
head(s)
length(s)
which(s)
s[34]

stats.tf[34,]

tmp.df <- stats.tf[which(!(as.character(stats.tf$tf) %in% tfs_to_remove)), ]
dim(tmp.df)


