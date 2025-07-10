library(rjson)

pErdosList = seq(0.15, 0.45, 0.05)

graphs = fromJSON(file="O10L5.json")
length(graphs)


table = matrix(0,length(pErdosList),4)
rowMatching = as.list(1:length(pErdosList))
names(rowMatching) = pErdosList

for (k in 1:length(graphs)){
  obj = graphs[[k]]
  row = rowMatching[[as.character(obj$pErdos)]]
  if (!any(obj$res3 == "NA") && !any(obj$resCan3 == "NA")){
    if (obj$res3$id){
      table[row, 1] <- table[row, 1]+1
    }
    if (obj$resCan3$id){
      table[row, 2] <- table[row, 2]+1
    }
    if (obj$res3$id && !obj$resCan3$id){
      table[row, 3] <- table[row, 3]+1
    }
    if (!obj$res3$id && obj$resCan3$id){
      table[row, 4] <- table[row, 4]+1
    }
  }
}

colnames(table) <- c("LSC", "CanLSC", "LSC\\CanLSC", "CanLSC\\ LSC")
rownames(table) <- pErdosList

print(table)
colSums(table)
write.table(table, file = "O10L5-can.txt", sep = "\t", row.names = TRUE, col.names = TRUE)