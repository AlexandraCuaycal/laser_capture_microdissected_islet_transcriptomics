## changing histological phenotype character + and - to words plus and minus for piano
change_to_words <- function(string){
  
  for(i in 1:length(string)){
    word <- string[i]
    word <- gsub("\\+", "_pos", word)
    word <- gsub("\\-", "_neg", word)
    word <- gsub("CD3", "_CD3", word)
    string[i] <- word
    
  }
  return(string)
}