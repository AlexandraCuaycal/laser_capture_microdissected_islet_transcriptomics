## change to groups
change_to_groups <- function(string){
  
  for(i in 1:length(string)){
    word <- string[i]
    word <- gsub("\\-", "vs", word)
    word <- gsub("_pos", "", word)
    word <- gsub("_neg", "", word)
    word <- gsub("_INS", "", word)
    word <- gsub("_CD3", "", word)
    string[i] <- word
    
  }
  return(string)
}
