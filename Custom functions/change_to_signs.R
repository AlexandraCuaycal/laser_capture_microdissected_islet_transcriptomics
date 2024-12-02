## change to signs
change_to_signs <- function(string){
  
  for(i in 1:length(string)){
    word <- string[i]
    word <- gsub("\\-", "vs", word)
    word <- gsub("_pos", "\\+", word)
    word <- gsub("_neg", "\\-", word)
    #word <- gsub("_INS", "CD3", word)
    word <- gsub("_CD3", "CD3", word)
    string[i] <- word
    
  }
  return(string)
}