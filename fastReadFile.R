faster.readfile <- function(x,y, ... ){
  require(data.table)
  a <- fread(x,colClasses = NULL,skip = y, fill = TRUE, stringsAsFactors = FALSE, header = TRUE)
  return(a)
}

