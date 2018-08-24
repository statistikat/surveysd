myfun <- function(x,w){
  return(sum(w*x))
}
myfun.char <- function(x,w){
  return(as.character(sum(w*x)))
}
myfun.mulval <- function(x,w){
  return(w*x)
}