
call_Z_maxNegCtrl = function(Ymtx, Neg_control){
  Z_mtx = sapply(1:ncol(Ymtx),function(i){
    ifelse(Ymtx[,i]>max(Neg_control[,i]),1,0)
  })
  return(Z_mtx)
}

call_Z_normNegCtrl = function(Ymtx, Neg_control){
  Z_mtx = sapply(1:ncol(Ymtx),function(i){
    ifelse(Ymtx[,i]>max(Neg_control[,i]),1,0)
  })
  return(Z_mtx)
}
