

#########  group_by_data  #################

####  2+2 +3 set up

theta <- cbind(c(0.2, 0.15, 0.1, 0.2, 0.10, 0.2, 0.05, 0),
               c(0.5, 0.2, rep(0,5), 0.3),
               c(0, rep(0.25, 4), 0, 0, 0)
)

signatures <- rbind.data.frame(c(0,0,0), c(0, 1, 0), c(1, 0, 0), c(1, 1, 0),
                         c(0,0,1), c(0, 1, 1), c(1, 0, 1), c(1, 1, 1))

new_theta <- matrix(0, dim(theta)[1], dim(theta)[2])
sig_list <- list()
for(k in 1:dim(theta)[2]){
  num_unique_sigs <- list()
  for(l in 1:dim(signatures)[2]){
    sig_list[[l]] <- tapply(theta[,k], factor(signatures[,l], levels=unique(signatures[,l])), sum)
    num_unique_sigs[[l]] <- 0:(length(unique(signatures[,l]))-1)
    if(l==1){
      f_array <- sig_list[[l]]
    }else{
      f_array <- outer(f_array, sig_list[[l]])
    }
  }
  signature_new <- as.numeric();
  vec <- numeric()
  grid <- expand.grid(num_unique_sigs)
  vec <- apply(grid, 1, function(x) return(f_array[matrix(as.numeric(x)+1,1)]))
  match(grid, signatures)
  vec <- vec[match(data.frame(t(signatures)), data.frame(t(grid)))]
  new_theta[,k] <-  vec
 
}



