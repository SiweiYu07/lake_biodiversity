#### calculate alpha_fundiversity
  
install.packages("fundiversity")
getwd()
setwd("/Users/siweiyu/SW/FD")

library(fundiversity)
abundance<-read.csv("abundance.csv",row.names = 1)
abundance<-as.matrix(abundance)
trait<-read.csv("trait.csv",row.names = 1)

# FRic trait
fd_fric(trait)

#calculate abundance
fd_fric(trait,abundance)

#standardize the trait with abundance
fd_fric(trait,abundance,stand=T)

#calculate FRic_intersect index
fd_fric_intersect(trait,abundance,stand=T)

#calculate FDiv index
fd_fdiv(trait,abundance)

#calculate FEve index
fd_feve(trait,abundance)

#calculate FDis index
fd_fdis(trait)
fd_fdis(trait,abundance)

#calculate Rao's Q index

fd_raoq(trait)
fd_raoq(trait,abundance)

--------------------------------------------------------
