##Parameters
quantsize=35000
#quantsize determines number of simulation runs
U=20
#we consider the uniform measure on 1/U, 2/U,...,(U-1)/U
brownupper=2000
#Number of L2-coefficients for the Brownian motion

##Program
unipoints=U-1
brownarray=array(0, dim=quantsize)
sinmat=array(0, dim=c(unipoints, brownupper))
for(i in 1:unipoints){
  for(k in 2:brownupper){
    sinmat[i,k]=sqrt(2)*sin(pi*(k-1)*i/U)/((k-1)*pi)
  }
}

quantfunction <- function(n){
  gaussvector=rnorm(brownupper, mean=0, sd=1)
  numerator=gaussvector[1]
  denominator=0
  for(i in 1:unipoints){
    denominator=denominator+sum((gaussvector*sinmat[i,]))^2*i^2/U^2
  }
  denominator=(denominator/unipoints)^(1/2)
  return(numerator/denominator)
}





simulate_w_quantiles <- function(alpha, n=1, plotHist=FALSE){
  
  for(i in 1:n){
    helparray=parallel::mclapply(brownarray, quantfunction, mc.cores=getOption("mc.cores", 6))
    for(quant in 1:quantsize){
      brownarray[quant]=helparray[[quant]]
    }
    
    ##Results
    if(plotHist){
      hist(brownarray, breaks=200) 
    }
    sortbrown =sort(brownarray, decreasing = FALSE)
    
    quants <- c()
    
    for(a in alpha){
      qnt =a*quantsize
      # print(paste(a, "-Quantil: ", sortbrown[qnt]))
      quants <- append(quants, sortbrown[qnt])    
    }
    
    return(quants)
  }
  # upperquantnn=(1-0.01)*quantsize
  # upperquantnf=(1-0.05)*quantsize
  # upperquantn=(1-0.1)*quantsize
  # print(sortbrown[upperquantnn])
  # print(sortbrown[upperquantnf])
  # print(sortbrown[upperquantn])
}




