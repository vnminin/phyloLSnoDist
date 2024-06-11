##################################################
#  R code written by:                            #
#                                                #
#  Jimmy A. Doi (jdoi@calpoly.edu)               #
#  Department of Statistics                      #
#  Cal Poly State Univ, San Luis Obispo          #
#  Web: www.calpoly.edu/~jdoi                    #
#                                                #
#  ............................................  #
#                                                #
#  If using please cite:                         #
#                                                #
#  Schilling, M., Doi, J.                        #
#  "A Coverage Probability Approach to Finding   #
#  an Optimal Binomial Confidence Procedure",    #
#  The American Statistician, 68, 133-145.       #
#                                                #
#  ............................................  #
#                                                #
#  Shiny app site:     jdoi.shinyapps.io/LCO-CI  #
#                                                #
#  ............................................  #
#                                                #
#                     Code updated on: 1AUG2014  #
##################################################


##############################################################################
#    The function LCO.CI() generates the LCO confidence intervals            #
#    for X = 0, 1, ..., n  for a specified n and confidence level.           #
#                                                                            #
#    Example: To generate all LCO confidence intervals at n=20,              #
#             level=.90, and 3rd decimal place accuracy, use                 #
#                                                                            #
#    > LCO.CI(20,.90,3)                                                      #
##############################################################################


LCO.CI <- function(n,level,dp)
{

  # For desired decimal place accuracy of dp, search on grid using (dp+1) 
  # accuracy then round final results to dp accuracy.
  iter <- 10**(dp+1)
  
  p <- seq(0,.5,1/iter)
  
  
  ############################################################################
  # Create initial cpf with AC[l,u] endpoints by choosing coverage 
  # probability from highest acceptance curve with minimal span.
  
  
  cpf.matrix <- matrix(NA,ncol=3,nrow=iter+1)
  colnames(cpf.matrix)<-c("p","low","upp")
  
  for (i in 1:(iter/2+1)){
    p <- (i-1)/iter
    
    bin <- dbinom(0:n,n,p)
    x   <- 0:n
    pmf <- cbind(x,bin)

    # Binomial probabilities ordered in descending sequence
    pmf <- pmf[order(-pmf[,2],pmf[,1]),] 
    pmf <- data.frame(pmf)
    
    # Select the endpoints (l,u) such that AC[l,u] will
    # be at least equal to LEVEL. The cumulative sum of
    # the ordered pmf will identify when this occurs.
    m.row  <- min(which((cumsum(pmf[,2])>=level)==TRUE))
    low.val <-min(pmf[1:m.row,][,1])
    upp.val <-max(pmf[1:m.row,][,1])

    cpf.matrix[i,] <- c(p,low.val,upp.val)
    
    # cpf flip only for p != 0.5
    
    if (i != iter/2+1){
      n.p <- 1-p
      n.low <- n-upp.val
      n.upp <- n-low.val
    
      cpf.matrix[iter+2-i,] <- c(n.p,n.low,n.upp)
    }
  }
    

  ############################################################################
  # LCO Gap Fix
  # If the previous step yields any violations in monotonicity in l for 
  # AC[l,u], this will cause a CI gap. Apply fix as described in Step 2 of 
  # algorithm as described in paper.

  # For p < 0.5, monotonicity violation in l can be determined by using the 
  # lagged difference in l. If the lagged difference is -1 a violation has 
  # occurred. The NEXT lagged difference of +1 identifies the (l,u) pair to 
  # substitute with. The range of p in violation would be from the lagged 
  # difference of -1 to the point just before the NEXT lagged difference of 
  # +1. Substitute the old (l,u) with updated (l,u) pair. Then make required
  # (l,u) substitutions for corresponding p > 0.5. 
  
  # Note the initial difference is defined as 99 simply as a place holder.
  
  diff.l <- c(99,diff(cpf.matrix[,2],differences = 1))
  
  if (min(diff.l)==-1){
  
    for (i in which(diff.l==-1)){
      j <- min(which(diff.l==1)[which(diff.l==1)>i])
      new.low <- cpf.matrix[j,2]
      new.upp <- cpf.matrix[j,3]
      cpf.matrix[i:(j-1),2] <- new.low
      cpf.matrix[i:(j-1),3] <- new.upp
      }
   
    # cpf flip
    pointer.1 <- iter - (j - 1) + 2
    pointer.2 <- iter - i + 2
    
    cpf.matrix[pointer.1:pointer.2,2]<- n - new.upp
    cpf.matrix[pointer.1:pointer.2,3]<- n - new.low
  }
  
    
  ############################################################################
  # LCO CI Generation
  
  ci.matrix <-  matrix(NA,ncol=3,nrow=n+1) 
  rownames(ci.matrix) <- c(rep("",nrow(ci.matrix)))
  colnames(ci.matrix) <- c("x","lower","upper")
  
  # n%%2 is n mod 2: if n%%2 == 1 then n is odd
  # n%/%2 is the integer part of the division: 5/2 = 2.5, so 5%/%2 = 2
  
  if (n%%2==1) x.limit <- n%/%2
  if (n%%2==0) x.limit <- n/2
  
  for (x in 0:x.limit)
  {
    num.row <- nrow(cpf.matrix[(cpf.matrix[,2]<=x & x<=cpf.matrix[,3]),])

    low.lim <- 
      round(cpf.matrix[(cpf.matrix[,2]<=x & x<=cpf.matrix[,3]),][1,1],
            digits=dp)

    upp.lim <- 
      round(cpf.matrix[(cpf.matrix[,2]<=x & x<=cpf.matrix[,3]),][num.row,1],
            digits=dp)

    ci.matrix[x+1,]<-c(x,low.lim,upp.lim)
    
    # Apply equivariance
    n.x <- n-x
    n.low.lim <- 1 - upp.lim
    n.upp.lim <- 1 - low.lim
    
    ci.matrix[n.x+1,]<-c(n.x,n.low.lim,n.upp.lim)
  }
  
  
  heading <- matrix(NA,ncol=1,nrow=1)    
  
  heading[1,1] <- 
    paste("LCO Confidence Intervals for n = ",n," and Level = ",level,sep="")
  
  rownames(heading) <- c("")
  colnames(heading) <- c("")
  
  print(heading,quote=FALSE)
  
  print(ci.matrix)
}


##############################################################################
#    The function LCO.CI() generates the LCO confidence intervals            #
#    for X = 0, 1, ..., n  for a specified n and confidence level.           #
#                                                                            #
#    Example: To generate all LCO confidence intervals at n=20,              #
#             level=.90, and 3rd decimal place accuracy, use                 #
#                                                                            #
#    > LCO.CI(20,.90,3)                                                      #
##############################################################################

