###########################################################################################################################
### This function constructs the O'Sullivan penalized splines design matrix (nonlinear term)                            ###
### based on the paper "On semiparametric regression with O'Sullivan penalized splines" (2008) by Wand and Ormerod      ###
### This function is called by bin_Bayes_penalized_splines, not directly called by the user                             ###
### Parameters:                                                                                                         ###
### data = a data frame with the observations arranged by row, and including the columns below (in the following order):###
###        1) y: the response outcome (=1 has response; = 0 no response) (int)                                          ###
###        2) x: the value of the continuous marker (double)                                                            ###
###        3) trt: the group indicator (=1 for experimental group; =0 for control group) (int)                          ###
### xrange = a vector of length 2 which gives the range of marker values (vector)                                       ###
### numIntKnots =  the number of interior knots (should scale with number of events) (int)                              ###
### next_x =  the marker value for the about to be assigned patients (double)                                           ###
###########################################################################################################################



Get_splines <- function(data,xrange,numIntKnots,next_x = NA)
{
  
  library(splines)
  
  # Create the Omega matrix
  formOmega <- function(a,b,intKnots)
  {
    allKnots <- c(rep(a,4),intKnots,rep(b,4))
    K <- length(intKnots) 
    L <- 3*(K+8)
    xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+rep(allKnots,each=3)[-c(1,2,L)])/2
    wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
    Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),outer.ok=TRUE)$design 
    Omega <- t(Bdd*wts)%*%Bdd
    # Omega <- (t(Bdd)%*%diag_w)%*%Bdd is equivalent
    return(Omega)
  }
  
  if(!is.na(next_x)){

    x <- c(data$x,next_x)
    data <- data.frame(cbind(y = c(data$y,NA),x = x, trt = c(data$trt,1)))

  } else{

  x <- data$x
    
  }

  
  # Set up the design matrix and related quantities
  a <- xrange[1]
  b <- xrange[2]
  plot_seq <- seq(a,b,length.out = 500)
  
  intKnots <- quantile(unique(x),seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))])
  
  names(intKnots) <- NULL
  #cubic B-splines
  B <- bs(x,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE)
  
  #cubic B-splines for plotting
  B_plot <-bs(plot_seq,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE)
  
  Omega <- formOmega(a,b,intKnots)
  
  
  # Obtain the spectral decomposition of Omega
  eigOmega <- eigen(Omega)
  
  
  # Obtain the matrix for linear transformation of $\\bB$ to $\\bZ$
  indsZ <- 1:(numIntKnots+2)
  UZ <- eigOmega$vectors[,indsZ]
  LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))
  
  
  # Perform a stability check
  indsX <- (numIntKnots+3):(numIntKnots+4)
  UX <- eigOmega$vectors[,indsX]   
  L <- cbind(UX,LZ)
  stabCheck <- t(crossprod(L,t(crossprod(L,Omega))))  
  stabCheck_idx <- sum(stabCheck^2) > 1.0001*(numIntKnots+2)
  
  
  # Form Z matrices:
  Z <- B%*%LZ
  
  # Form Z matrices for plotting
  Z_plot <-B_plot%*%LZ
  
  # Interaction between trt and Z matrices:
  trtZ <- data$trt*Z
  
  #intercept vector
  #intercept = rep(1, nrow(data))
  
  # Add both matrices to the data frame:
  # outcome, intercept(1), x, Z matrix, trt, interaction btw x and trt, trtZ
  data2 <- as.data.frame(cbind(data[,1],  data[,2], Z, data$trt, 
    (data$x)*(data$trt), trtZ))
  colnames(data2) <- c(colnames(data)[1], colnames(data)[2],
    as.vector(outer("z", 1:ncol(Z), paste0)), colnames(data)[3], 
    "xtrt", paste0(as.vector(outer("z", 1:ncol(Z), paste0)),"trt") )
  
  formula <- as.formula(paste("y ~ ", paste(colnames(data2)[-c(1)], collapse="+") )) 
  
  #for plotting only need one design matrix:
  data_plot <- as.data.frame(cbind(rep(1,length(plot_seq)),  plot_seq, Z_plot))
  colnames(data_plot) <- c(colnames(data)[1], colnames(data)[2],
    as.vector(outer("z", 1:ncol(Z), paste0)))
  
    return(list(data=data2, formula=formula, stabCheck=stabCheck_idx, 
      ncolz=ncol(Z), next_pt = NA,data_plot = data_plot, B_matrix = B, LZ_matrix = LZ ))
    

}








