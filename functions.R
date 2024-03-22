mNW <- function(x, X, Y, h, K = dnorm) {
  
  # Arguments
  # x: evaluation points
  # X: vector (size n) with the predictors
  # Y: vector (size n) with the response variable
  # h: bandwidth
  # K: kernel
  
  # Matrix of size n x length(x)
  Kx <- sapply(X, function(Xi) K((x - Xi) / h) / h)
  
  # Weights
  W <- Kx / rowSums(Kx) # Column recycling!
  
  # Means at x ("drop" to drop the matrix attributes)
  drop(W %*% Y)
  
}

scalar.prod <- function(tt,X,phi){
  # give an approximation of <X,phi> with the trapezoidal rule at time grid tt
  ## tt: vector of size p
  ## X: vector of size p
  ## phi: value of the function phi at the points tt
  p = length(tt)
  sum((tt[2:p]-tt[1:(p-1)])*(X[2:p]*phi[2:p]+X[1:(p-1)]*phi[1:(p-1)]))/2
}  

  
approx.histo <- function(tt,X,D,q=500){
  # give an approximation of X observed at time tt with a regular histogram basis
  ## tt: vector of size p
  ## X: vector of size p
  ## D: number of bins
  ##q: number of points of the final grid
  tildeX = rep(0,q)
  grid = seq(0,1,length.out=q)
  for (d in 1:D){
    phi = sqrt(D)*(((d-1)/D)<tt & (d/D)>=tt)
    tildeX = tildeX + scalar.prod(tt,X,phi)*sqrt(D)*(((d-1)/D)<grid & (d/D)>=grid)
  }
  tildeX
}

approx.fourier <- function(tt,X,D,q=500){
  # give an approximation of X observed at time tt with a regular histogram basis
  ## tt: vector of size p
  ## X: vector of size p
  ## D: number of functions (odd)
  ## q: number of points of the final grid
  grid = seq(0,1,length.out=q)
  phi1 = rep(1,length(tt))
  tildeX = rep(scalar.prod(tt,X,phi1),q)
  
  if(D>1){
    Dp = (D-1)/2
    for (d in 1:Dp){
      phi2d = sqrt(2)*cos(2*pi*d*tt)
      tildeX = tildeX + scalar.prod(tt,X,phi2d)*sqrt(2)*cos(2*pi*d*grid)
      phi2dp = sqrt(2)*sin(2*pi*d*tt)
      tildeX = tildeX + scalar.prod(tt,X,phi2dp)*sqrt(2)*sin(2*pi*d*grid)
    }
  }
  tildeX
}
