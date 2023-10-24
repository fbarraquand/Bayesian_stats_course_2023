### 2021-10-19 F Barraquand -- Bayesian Statistics Course for PhD Students

################# TD1 - estimation of a proportion #######################

### 1. Plotting the three prior distributions

# First param set = Uniform distribution
a1=1
b1=1

# Second param set = 10 times less obs than in the 18th century
#n=4848
#y=530
n=485
y=53
a2=y+1
b2=n-y+1
# we could have used the original, but it is a bit caricatural 
# (I initially discarded it because I computed the formula analytically and had difficulties with very large factorials)

# Third param set
mean_prop = y / n
kappa = 15
a3 = kappa*mean_prop
b3 = kappa-a3
# or
#a3 = 1.63
#b3 = 13.36

pdf(width=8,height=6,file="beta_priors.pdf")
curve(dbeta(x,a1,b1),from=0,to=0.6,lwd=3,ylim=c(0,30),ylab="Prior density",xlab = "Proportion (theta)")
curve(dbeta(x,a2,b2),from=0,to=0.6,col="orange",lwd=3,ylim=c(0,15),add=T)
curve(dbeta(x,530+1,4848-530+1),from=0,to=0.6,col="red",lwd=3,add=T)
abline(v=mean_prop,lwd=2,col="red", lty=2) # Previous estimate
curve(dbeta(x,a3,b3),from=0,to=0.6,col="violet",lwd=3,add=T)
dev.off()



### 2. Computing the corresponding posteriors

# New observed data from present year 
nnew = 79
ynew = 5

pdf(width=8,height=6,file="beta_posteriors.pdf")
curve(dbeta(x,ynew+a1,nnew-ynew+b1),from=0,to=0.6,lwd=3,ylim=c(0,30),ylab="Posterior density",xlab = "Proportion (theta)")
curve(dbeta(x,y+1+ynew,n-y+1+nnew-ynew),from=0,to=0.6,col="orange",lwd=3,add=T)
curve(dbeta(x,ynew+a3,nnew-ynew+b3),from=0,to=0.6,col="violet",lwd=3,add=T)
abline(v=mean_prop,lwd=2,col="red")
abline(v=ynew/nnew,lwd=2,col="blue")
dev.off()

### 3. Now computing the CRedibility Intervals at 95%
?qbeta
qbeta(c(0.025,0.975),ynew+a1,nnew-ynew+b1) # broader because of uniform prior
qbeta(c(0.025,0.975),ynew+a2,nnew-ynew+b2) # the new max likelihood estimate (MLE) at 6.7% is below that CRI // strongly influenced by the previous data
# if used the actual observed data (red curve, 10 times more data), we would have
qbeta(c(0.025,0.975),ynew+530+1,nnew-ynew+4848-530+1) # even more narrowly distributed around previous param. 
qbeta(c(0.025,0.975),ynew+a3,nnew-ynew+b3) # from the prior distribution that makes sense // closer to scenario 1, includes the MLE 


