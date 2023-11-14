Gamma = matrix(c(0.9,0.05,0.05,
                 0.7,0,0.3,
                 0.8,0,0.2),
                 byrow=TRUE,nrow=3,ncol=3)

library(expm)
p0 = c(1,0,0)
p2 = p0 %*% (Gamma%^%2)
p2
p10 = p2 %*% (Gamma%^%8)
p10

p0 = c(0,0,1)
p10 = p0 %*% (Gamma%^%10)
p10



