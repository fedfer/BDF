### Libraries
library(bayesSurv)

### Generate Data
p = 20;k = 12; n = 700
L = matrix(rnorm(p*k,0,1),p,k)
X = scale(matrix(rnorm(n*k,0,1),n,k)%*%t(L))
#X = matrix(rnorm(n*p,0,1),n,p)

beta_true = numeric(p);beta_true[c(1,2,3)] = c(1,-1,1)
sigma_true = 0.5
y=as.vector(X%*%beta_true+rnorm(n,0,sigma_true))
df = as.data.frame(cbind(y,X))
summary(lm(y~. , data = df))
### 

X2 = t(X) %*% X

### Prior parameters
as = 1
bs = 0.3

### Initialition
Lambda = matrix(0, nrow = p, ncol = k)
Gamma = matrix(0, nrow = k, ncol = p)
psi = rgamma(p,as,bs)
csi = rgamma(k,as,bs)
Psi = diag(1/psi)
Csi = diag(1/csi)
eta = matrix(0,n,k)
theta = numeric(k)
sigmasq = 1

### Storage 
nrun = 1000
burn = 500
Lambda_st = array(0, dim=c(nrun-burn,p,k))
Gamma_st = array(0, dim=c(nrun-burn,k,p))
Psi_st = matrix(0, nrow = nrun-burn, ncol = p)
Csi_st = matrix(0, nrow = nrun-burn, ncol = k)
sigmasq_st = numeric(nrun-burn)
beta = matrix(0, nrow = nrun-burn, ncol = p)


for(s in 1:nrun){
   
   # --- update eta ---
   V = solve(diag(csi) + theta%*%t(theta)/sigmasq)
   for(i in 1:n){
      
      #m = V%*%(theta*y[i]/sigmasq + t(Lambda)%*%diag(psi)%*%X[i,])
      m = V%*%(theta*y[i]/sigmasq + diag(csi)%*%Gamma%*%X[i,])
      eta[i,] = bayesSurv::rMVNorm(1, 
                               mean = m, 
                               Sigma = V)
      
      
      # update without the regression
      #eta[i,] = bayesSurv::rMVNorm(1, 
      #                             mean = Gamma%*%X[i,] , 
      #                             Sigma = Csi)
   }
   
   # --- update Lambda 
   # No constraint
   #eta2 = t(eta) %*% eta    # prepare eta crossproduct before the loop
   #for(j in 1:p) {
   #   
   #   V = solve(eta2*psi[j] + diag(k))
   #   m = psi[j]*V%*%t(eta)%*%X[,j]
   #   Lambda[j,] = bayesSurv::rMVNorm(1, 
   #                                   mean = m, 
   #                                   Sigma = V)
   #}
   
   # Constraint
   Lambda = Psi%*%t(Gamma)%*%diag(csi)
   
   # --- update Gamma ---
   # No constraint
   for(j in 1:k) {
      
      V = solve(X2*csi[j] + diag(p))
      m = csi[j]*V%*%t(X)%*%eta[,j]
      Gamma[j,] = bayesSurv::rMVNorm(1, 
                                      mean = m, 
                                      Sigma = V)
   }

   # Constraint
   #Gamma = Csi%*%t(Lambda)%*%diag(psi)
   
   # --- update Psi ---
   X_til = X - eta %*% t(Lambda)
   psi = rgamma(p, as + 0.5*n, bs+0.5*colSums(X_til^2))
   Psi = diag(1/psi)
   
   
   # --- update Csi ---
   # No constraint
   eta_til = eta - X %*% t(Gamma)
   csi = rgamma(k, as + 0.5*n, bs+0.5*colSums(eta_til^2))
   Csi = diag(1/csi)
   
   # Constraint 
   #Csi = solve(Gamma%*%diag(psi)%*%Lambda + 0.001*diag(k))
   #csi = 1/diag(Csi)
   
   # --- update theta ---
   Vtheta = solve(eta2/sigmasq + diag(rep(1,k))/100)
   Mtheta = Vtheta%*%t(eta)%*%y/sigmasq
   theta = as.numeric(bayesSurv::rMVNorm(1,Mtheta,Vtheta))
   
   # --- update sigmasq ---
   y_til = y - eta %*% theta
   sigmasq = 1/rgamma(1, as + 0.5*n, bs+0.5*colSums(y_til^2))
   
   # Storage
   
   if(s > burn){
      Lambda_st[s-burn,,] = Lambda
      Gamma_st[s-burn,,] = Gamma
      Psi_st[s-burn,] = psi
      Csi_st[s-burn,] = csi
      #beta[s-burn,] = as.numeric(t(theta)%*%Gamma)
      beta[s-burn,] = as.numeric(t(theta)%*%Gamma)
      sigmasq_st[s-burn] = sigmasq
   }
   
   if(s%%200 == 0){
      print(s)
      print(mean((diag(psi)%*%Lambda - t(Gamma)%*%diag(csi))^2))
      print(mean((diag(psi)%*%Lambda - t(Gamma)%*%diag(csi))))
      print(diag((diag(psi)%*%Lambda - t(Gamma)%*%diag(csi))))
      
   }

   
}

apply(beta, 2, mean)
plot(apply(beta, 2, mean), beta_true)


