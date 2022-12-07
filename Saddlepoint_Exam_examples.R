### -------------------------------------------- ###
### Saddlepoint approximation: general functions ###
### -------------------------------------------- ###



# Saddlepoint approximation: density
# sx: the saddlepoint times x
# K: the CGF evaluated at s - the saddlepoint
# ddK; the twice derrivative of K evaluated in the saddlepoint
# does not compute for values of x that lies on the edge of the support
simple_sp <- function(sx,K,ddK){
  return((1/sqrt(2*pi*ddK))*exp(K-sx))
}


# Saddlepoint normalization
# y_sp: a vector of values for a saddlepoint density,
# preferrably calculated up to as many points as possible.
# Must include some approximation to edge-values. 
sp_norm <- function(y_sp){
  return(y_sp/sum(y_sp))
}


#################### Binomial Example density ####################


# CGF
K_binom <- function(s,p,n){
  return(n*log(p*(exp(s)-1)+1))
}

# Saddlepoint
s_binom <- function(x,p,n){
  return(log(x*(1-p))-log((n-x)*p))
}

# Double derivative of CGF
ddK_binom <- function(x,n){
  return(x*(n-x)/n)
}


### -------------------------------------------------------- ###


p = 0.6                 # set p
n = 15                    # set n
x = seq(1,n-1,1)         # values to compute at (ignore 0 and n)


s = s_binom(x,p,n)
K = K_binom(s,p,n)
ddK = ddK_binom(x,n)


# compute the true density from R-studio
true_binom = dbinom(x,n,p)  

#compute the saddlepoint density
sp_reg_binom = simple_sp(s*x,K,ddK)

# add the true value at 0 for correct normalization
sp_normal_binom = c(dbinom(0,n,p),sp_reg_binom,dbinom(n,n,p))

#return a normalized saddlepoint without a value at 0
sp_normal_binom = sp_norm(sp_normal_binom)[-c(1,n)]

#normal approximation for comparison - continuity corrected
norm_approx2 = dnorm(x,n*p,sqrt(n*p*(1-p)))

#Poisson approximation for comparison - continuity corrected
pois_approx = dpois(x,n*p)

df_binom = data.frame(x = x, binomial = true_binom, sadpoint = sp_reg_binom, sadpoint_norm = sp_normal_binom, normal = norm_approx2, poisson = pois_approx)


# Overview plot
l_plot_x = 1
r_plot_x = 10

ggplot(df_binom, aes(x=x)) + 
  geom_step(aes(y=binomial, color = 'Binomial'), size = 1.1) +
  geom_step(aes(y=sadpoint, color = 'Saddlepoint')) +
  geom_step(aes(y=sadpoint_norm, color = 'Normalized saddlepoint')) +
  geom_step(aes(y=normal, color = 'Normal approximation')) +
  geom_step(aes(y=poisson, color = 'Poisson approximation')) +
  xlim(l_plot_x,r_plot_x) +
  labs(x = 'x',
       y = 'f(x)',
       color = 'Legend') + 
  theme_bw()

# Not gg_plot
#plot(x[l_plot_x:r_plot_x],true_p[l_plot_x:r_plot_x], type = 's') # density of Poisson
#lines(x[(l_plot_x):r_plot_x],sp_reg[l_plot_x:(r_plot_x)], col = 'red', type = 's') # sadpod approx
#lines(x[l_plot_x:r_plot_x],sp_normal[l_plot_x:r_plot_x],col = 'blue', type = 's') # sadpod + normalization
#lines(x[l_plot_x:r_plot_x], norm_approx[l_plot_x:r_plot_x], col = 'green', type = 's') # normal approximation


# Compute absolute error
abs_err_reg_sp = sum(abs(sp_reg_binom-true_binom))           # sadpod
abs_err_normalized_sp = sum(abs(sp_normal_binom-true_binom)) # normalized sadpod
abs_err_norm = sum(abs(norm_approx2-true_binom))             # normal approx
abs_err_pois = sum(abs(pois_approx-true_binom))              # poisson approx


# Look at the tails
tail_v = 5
true_binom[tail_v]
sp_reg_binom[tail_v]
sp_normal_binom[tail_v]
norm_approx2[tail_v]
pois_approx[tail_v]






#################### Poisson Example density ####################


library(ggplot2)


# CGF
K_pois <- function(s,lambda){
  return(lambda*(exp(s)-1))
}
# Saddlepoint
s_pois <- function(x,lambda){
  return(log(x/lambda))
}

# Double derivative of CGF
ddK_pois <- function(x){
  return(x)
}


l = 35                   # define lambda
x = seq(1,500,1)         # values to compute at (ignore 0)

s = s_pois(x,l)
K = K_pois(s,l)
ddK = ddK_pois(x)

# compute the true density from R-studio
true_p = dpois(x,l)    

#compute the saddlepoint density
sp_reg = simple_sp(s*x,K,ddK)

# add the true value at 0 for correct normalization
sp_normal = c(dpois(0,l),sp_reg)

#return a normalized saddlepoint without a value at 0
sp_normal = sp_norm(sp_normal)[-1]

#normal approximation for comparison - continuity corrected
norm_approx = dnorm(x,l,sqrt(l))


df_pois = data.frame(x = x, poisson = true_p, sadpoint = sp_reg, sadpoint_norm = sp_normal, normal = norm_approx)


# Overview plot
l_plot_x = 5
r_plot_x = 80

ggplot(df_pois, aes(x=x)) + 
  geom_step(aes(y=poisson, color = 'Poisson'), size = 1.1) +
  geom_step(aes(y=sadpoint, color = 'Saddlepoint')) +
  geom_step(aes(y=sadpoint_norm, color = 'Normalized saddlepoint')) +
  geom_step(aes(y=normal, color = 'Normal approximation')) +
  xlim(l_plot_x,r_plot_x) +
  labs(x = 'x',
       y = 'f(x)',
       color = 'Legend') + 
  theme_bw()

# Not gg_plot
#plot(x[l_plot_x:r_plot_x],true_p[l_plot_x:r_plot_x], type = 's') # density of Poisson
#lines(x[(l_plot_x):r_plot_x],sp_reg[l_plot_x:(r_plot_x)], col = 'red', type = 's') # sadpod approx
#lines(x[l_plot_x:r_plot_x],sp_normal[l_plot_x:r_plot_x],col = 'blue', type = 's') # sadpod + normalization
#lines(x[l_plot_x:r_plot_x], norm_approx[l_plot_x:r_plot_x], col = 'green', type = 's') # normal approximation


# Compute absolute error
abs_err_reg_sp = sum(abs(sp_reg-true_p)) # sadpod
abs_err_normalized_sp = sum(abs(sp_normal-true_p)) # normalized sadpod
abs_err_norm = sum(abs(norm_approx-true_p)) # normal approx



# Look at the tails
tail_v = 250
true_p[tail_v]
sp_reg[tail_v]
sp_normal[tail_v]
norm_approx[tail_v]



#################### Poisson Example p-value ####################




############## Poisson p-val functions #########################


w1 <- function(k,s,lambda){
  a = s*k-K_pois(s,lambda)
  return(sign(s)*sqrt(2*a))
}

# k must cannot be equal to lambda
p_approx1 <- function(k,lambda=1){
  y = numeric(length(k))
  k = k+1
  for (i in (1:length(y))){
    if (k[i] == lambda){
      y[i] = 1/2 - 1/(sqrt(2*pi))*(1/(6 * sqrt(lambda))-1/(2*sqrt(lambda)))
    }
    else {
    s = s_pois(k[i],lambda)
    w = w1(k[i],s,lambda)
    u = (1-(lambda/k[i]))*sqrt(k[i])
    y[i] = 1 - pnorm(w)-dnorm(w)*((1/w)-(1/u))
    }
  }
  return(y)
}

p_approx2 <- function(k,lambda=1){
  y = numeric(length(k))
  k = k+(1/2)
  for (i in (1:length(y))){
    s = s_pois(k[i],lambda)
    w = w1(k[i],s,lambda)
    u = 2*sinh(s/2)*sqrt(k[i])
    y[i] = 1 - pnorm(w)-dnorm(w)*((1/w)-(1/u))
  }
  return(y)
}

p_approx3 <- function(k,lambda=1){
  y = numeric(length(k))
  k = k+(1/2)
  for (i in (1:length(y))){
    s = s_pois(k[i],lambda)
    w = w1(k[i],s,lambda)
    u = s*sqrt(k[i])
    y[i] = 1 - pnorm(w)-dnorm(w)*((1/w)-(1/u))
  }
  return(y)
}



x = seq(10,450,1) 
x = c(10,50,125)

p1 = p_approx1(x,l)
p2 = p_approx2(x,l)
p3 = p_approx3(x,l)
pt = ppois(x,l,lower.tail = F)
pn = pnorm(x+1/2,l,sqrt(l),lower.tail = F)

p_num_int = numeric(500)

for (i in 1:500){
  p_num_int[i] = sum(sp_normal[(500-i):500])
}

p_num_int = rev(p_num_int)[-c(1,2)]


df_pois = data.frame(x = x, Poisson = pt, Pr1 = p1, Pr2 = p2, Pr3 = p3, summed = p_num_int[x], Normal = pn)


# Overview plot

ggplot(df_pois, aes(x=x)) + 
  geom_step(aes(y=Poisson, color = 'Poisson'), size = 1.1) +
  geom_step(aes(y=Pr1, color = 'Pr1')) +
  geom_step(aes(y=Pr2, color = 'Pr2')) +
  geom_step(aes(y=Pr3, color = 'Pr3')) +
  geom_step(aes(y=summed, color = 'Sum of density')) +
  geom_step(aes(y=Normal, color = 'Normal approximation')) +
  labs(x = 'x',
       y = 'P(X>x)',
       color = 'Legend') + 
  theme_bw()

#plot(x,pt, type = 's')
#lines(x,p1,col = 'red',  type = 's')
#lines(x,p2,col = 'blue', type = 's')
#lines(x,p3,col = 'orange', type = 's')
#lines(x,pn,col = 'green', type = 's')
#lines(x,p_num_int[x+2], col = 'purple', type = 's')

sum(abs(pt-p1))
sum(abs(pt-p2))
sum(abs(pt-p3))
sum(abs(pt-pn))
sum(abs(pt-p_num_int[x]))

b = 0.01
which(pt < b)[1]
which(p1 < b)[1]
which(p2 < b)[1]
which(p3 < b)[1]
which(pn < b)[1]
which(p_num_int < b)[1]
#data.frame(pt,p1,p2,p3,pn)

