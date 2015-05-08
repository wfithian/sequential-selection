set.seed(0)
X = read.table("X_100_200.csv", sep=',')
X = as.matrix(X) 
Y = as.matrix(read.table("y_100.csv", sep=','))
variables = as.matrix(read.table("variables.csv", sep=','))

print("step 7")
print(variables[1:7])
X7 = X[,variables[1:6]]

print("step 8")
print(variables[1:8])
X8 = X[,variables[1:7]]

P7 = X7 %*% solve(t(X7) %*% X7) %*% t(X7)
P8 = X8 %*% solve(t(X8) %*% X8) %*% t(X8)

A7 = as.matrix(read.table("A_step7.csv", sep=','))
A8 = as.matrix(read.table("A_step8.csv", sep=','))

sigma = as.numeric(read.table("sigma.csv"))

n = nrow(X)
p = ncol(X)
R7 = diag(rep(1,n)) - P7
R8 = diag(rep(1,n)) - P8

P7Y = P7 %*% Y
P8Y = P8 %*% Y

# here is the null distribution

null_sample = function(R, PY, sigma) {
     n = nrow(R)
     Z = PY + R %*% rnorm(n) * sigma
     return(Z)
}

# how many tries until we accept

wait_until_accept = function(A, R, PY, sigma) {
    count = 0
    while(TRUE) {
        count = count + 1
        Z = null_sample(R, PY, sigma)
	print(max(A %*% Z))
        if(max(A %*% Z) < 0) {
            return(count)
        }
   }
}

count7_sample = c()
for (i in 1:100) {
    count7_sample = c(count7_sample, wait_until_accept(A7, R7, P7Y, sigma))
    print(1. / mean(count7_sample))
}

count8_sample = c()
for (i in 1:100) {
    count8_sample = c(count8_sample, wait_until_accept(A8, R8, P8Y, sigma))
    print(1. / mean(count8_sample))
}
