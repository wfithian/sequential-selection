set.seed(0)
X = read.table("X_100_200.csv", sep=',')
D = X
X = as.matrix(X) 
Y = as.matrix(read.table("y_100.csv", sep=','))
variables = as.matrix(read.table("variables.csv", sep=','))
D$Y = Y

BICmodel = step(lm(Y ~ 1, data=D), scope=list(upper= ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + 
                                                   V21 + V22 + V23 + V24 + V25 + V26 + V27 + V28 + V29 + V30 + V31 + V32 + V33 + V34 + V35 + V36 + V37 + V38 + V39 + 
                                                   V40 + V41 + V42 + V43 + V44 + V45 + V46 + V47 + V48 + V49 + V50 + V51 + V52 + V53 + V54 + V55 + V56 + V57 + V58 + V59 + 
                                                   V60 + V61 + V62 + V63 + V64 + V65 + V66 + V67 + V68 + V69 + V70 + V71 + V72 + V73 + V74 + V75 + V76 + V77 + V78 + V79 + 
                                                   V80 + V81 + V82 + V83 + V84 + V85 + V86 + V87 + V88 + V89 + V90 + V91 + V92 + V93 + V94 + V95 + V96 + V97 + V98 + V99 + 
						   V100 + V101 + V102 + V103 + V104 + V105 + V106 + V107 + V108 + V109 + V110 + V111 + V112 + V113 + V114 + V115 + V116 + 
						   V117 + V118 + V119 + V120 + V121 + V122 + V123 + V124 + V125 + V126 + V127 + V128 + V129 + V130 + V131 + V132 + V133 + 
						   V134 + V135 + V136 + V137 + V138 + V139 + V140 + V141 + V142 + V143 + V144 + V145 + V146 + V147 + V148 + V149 + V150 + 
						   V151 + V152 + V153 + V154 + V155 + V156 + V157 + V158 + V159 + V160 + V161 + V162 + V163 + V164 + V165 + V166 + V167 + 
						   V168 + V169 + V170 + V171 + V172 + V173 + V174 + V175 + V176 + V177 + V178 + V179 + V180 + V181 + V182 + V183 + V184 + 
						   V185 + V186 + V187 + V188 + V189 + V190 + V191 + V192 + V193 + V194 + V195 + V196 + V197 + V198 + V199 + V200), direction='forward', 
                                              k=log(100), trace=FALSE)

print("BIC")
print(BICmodel)

AICmodel = step(lm(Y ~ 1, data=D), scope=list(upper= ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + 
                                                   V21 + V22 + V23 + V24 + V25 + V26 + V27 + V28 + V29 + V30 + V31 + V32 + V33 + V34 + V35 + V36 + V37 + V38 + V39 + 
                                                   V40 + V41 + V42 + V43 + V44 + V45 + V46 + V47 + V48 + V49 + V50 + V51 + V52 + V53 + V54 + V55 + V56 + V57 + V58 + V59 + 
                                                   V60 + V61 + V62 + V63 + V64 + V65 + V66 + V67 + V68 + V69 + V70 + V71 + V72 + V73 + V74 + V75 + V76 + V77 + V78 + V79 + 
                                                   V80 + V81 + V82 + V83 + V84 + V85 + V86 + V87 + V88 + V89 + V90 + V91 + V92 + V93 + V94 + V95 + V96 + V97 + V98 + V99 + 
						   V100 + V101 + V102 + V103 + V104 + V105 + V106 + V107 + V108 + V109 + V110 + V111 + V112 + V113 + V114 + V115 + V116 + 
						   V117 + V118 + V119 + V120 + V121 + V122 + V123 + V124 + V125 + V126 + V127 + V128 + V129 + V130 + V131 + V132 + V133 + 
						   V134 + V135 + V136 + V137 + V138 + V139 + V140 + V141 + V142 + V143 + V144 + V145 + V146 + V147 + V148 + V149 + V150 + 
						   V151 + V152 + V153 + V154 + V155 + V156 + V157 + V158 + V159 + V160 + V161 + V162 + V163 + V164 + V165 + V166 + V167 + 
						   V168 + V169 + V170 + V171 + V172 + V173 + V174 + V175 + V176 + V177 + V178 + V179 + V180 + V181 + V182 + V183 + V184 + 
						   V185 + V186 + V187 + V188 + V189 + V190 + V191 + V192 + V193 + V194 + V195 + V196 + V197 + V198 + V199 + V200), direction='forward', 
                                              k=2, trace=FALSE)

print("AIC")
print(AICmodel)

print("step 8")
print(variables[1:8])
X8 = X[,variables[1:7]]
P8 = X8 %*% solve(t(X8) %*% X8) %*% t(X8)
A7 = as.matrix(read.table("A_step7.csv", sep=','))
A8 = as.matrix(read.table("A_step8.csv", sep=','))

print("step 12")
print(variables[1:12])
X12 = X[,variables[1:11]]
P12 = X12 %*% solve(t(X12) %*% X12) %*% t(X12)
A11 = as.matrix(read.table("A_step11.csv", sep=','))

sigma = as.numeric(read.table("sigma.csv"))

n = nrow(X)
p = ncol(X)
R8 = diag(rep(1,n)) - P8
R12 = diag(rep(1,n)) - P12

P8Y = P8 %*% Y
P8Y0 = P8 %*% rnorm(n) * sigma
P12Y = P12 %*% Y

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
        if(max(A %*% Z) < 0) {
            return(count)
        }
	if(count > 100000) {
            return(count) # give up...
        }
   }
}

maxT_sample = c()
for (i in 1:100) {
    maxT_sample = c(maxT_sample, wait_until_accept(A7, R8, P8Y, sigma))
}
print('maxT step8')
print(1. / mean(maxT_sample))

maxT_identify_sample = c()
for (i in 1:100) {
    maxT_identify_sample = c(maxT_identify_sample, wait_until_accept(A8, R8, P8Y, sigma))
}
print('maxT step 8, identifying variable')
print(1. / mean(maxT_identify_sample))

maxT12 = c()
for (i in 1:100) {
    maxT12 = c(maxT12, wait_until_accept(A11, R12, P12Y, sigma))
}
print('maxT step 12')
print(1. / mean(maxT12))

print('maxT step8 under global null with same model, to get one point...')
print(wait_until_accept(A7, R8, P8Y0, sigma))