pdf('maxT_view.pdf')
D1 = read.table('test1.csv', sep=',', header=TRUE)
D2 = read.table('test2.csv', sep=',', header=TRUE)
D3 = read.table('test3.csv', sep=',', header=TRUE)
D4 = read.table('test4.csv', sep=',', header=TRUE)
D = rbind(D1,D2,D3,D4)
write.table(D, "test.csv", row.names=FALSE, col.names=TRUE, sep=',')

par(mfrow=c(2,2))

# mostly signal

plot(ecdf(D$maxT_pvalue_4), main='Step 4', xlim=c(0,1))
plot(ecdf(D$maxT_identify_pvalue_4), add=TRUE, col='blue')
plot(ecdf(D$maxT_identify_unknown_pvalue_4), add=TRUE, col='yellow')
plot(ecdf(D$maxT_unknown_pvalue_4), add=TRUE, col='green')
plot(ecdf(D$saturated_pvalue_4), add=TRUE, col='red')

plot(ecdf(D$maxT_pvalue_5), main='Step 5', xlim=c(0,1))
plot(ecdf(D$maxT_identify_pvalue_5), add=TRUE, col='blue')
plot(ecdf(D$saturated_pvalue_5), add=TRUE, col='red')
plot(ecdf(D$maxT_identify_unknown_pvalue_5), add=TRUE, col='yellow')
plot(ecdf(D$maxT_unknown_pvalue_5), add=TRUE, col='green')

plot(ecdf(D$maxT_pvalue_6), main='Step 6', xlim=c(0,1))
plot(ecdf(D$maxT_identify_pvalue_6), add=TRUE, col='blue')
plot(ecdf(D$saturated_pvalue_6), add=TRUE, col='red')
plot(ecdf(D$maxT_identify_unknown_pvalue_6), add=TRUE, col='yellow')
plot(ecdf(D$maxT_unknown_pvalue_6), add=TRUE, col='green')

plot(ecdf(D$maxT_pvalue_7), main='Step 7', xlim=c(0,1))
plot(ecdf(D$maxT_identify_pvalue_7), add=TRUE, col='blue')
plot(ecdf(D$saturated_pvalue_7), add=TRUE, col='red')
plot(ecdf(D$maxT_identify_unknown_pvalue_7), add=TRUE, col='yellow')
plot(ecdf(D$maxT_unknown_pvalue_7), add=TRUE, col='green')

# mixed

plot(ecdf(D$maxT_pvalue_8), main='Step 8', xlim=c(0,1))
plot(ecdf(D$maxT_identify_pvalue_8), add=TRUE, col='blue')
plot(ecdf(D$saturated_pvalue_8), add=TRUE, col='red')
plot(ecdf(D$maxT_identify_unknown_pvalue_8), add=TRUE, col='yellow')
plot(ecdf(D$maxT_unknown_pvalue_8), add=TRUE, col='green')

plot(ecdf(D$maxT_pvalue_9), main='Step 9', xlim=c(0,1))
plot(ecdf(D$maxT_identify_pvalue_9), add=TRUE, col='blue')
plot(ecdf(D$saturated_pvalue_9), add=TRUE, col='red')
plot(ecdf(D$maxT_identify_unknown_pvalue_9), add=TRUE, col='yellow')
plot(ecdf(D$maxT_unknown_pvalue_9), add=TRUE, col='green')

plot(ecdf(D$maxT_pvalue_10), main='Step 10', xlim=c(0,1))
plot(ecdf(D$maxT_identify_pvalue_10), add=TRUE, col='blue')
plot(ecdf(D$saturated_pvalue_10), add=TRUE, col='red')
plot(ecdf(D$maxT_identify_unknown_pvalue_10), add=TRUE, col='yellow')
plot(ecdf(D$maxT_unknown_pvalue_10), add=TRUE, col='green')

plot(ecdf(D$maxT_pvalue_11), main='Step 11', xlim=c(0,1))
plot(ecdf(D$maxT_identify_pvalue_11), add=TRUE, col='blue')
plot(ecdf(D$saturated_pvalue_11), add=TRUE, col='red')
plot(ecdf(D$maxT_identify_unknown_pvalue_11), add=TRUE, col='yellow')
plot(ecdf(D$maxT_unknown_pvalue_11), add=TRUE, col='green')

# mostly nulls

{ 
plot(ecdf(D$maxT_pvalue_20), main='Step 20')
plot(ecdf(D$maxT_identify_pvalue_20), add=TRUE, col='blue')
plot(ecdf(D$maxT_identify_unknown_pvalue_20), add=TRUE, col='yellow')
plot(ecdf(D$maxT_unknown_pvalue_20), add=TRUE, col='green')

plot(ecdf(D$maxT_pvalue_21), main='Step 21')
plot(ecdf(D$maxT_identify_pvalue_21), add=TRUE, col='blue')
plot(ecdf(D$maxT_identify_unknown_pvalue_21), add=TRUE, col='yellow')
plot(ecdf(D$maxT_unknown_pvalue_21), add=TRUE, col='green')

plot(ecdf(D$maxT_pvalue_22), main='Step 22')
plot(ecdf(D$maxT_identify_pvalue_22), add=TRUE, col='blue')
plot(ecdf(D$maxT_identify_unknown_pvalue_22), add=TRUE, col='yellow')
plot(ecdf(D$maxT_unknown_pvalue_22), add=TRUE, col='green')

plot(ecdf(D$maxT_pvalue_23), main='Step 23')
plot(ecdf(D$maxT_identify_pvalue_23), add=TRUE, col='blue')
plot(ecdf(D$maxT_identify_unknown_pvalue_23), add=TRUE, col='yellow')
plot(ecdf(D$maxT_unknown_pvalue_23), add=TRUE, col='green')

# later in the path

plot(ecdf(D$maxT_pvalue_30), main='Step 30')
plot(ecdf(D$maxT_identify_pvalue_30), add=TRUE, col='blue')
plot(ecdf(D$maxT_pvalue_31), main='Step 31')
plot(ecdf(D$maxT_identify_pvalue_31), add=TRUE, col='blue')
plot(ecdf(D$maxT_pvalue_33), main='Step 33')
plot(ecdf(D$maxT_identify_pvalue_33), add=TRUE, col='blue')
plot(ecdf(D$maxT_pvalue_33), main='Step 33')
plot(ecdf(D$maxT_identify_pvalue_33), add=TRUE, col='blue')

# later in the path

plot(ecdf(D$maxT_pvalue_36), main='Step 36')
plot(ecdf(D$maxT_identify_pvalue_36), add=TRUE, col='blue')
plot(ecdf(D$maxT_pvalue_37), main='Step 37')
plot(ecdf(D$maxT_identify_pvalue_37), add=TRUE, col='blue')
plot(ecdf(D$maxT_pvalue_38), main='Step 38')
plot(ecdf(D$maxT_identify_pvalue_38), add=TRUE, col='blue')
plot(ecdf(D$maxT_pvalue_39), main='Step 39')
plot(ecdf(D$maxT_identify_pvalue_39), add=TRUE, col='blue')

}

dev.off()