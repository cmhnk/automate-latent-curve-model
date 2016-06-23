# rm(list=ls())
require(MASS)
require(lavaan)
require(tidyr)

## generate data for 10 centers, each with N = 100 
## some centers have significant tx effect
## individuals randomly assigned to either control or placebo 

## Users must input:
center = 10
n = 200
t = 5
timecode = 0:(t-1)


N = center*n

## correlation matrices for latent intercept and slope
 mycor = c(-.4, -.3, -.2, -.1, .1, .2, .3, .4, .5, .6)  # hard coded

sig = list()
for (i in 1:center){
	sig[[i]] <- matrix(c(1, mycor[i],
					   mycor[i], 1), byrow=TRUE, ncol=2)
}

## generate values for the factors 
set.seed(34)
myslope.tx = rnorm(center, mean=-1, sd=1)
myslope.p = rnorm(center, mean=-0.5, sd=.5)
slope.dif = myslope.p - myslope.tx
mint = rnorm(center, mean=10, sd=1)

fac = list()
tx = c(rep(0, n/2), rep(1, n/2))
id = seq(1:N)


for (i in 1:center){
	
	a = mvrnorm(n = n/2, mu = c(mint[i], myslope.tx[i]), 
				 Sigma = sig[[i]], empirical = TRUE)
	b = mvrnorm(n = n/2, mu = c(mint[i], myslope.p[i]), 
				 Sigma = sig[[i]], empirical = TRUE)
	a.df = as.data.frame(a)
	b.df = as.data.frame(b)
	
	ab.df <- rbind(a.df, b.df)
	fac[[i]] <- ab.df
	fac[[i]] <- data.frame(tx, fac[[i]])
	names(fac[[i]]) <- c("tx", "alpha", "beta")
	fac[[i]] <- fac[[i]][with(fac[[i]], order(fac[[i]]$tx)),]
	
}

## add id's & center to each dataframe in fac
some.id <- list()
k = 1
m = 1
some.id[[1]] <- data.frame(id[k:n])
for (z in 2:center){
			k = k+1
			j = k*n
	some.id[[z]] <- data.frame(id[(m*n+1):j])
	m = m + 1 	
}

for (i in 1:center){
	fac[[i]] <- cbind(some.id[[i]], rep(i, n), fac[[i]])
	names(fac[[i]])[1:2] <- c("id", "center")
}

## sanity check
# ddply(fac[[1]], .(tx), summarize, Slope = mean(beta))


## generate observed variables 

ob <- c()

for (j in 1:center){
	for (i in 1:t){
 		ob[i] = paste("dep", i, sep="")
 		fac[[j]]$new <- c()
 		fac[[j]]$new <- fac[[j]]$alpha + timecode[i]*fac[[j]]$beta + rnorm(n=n, sd=.1) #was .7
 		replace <- which(names(fac[[j]])=="new")
 		names(fac[[j]])[replace]<-ob[i]	
	}  
}

# lapply(fac, summary) ## excellent, we have data! 

## now fit the model 
# hard coded example model
# model2 <-   'i =~ 1*dep1 + 1*dep2 + 1*dep3 + 1*dep4 + 1*dep5
# 			 s =~ 0*dep1 + 1*dep2 + 2*dep3 + 3*dep4 + 4*dep5'


## flexibly generated model syntax 

modsyn2 <- paste('i =~ ')
for (i in 1:length(ob)){
	if (i < length(ob)){
		modsyn2 <- strwrap(paste(modsyn2, '1*', ob[i], " + ", sep=""))	
	} else 
		modsyn2 <- strwrap(paste(modsyn2, '1*', ob[i], sep=""))	
}

modsyn2.2 <- paste('s =~ ')
for (i in 0:(length(ob)-1)){
	if (i < (length(ob)-1)){
		modsyn2.2 <- strwrap(paste(modsyn2.2, i, "*", ob[i+1], " +", " ", sep=""))	
	} else 
		modsyn2.2 <- strwrap(paste(modsyn2.2, i, '*', ob[i+1], sep=""))		
}

modsynfin <- paste(modsyn2, "\n", modsyn2.2, sep="")


## call model syntax to run model on each of the datasets in the list 

## run the model with everything fixed to equality between groups 
myfit1 <- list()
for (i in 1:center){
	myfit1[[i]] <- growth(modsynfin, data=fac[[i]], group="tx", group.equal=c("lv.variances", "means", "lv.covariances"))	
}

myparam1 <- list()
for (i in 1:center){
	myparam1[[i]] <- parameterEstimates(myfit1[[i]])[c(16:18, 24:25),]
	model <- "1"
	myparam1[[i]] <- data.frame(model, myparam1[[i]])
}

## run the model with the slopes freed between groups to test treatment effect 
myfit2 <- list()
for (i in 1:center){
	myfit2[[i]] <- growth(modsynfin, data=fac[[i]], group="tx", 
			  group.equal=c("lv.variances", "means", "lv.covariances"), 
			  group.partial=c("i~~s", "s~ 1", "s~~s"))	
}

myparam2 <- list()
for (i in 1:center){
	myparam2[[i]] <- parameterEstimates(myfit2[[i]])[c(16, 24, 17, 18, 25, 42, 43, 50),]
	model <- "2"
	myparam2[[i]] <- data.frame(model, myparam2[[i]])
	
}

## for each dataset, conduct LRT to determine if Mod 1 fits significantly worse than Mod 2
lrt <- list()
tx.difference <- list()
for (i in 1:center){
	lrt[[i]] <- anova(myfit2[[i]], myfit1[[i]])
	if (lrt[[i]][2,7] < .05){
		tx.difference[[i]] <- "TRUE"
	} else {
		tx.difference[[i]] <- "FALSE"
	}
}

## create new dataframe ready for Wordsmith

wordsmith <- list()

for (i in 1:center){
	center.id <- i
	significant.tx <- tx.difference[[i]]
	wordsmith[[i]] <- data.frame(center.id, significant.tx, rbind(myparam1[[i]], myparam2[[i]]))
	if (significant.tx=="TRUE"){
		wordsmith[[i]]<-wordsmith[[i]][which(wordsmith[[i]]$model=="2"),]
	} else {
		wordsmith[[i]]<-wordsmith[[i]][which(wordsmith[[i]]$model=="1"),]
	}
	wordsmith[[i]][,9:14] <- round(wordsmith[[i]][,9:14], 2)
	wordsmith[[i]]$new <- paste(wordsmith[[i]]$lhs, wordsmith[[i]]$op, wordsmith[[i]]$rhs, wordsmith[[i]]$group)
	# get rid of variables we don't want 
	wordsmith[[i]] <- wordsmith[[i]][,-c(3:8,11,13,14)]
}

## change list to one dataframe in long format 

wsdf.long <- do.call(rbind, wordsmith)
wsdf.long$center.id <- as.factor(wsdf.long$center.id)

wsdf.wide <- reshape(wsdf.long, timevar = "new", 
					 idvar = c("center.id", "significant.tx"), 
					 direction = "wide")
dim(wsdf.wide)
stder <- grep("se", names(wsdf.wide))

## get rid of standard error columns to clean things up
wsdf.wide <- wsdf.wide[,-stder]

## make names nice
thenames <- names(wsdf.wide)
thenames <- gsub("i ~~ i", "intercept.var", thenames)
thenames <- gsub("s ~~ s", "slope.var", thenames)
thenames <- gsub("i ~1", "intercept.mean", thenames)
thenames <- gsub("s ~1", "slope.mean", thenames)
thenames <- gsub("i ~~ s", "int.slope.cov", thenames)
thenames <- gsub(" 1", ".tx", thenames)
thenames <- gsub(" 2", ".control", thenames)
thenames <- gsub(" ", ".", thenames)
names(wsdf.wide) <- thenames

## save file as csv for Wordsmith! 

# write.table(wsdf.wide, "//Users//macuser//Dropbox//Automated Insights//10centers.csv", sep=",", row.names=FALSE, col.names=TRUE)

## show original data 
originaldata <- do.call(rbind, fac)
originaldata
wsdf.wide
