########################################################################
####### This is file "euro11.R"
####### This version: 2014-02-04
####### I run the Synthetic control analyses for the Euro 11 countries
########################################################################




#### Let's begin with defining some predictors:

pred <- names(scdata)[c(#
	7, # pop65+
	8, # pop15-
	#9,# unemployment (1)
	#10,# system
	11,# yrcurnt
	12,# allhouse
	13,# legelec
	14, # pluralty
	15, # pr
	#16, # checks
	#17, # fri
	#18,# GDP growth
	#19,# growth in multi factor productivity
	20,# labor productivity annual growth
	21, # health expenditure/GDP
	22,# GDP expenditure approach
	#23,# tax revenue %GDP (general)
	24,# taxrev %GDP, no Social Security
	25, # CO2 emissions
	#26, # FDI
	#27, # GDP growth
	#28,# Gini index
	#29,#,# Inflation (Consumer Prices)
	#30, # Poverty
	31,#, # unemployment (World Bank)
	##32, #Population
	#33,#, #openness (PWT)
	34,#,# openness (expenditure)
	#35, # Expenditure on Families %GDP
	36, # PolconIII
	#37, # PolconV
	38,  # Potrafke ideology
	39, # Majority margin
	#40, # Herfindahl Index Government
	41, #lag debt/gdp (RR)
	42#,# Rae Fractionalisation index (government)
	#43 # Rae Fractionalisation Index (total)
	)]

#### define countries used for synthetic control group

contr <- sort(unique(scdata$ccode[is.element(scdata$country,setdiff(scdata$country,c(Euro12,"Euro 11","Slovenia")))]))

#### The following countries have to be excluded due to data constraints (missing values)
contr <- setdiff(contr, c(1111,2222,70,155,225,269,290,310,316,317,349,355,360,366,666,732,sort(unique(scdata$ccode[scdata$eu==0 & scdata$oecd==0]))))


#### Show countries that are available for Synthetic Greece
country <- sort(unique(scdata$country[scdata$ccode %in% contr]))

################ Quang's work for in-space placebo tests ##############
gaps.storage <- matrix(1:308, nrow=11, ncol=28)
for (ind in 1:length(c(0, contr))){
  id = c(0, contr)[ind]
  #### Generate Synth object (debt/gdp ratio) to be used to run the analysis
  sdata <- dataprep(foo = scdata[scdata$ccode %in% contr | scdata$country == "Euro 11",],
                    predictors = pred,
                    dependent = names(scdata[6]),
                    unit.variable = "ccode",
                    time.variable = "Year",
                    treatment.identifier = id,
                    controls.identifier = setdiff(contr, c(id)),
                    time.predictors.prior = c(1983:1998),
                    time.optimize.ssr = c(1983:1999),
                    unit.names.variable = "country",
                    time.plot = 1983:2010
  )

  #### Run the synthetic control analysis:
  synth.out <- synth(data.prep.obj = sdata, method = "BFGS")

  #### calculate output gaps from the results
  gaps <- sdata$Y1plot - (sdata$Y0plot %*% synth.out$solution.w)
  gaps.storage[ind,] <- gaps
}

# export results
write.csv(gaps.storage, file='placebo_test.csv')

# to read matrix from 'placebo_test.csv'
gaps.storage <- read.csv("placebo_test.csv")
gaps.storage <- as.matrix(gaps.storage)
gaps.storage <- gaps.storage[, 2:29]


# get p-values for each year
p.storage <- c()
for (i in 1:28){
	percentile = ecdf(-abs(gaps.storage[,i]))
  p.storage[i] <- percentile(-abs(gaps.storage[1,i]))
}

# plotting
Ylim <- c(
  -(.3*max(abs(gaps.storage[1,])) + max(abs(gaps.storage[1,]))),
  (.3*max(abs(gaps.storage[1,])) + max(abs(gaps.storage[1,])))
)
plot(1983:2010, gaps.storage[1,], t='l',
     col='black', lwd=2, main=c("Debt/GDP gaps in euro 11 and placebo gaps in all 10 control countries"),
     ylab=c('Gap in Debt/GDP (percentage points, 1983-2010)'),
     xlab=c('Year'), ylim=Ylim)
for (i in 2:10){
  lines(1983:2010, gaps.storage[i,], lwd=1, col='grey')
}
abline(h=0, col='black', lty='dashed',lwd=2)
abline(v=1999, col='black', lty='dotted', lwd=2)
lines(1983:2010, gaps.storage[1,], lwd=2, col='black')
legend(1985, 30, legend=c('Euro 11',
                           'Control countries'),
       col=c("black", "grey"), lty=c(1, 1), lwd=c(2,1),cex=0.8)

# a variation of the plot where MSPE > 2 times MSPE of euro 11 is eliminated
# calculate MSPE
mspe.storage <- c()
for (i in 1:11){
  pre.treat <- gaps.storage[i,1:16]
  mspe.storage[i] <- sum(pre.treat**2)/16
}

not.cutoff <- mspe.storage <= 20*mspe.storage[1]
Ylim <- c(
  -(.3*max(abs(gaps.storage[1,])) + max(abs(gaps.storage[1,]))),
  (.3*max(abs(gaps.storage[1,])) + max(abs(gaps.storage[1,])))
)
plot(1983:2010, gaps.storage[1,], t='l',
     col='black', lwd=2, main=c("Debt/GDP gaps in euro 11 and
                                placebo gaps in all 10 control countries"),
     ylab=c('Gap in Debt/GDP (percentage points, 1983-2010)'),
     xlab=c('Year'), ylim=Ylim)
for (country in (1:11)[not.cutoff][2:sum(not.cutoff)]){
  lines(1983:2010, gaps.storage[country,], lwd=1, col='grey')
}
abline(h=0, col='black', lty='dashed',lwd=2)
abline(v=1999, col='black', lty='dotted', lwd=2)
lines(1983:2010, gaps.storage[1,], lwd=2, col='black')
legend(1985, 30, legend=c('Euro 11',
                           'Control countries'),
       col=c("black", "grey"), lty=c(1, 1), lwd=c(2,1), cex=0.8)

################# End of Quang's work #####################




################ Khoi's work for p-value estimates using maximum likelihood ##############
mle.model <- function(ctrl.outcomes) { # fits a ditribution on data
   mean.est <- sum(ctrl.outcomes) / length(ctrl.outcomes)
   std.est <- sqrt(1 / length(ctrl.outcomes) * sum((ctrl.outcomes - mean.est)^2))
   return(c(mean.est, std.est))
}

calc.prob <- function(mean, std, xs) {
  ys <- 1 / (std * sqrt(2*pi)) * exp(-(xs - mean)^2 / (2 * std^2))
  return(ys)
}

#### MLE for all outcomes ####

act.trt <- abs(gaps.storage[1, 28]) # actual treated unit 2010 effect
synth.trt.eff <- abs(gaps.storage[2:11, 28]) # 2010 "treatment" effect magnitudes for control units
fitted.model <- mle.model(synth.trt.eff)
synth.trt.prob <- calc.prob(fitted.model[1], fitted.model[2], synth.trt.eff)
plot(density(x = synth.trt.eff, weights = synth.trt.prob))
fitted.model

synth.trt.eff
synth.trt.prob
act.trt

# create a sequence of numbers between 0 and 150 for magnitude scale
x <- seq(0, 150, by = .1)

# values of normal distribution
y <- dnorm(x, mean = fitted.model[1], sd = fitted.model[2])

# set significance threshod
alpha <- 0.05

# control unit results
plot(synth.trt.eff, synth.trt.prob, pch = 19, col = "blue",
     ylab = "Density", xlab = "Magnitude of Debt/GDP gap (percentage points)")
lines(x,y, type = "l")

# actual treatment unit results
abline(v = act.trt, col = 'red', lwd = 1)
abline(v = qnorm(1 - alpha, fitted.model[1], fitted.model[2]), # 95th percentile
       col = 'black', lty='dotted', lwd = 1) # one-sided test
legend(50, 0.002, legend=c(paste("ATE = ", as.character(round(act.trt, digits = 2))),
                           paste("alpha = ", as.character(round(alpha, digits = 2)))),
       col=c("red", "black"), lty=c(1, 3), cex=0.8)



# p-value for treatment unit = 0.673
1 - pnorm(act.trt, fitted.model[1], fitted.model[2])

#### End of MLE 1 ####


#### MLE for trimmed outcomes ####
nocut.idx <- which(not.cutoff)[-1] # remove treatment unit
nocut.idx
synth.trt.eff.nocut <- abs(gaps.storage[nocut.idx, 28]) # 2010 "treatment" effect magnitudes for non-outlier control units 
fitted.model2 <- mle.model(synth.trt.eff.nocut)
synth.trt.prob.nocut <- calc.prob(fitted.model2[1], fitted.model2[2], synth.trt.eff.nocut)
plot(density(x = synth.trt.eff.nocut, weights = synth.trt.prob.nocut))
fitted.model2

<<<<<<< HEAD
synth.trt.eff.nocut
synth.trt.prob.nocut
=======
synth.trt.eff <- abs(gaps.storage[2:11, 28]) # 2010 "treatment" effect magnitudes for control units
fitted.model <- mle.model(synth.trt.eff)
synth.trt.prob <- calc.prob(fitted.model[1], fitted.model[2], synth.trt.eff)
plot(density(x = synth.trt.eff, weights = synth.trt.prob))
fitted.model

synth.trt.eff
synth.trt.prob
>>>>>>> 40180a31154808d9941389e0872674f7a62e580d

# create a sequence of numbers between 0 and 150 for magnitude scale
x <- seq(0, 70, by = .1)

# values of normal distribution
y <- dnorm(x, mean = fitted.model2[1], sd = fitted.model2[2])

# control unit results
plot(synth.trt.eff.nocut, synth.trt.prob.nocut, pch = 19, col = "blue",
     ylab = "Density", xlab = "Magnitude of Debt/GDP gap (percentage points)"
     , xlim = c(0, 70), ylim = c(0, 0.03))
lines(x,y, type = "l")

abline(v = act.trt, col = 'red', lwd = 1)
abline(v = qnorm(1 - alpha, fitted.model2[1], fitted.model2[2]), # 95th percentile
       col = 'black', lty='dotted', lwd = 1) # one-sided test
legend(35, 0.005, legend=c(paste("ATE = ", as.character(round(act.trt, digits = 2))),
                           paste("alpha = ", as.character(round(alpha, digits = 2)))),
       col=c("red", "black"), lty=c(1, 3), cex=0.8)

# p-value for treatment unit = 0.657
1 - pnorm(act.trt, fitted.model2[1], fitted.model2[2])

#### End of MLE 2 ####

################# End of my work #####################




#### Generate Synth object (debt/gdp ratio) to be used to run the analysis
sdata <- dataprep(foo = scdata[scdata$ccode %in% contr | scdata$country == "Euro 11",],
		predictors = pred,
		dependent = names(scdata[6]),
		unit.variable = "ccode",
		time.variable = "Year",
		treatment.identifier = 0,
		controls.identifier = contr,
		time.predictors.prior = c(1983:1998),
		time.optimize.ssr = c(1983:1999),
		unit.names.variable = "country",
		time.plot = 1983:2010
)


#### Run the synthetic control analysis:
synth.out <- synth(data.prep.obj = sdata, method = "BFGS")

#### calculate output gaps from the results
gaps <- sdata$Y1plot - (sdata$Y0plot %*% synth.out$solution.w)


####
synth.tables <- synth.tab(dataprep.res = sdata, synth.res=synth.out)





#### Plot the Path of the Debt to GDP ratio for the Euro11 and the Synthetic control
path.plot(synth.res = synth.out,
	dataprep.res = sdata,
	Ylab="Debt/GDP (Nominal)",
	Xlab="Year",
	Legend=c("Euro 11","Synthetic Euro 11"),
	Legend.position="bottomright", abline(v=1999,lty="dashed")
)



#### Plot the gap of the Debt to GDP ratio for the Euro11 and the Synthetic control
gaps.plot(synth.res = synth.out, dataprep.res = sdata,
	Ylab= "Gap in Debt/GDP (percentage points, 1983-2010)", Xlab="Year",
	Main=NA, abline(v=1999,lty="dashed"))



#### extract Country weights from Synthetic control
a <- cbind(synth.tables$tab.w[,1],synth.tables$tab.w[,2])
row.names(a) <- synth.tables$tab.w[,2]

#### Plot Country weights
dotchart(a[,1],pch=16)

#### store estimation results in list
scresults[[1]] <- synth.out

#### Print estimation results on screen
print(synth.out)

# generate table 3
xtable(synth.tables$tab.pred)
