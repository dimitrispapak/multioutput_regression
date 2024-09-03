# Packages :
# R version 4.1.3
# interp_1.1-4
# knitr_1.43
# readxl_1.4.3
# dplyr_1.1.0
# randomForestSRC_3.2.2
# stringr_1.5.0
library(dplyr)
library(stringr)
library(knitr)
library(readxl)
library(gplots)
library(randomForestSRC)

perf_metrics <- function(actual,pred){
	na_indeces <- which(is.na(actual))
	if (length(na_indeces) > 0){
		actual <- actual[-na_indeces]
		pred <- pred[-na_indeces]
	}
	print(paste("mean actual",mean(actual)))
	mae <- sum(abs(actual - pred))/length(actual)
	rmse <- sqrt(mean((actual - pred)^2)) 
	rho <- cor.test(actual, pred, method = 'spearman')$estimate
	r_squared <- (cor(actual,pred))^2
	return(c('mae' = mae,'rmse' = rmse,'rho'= rho,'r^2'= r_squared))
}


qq_plot <- function(actual,pred,target,log=F){
	na_indeces <- which(is.na(actual))
	if (length(na_indeces) > 0){
		actual <- actual[-na_indeces]
		pred <- pred[-na_indeces]
	}
	if (log == T){
		actual = log(actual)
		pred = log(pred)
	}

	min_ = min(min(actual),min(pred))
	max_ = max(max(actual),max(pred))


	plot(x = pred, y = actual,
	     xlab = 'Predicted ',
	     ylab = 'Actual',
	     ylim = c(min_,max_),
	     xlim = c(min_,max_),
	     main = paste(target,'Predicted vs. Actual'))
	# Adding a diagonal line for the ideal prediction
	abline(a = 0, b = 1)
	# Fit the linear model
	df = data.frame(y = actual, x = pred)
	model <- lm(y ~ x ,df)
	# Add the least squares line to the plot
	abline(model, col = "red")
	# Get the coefficients of the line
	coefficients <- coef(model)
	# Calculate R-squared
	r_squared <- summary(model)$r.squared
	# spearman correlation rho value
	rho <- cor.test(actual, pred, method = 'spearman')$estimate
	# Create the legend text
	#legend_text <- paste("y = ", round(coefficients[2], 2), "x + ", round(coefficients[1], 2), "\nR^2 = ", round(r_squared, 2))
	legend_text <- paste("y = ", round(coefficients[2], 2), "x + ", round(coefficients[1], 2), "\nR^2 = ", round(r_squared, 2),",\tSpearman correlation rho=",round(rho,2))
	# Add the legend to the plot
	legend("topleft", legend = legend_text, col = "red", lty = 1,bg='grey')
}

train_test_ratio = 0.9
current_date <- format(Sys.Date(), "%d_%m_%Y")
train_percent <- train_test_ratio * 100
test_percent <- (1 - train_test_ratio)* 100
pdf_file = paste0(train_percent,'-',test_percent,'_results_',current_date,'.pdf')
pdf(pdf_file,width=12,height=9)
#data <- read.csv('data_29_08_23.csv')
#data <- read.csv('data_29_03_24_cleaned.csv')
data <- as.data.frame(read_excel('/mnt/beegfs/userdata/d_papakonstantinou/damage/TelikogiaDim.xlsx'))
########################## Data Preparation  ###########################
print(colnames(data))
data$alpha <- as.numeric(data$alpha)
data$beta <- as.numeric(data$beta)
data$LET <- as.numeric(data$LET)
data$RBE <- as.numeric(data$RBE)
data$DoseRate <- as.numeric(data$DoseRate)
data$Energy <- as.numeric(data$Energy)
data$DSBs <- as.numeric(data$DSBs)
data$nonDSBClusters <- as.numeric(data$nonDSBClusters)

hist(data$alpha)
hist(data$beta)
hist(data$LET)
hist(data$RBE)
hist(data$DoseRate)
hist(data$Energy)

# Removal of mostly empty features
data$RBE <-NULL
data$Energy <-NULL
data$DoseRate <-NULL
data$X <-NULL
data$`#ExpID` <-NULL
data$PMID <-NULL
data$`#Exp` <-NULL

kable(head(data))
data <- data %>% mutate(CellLine = str_trim(CellLine))
data <- data %>% mutate(CellLine = toupper(CellLine))
data <- data %>% mutate(CellClass = str_trim(CellClass))
data <- data %>% mutate(IrradiationConditions = str_trim(IrradiationConditions))
data <- data %>% mutate(Tissue = str_trim(Tissue))

data <- data %>% mutate(Tissue = toupper(Tissue))
print(dim(data))
data <- data[!is.na(data$LET),]
print(dim(data))

data$RadiationType[data$RadiationType == "α-particles"] <- "alpha-particles" # remove greek letters
data$RadiationType[data$RadiationType == "γ-rays"] <- "gamma-rays" # remove greek letters
data$RadiationType         <- as.factor(data$RadiationType)
data$CellLine              <- as.factor(data$CellLine)
data$Tissue                <- as.factor(data$Tissue)
data$CellClass             <- as.factor(data$CellClass)
data$CellCycle             <- as.factor(data$CellCycle)
data$IrradiationConditions <- as.factor(data$IrradiationConditions)
sapply(data,class)

barplot(table(data$CellCycle)/length(data$CellCycle))
barplot(table(data$IrradiationConditions)/length(data$IrradiationConditions))
barplot(table(data$RadiationType)/length(data$RadiationType))
barplot(table(data$CellClass)/length(data$CellClass))
print(dim(data))
#######################################################################
smp_size <- floor(train_test_ratio * nrow(data))
## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(data)), size = smp_size)

train <- data[train_ind, ]
test <- data[-train_ind, ]

print(colnames(test))


o <- tune(Multivar(alpha, beta) ~., train)
print('optimal parameters:')
print(o$optimal)
## visualize the nodesize/mtry OOB surface
if (library("interp", logical.return = TRUE)) {

	## nice little wrapper for plotting results
	plot.tune <- function(o, linear = TRUE) {
		x <- o$results[,1]
		y <- o$results[,2]
		z <- o$results[,3]
		so <- interp(x=x, y=y, z=z, linear = linear)
		idx <- which.min(z)
		x0 <- x[idx]
		y0 <- y[idx]
		filled.contour(x = so$x,
			       y = so$y,
			       z = so$z,
			       xlim = range(so$x, finite = TRUE) + c(-2, 2),
			       ylim = range(so$y, finite = TRUE) + c(-2, 2),
			       color.palette =
			       colorRampPalette(c("yellow", "red")),
			       xlab = "nodesize",
			       ylab = "mtry",
			       main = "error rate for nodesize and mtry",
			       key.title = title(main = "OOB error", cex.main = 1),
			       plot.axes = {axis(1);axis(2);points(x0,y0,pch="x",cex=1,font=2);
			       points(x,y,pch=16,cex=.25)})
	}
	## plot the surface
	plot.tune(o)
}
# Fit model with optimal parameters
model <- rfsrc(Multivar(alpha, beta) ~., train, importance = TRUE, mtry = o$optimal['mtry'], nodesize=o$optimal['nodesize'])

par(mfrow=c(2,2))
print(model)
plot(model)
## performance error and vimp
par(mfrow=c(1,1))
vmp <- get.mv.vimp(model)
print(vmp)
print('alpha')
print(vmp[,'alpha'])
# Create a bar plot without x-axis labels
par(mar = c(5, 10, 4, 1))
barplot(sort(vmp[,'alpha'],decreasing=F), main = "Variable Importance (alpha)", col = "lightblue",horiz=T, las = 1, border = "black", names.arg = NULL)
par(mar = c(5, 10, 4, 1))
barplot(sort(vmp[,'beta'],decreasing=F), main = "Variable Importance (beta)", col = "lightblue",horiz=T, las = 1, border = "black", names.arg = NULL)

pred <- predict.rfsrc(model, test,na.action='na.impute')
pred_alpha = pred$regrOutput$alpha$predicted
pred_beta = pred$regrOutput$beta$predicted

## standardized error and vimp
err.std <- get.mv.error(model, standardize = TRUE)
vmp.std <- get.mv.vimp(model, standardize = TRUE)

print(err.std)
print(vmp.std)
## qqplots 
par(mfrow=c(1,1))
qq_plot(test[,'alpha'],pred_alpha,'alpha')
res1 <- perf_metrics(test[,'alpha'],pred_alpha)

par(mfrow=c(1,1))
qq_plot(test[,'alpha'],pred_alpha,'log(alpha)',log=T)

indices <- which(test[,'alpha'] >= 0 & test[,'alpha'] <= 2)
par(mfrow=c(1,1))
qq_plot(test[,'alpha'][indices],pred_alpha[indices],'alpha (0 <-> 2)')
res2 <- perf_metrics(test[,'alpha'][indices],pred_alpha[indices])

par(mfrow=c(1,1))
qq_plot(test[,'beta'],pred_beta,'beta')
res3 <- perf_metrics(test[,'beta'],pred_beta)

par(mfrow=c(1,1))
indices <- which(test[,'beta'] >= 0 & test[,'beta'] <= 2)
qq_plot(test[,'beta'][indices],pred_beta[indices],'beta (0 <-> 2)')
res4 <- perf_metrics(test[,'beta'][indices],pred_beta[indices])

list_of_vectors <- list(res1,res2,res3,res4)
df <- do.call(rbind,list_of_vectors)
row.names(df) <- c("alpha","alpha 0 <-> 2","beta","beta 0 <-> 2")
df <- t(df)
row.names(df) <- c("Mean Absolute Error","RMSE","Spearman Rho","R^2")
print(textplot(round(df,3)))

# quantiles
quantiles <- quantreg(Multivar(alpha, beta) ~., data = train,newdata=test,mtry =o$optimal['mtry'],nodesize= o$optimal['nodesize'])

plot.quantreg(quantiles,m.target='alpha')
plot.quantreg(quantiles,m.target='beta')

# Partial Dependence plots
variables = c('LET', 'IrradiationConditions', 'RadiationType', 'nonDSBClusters', 'DSBs', 'CellClass', 'CellCycle')
for (variable in variables){
	par(mfrow=c(1,1))
	plot.variable(model,
		      variable,
		      partial = T,
		      m.target = "alpha",
		      )
}

variables = c('Tissue','CellLine')
for (variable in variables){
	partial.obj <- partial(model,
			       partial.xvar = variable,
			       partial.values = model$xvar[,variable]
			       )
	pdta1 <- get.partial.plot.data(partial.obj, m.target = "alpha")
	n <- 5  # number of max elements to select
	numbers <-  unlist(pdta1$yhat)
	names <- unlist(pdta1$x)
	df <- data.frame(numbers = numbers,names=names )
	print(dim(df))
	top_n <- df %>% top_n(5,numbers)
	top_names <- top_n[,"names"]
	df$top_names <- ifelse(df$names %in% top_names,as.character(df$names),"")	

	p <- barplot(pdta1$yhat, 
		     names.arg="", 
		     ylim = c(0.95*min(pdta1$yhat),1.1*max(pdta1$yhat)),
		     ylab = expression(hat(y)),
		     xlab = variable,
		     xpd = F
		     )
	box()
	text(x = p, y = pdta1$yhat + 0.01,labels=df$top_names,srt = -45)

}
print(paste('ratio: ',train_test_ratio))
print('dim(train):')
print(dim(train))
print('dim(test)')
print(dim(test))
print('optimal parameters')
print(o$optimal)
print(quantiles)
sessionInfo()
q()
