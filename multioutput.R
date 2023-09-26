#Attached Packages:
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
library(randomForestSRC)
qq_plot <- function(actual,pred,target){
	na_indeces <- which(is.na(actual))
	actual <- actual[-na_indeces]
	pred <- pred[-na_indeces]

	print(paste("RMSE: ",target))
	print(sqrt(mean((actual - pred)^2)))

	min = min(min(actual),min(pred))
	max = max(max(actual),max(pred))
	plot(x = pred, y = actual,
	     xlab = 'Predicted ',
	     ylab = 'Actual',
	     ylim = c(min,max),
	     xlim = c(min,max),
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
	# Create the legend text
	legend_text <- paste("y = ", round(coefficients[2], 2), "x + ", round(coefficients[1], 2), "\nR^2 = ", round(r_squared, 2))
	# Add the legend to the plot
	legend("topleft", legend = legend_text, col = "red", lty = 1,bg='grey')
}

pdf('multioutput.pdf',width=12,height=9)
data <- read.csv('data_29_08_23.csv')
########################## Data Preparation  ###########################
print(colnames(data))
hist(data$alpha)
hist(data$beta)

# Removal of mostly empty features
data$RBE <-NULL
data$Energy <-NULL
data$DoseRate <-NULL
data <- data %>% mutate(CellLine = str_trim(CellLine))
data <- data %>% mutate(CellLine = toupper(CellLine))
data <- data %>% mutate(CellClass = str_trim(CellClass))
data <- data %>% mutate(IrradiationConditions = str_trim(IrradiationConditions))
data <- data %>% mutate(Tissue = str_trim(Tissue))
data <- data %>% mutate(Tissue = toupper(Tissue))
data <- subset(data, LET < 1000)

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
## 90% of the sample size
smp_size <- floor(0.9 * nrow(data))
## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(data)), size = smp_size)

train <- data[train_ind, ]
test <- data[-train_ind, ]

o <- tune(Multivar(alpha, beta) ~., train)
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
par(mfrow=c(1,1))
qq_plot(test[,'beta'],pred_alpha,'beta')

# quantiles
o <- quantreg(Multivar(alpha, beta) ~., data = train,newdata=test,mtry =o$optimal['mtry'],nodesize= o$optimal['nodesize'])

plot.quantreg(o,m.target='alpha')
plot.quantreg(o,m.target='beta')

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
sessionInfo()
q()
