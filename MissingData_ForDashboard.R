# Set-up ----------------------------------------------------------------
# Packages
library(dplyr) # For data manipulation
library(magrittr) # For piping
library(ggplot2) # For plots
library(grid) # To organize multiple subplots in one large plot
library(gridExtra) # To organize multiple subplots in one large plot
library(msm) # For use of rtnorm function to simulate truncated normal data
library(mice) # For imputation
library(miceadds) # For multiple imputation correlation and descriptives

# Set seed
set.seed(123)

# Simulate data -----------------------------------------------------------
# Induce correlation
n     <- 100                                   # length of vector
rho   <- 0.6                                  # desired correlation = cos(angle)
theta <- acos(rho)                            # corresponding angle
weight <- rtnorm(n = n, mean = 72.55, sd = 11.86, lower = 50, upper = 100) %>% round(2)   # Truncated normal data
height <- rtnorm(n = n, mean = 1.70, sd = 0.06, lower = 1.55, upper = 1.95) %>% round(4) # Truncated normal data
X     <- cbind(weight, height)                # matrix
Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)
Id   <- diag(n)                               # identity matrix
Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
height <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]# final new vector
cor(height, weight)                           # check correlation = rho
height <- height + 1.70                       # Add original mean of height to scaled new variable

# Inspect scatterplot
plot(height, weight)
# Add to dataframe
df <- data.frame(height, weight)
# Calculate BMI
df$bmi <- df$weight / df$height^2

# #remove unnecessary variables
remove(n, rho, theta, weight, height, X, Xctr, Id, Q, P, x2o, Xc2, Y)

# Induce NA's -------------------------------------------------------------
# Create new dataframe
df_NA <- df
# Induced biased missing pattern, by deleting half of all weight values that are associated with a BMI > 25 (overweight)
for (i in which(df$bmi > 25)) {
  if ((i %% 2) == 0) {
    df_NA[i, "weight"] <- NA
  } else {
    df_NA[i, "weight"] <- df_NA[i, "weight"]
  }
}
# Also delete all BMI values that have NA's in weight
df_NA[which(is.na(df_NA$weight)), "bmi"] <- NA

# #remove unnecessary variables
remove(i)

# Complete-case analysis --------------------------------------------------
# Create new dataframe with only complete observations
df_CCA <- na.omit(df_NA)

## Plots
data <- df_CCA

# Scatterplot
f1_CCA <- ggplot(data, aes(weight, height)) + 
  geom_point(shape=1, size=1.5, stroke=1.5, colour = "black") + # Use hollow circles, and make bigger and thicker
  geom_smooth(method=lm,                            # Add linear regression line
              se=FALSE,                             # Don't add shaded confidence region
              color="black") +                      # Make regression line black
  labs(x = "Weight (kg)",                           # Change x-axis label
       y = "Height (m)") +                          # Change y-axis label
  theme_classic() +                                 # Set theme to classic
  theme(axis.title = element_text(size = 14),       # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black") # Increase size of axis labels and make black
  )
# Barplot weight
f2_CCA <- ggplot(data, aes(weight)) + 
  geom_histogram(fill="white", bins = 40, colour = "black") + # Use hollow circles, and make bigger and thicker
  labs(x = "Weight (kg)",                   # Change x-axis label
       y = "Count") +                          # Change y-axis label
  theme_classic() +                         # Set theme to classic
  theme(axis.title = element_text(size = 14), # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black") # Increase size of axis labels and make black
  )
# Barplot BMI
f3_CCA <- ggplot(data, aes(bmi)) + 
  geom_histogram(fill="white", bins = 40, colour = "black") + # Use hollow circles, and make bigger and thicker
  labs(x = "Body-mass index",                   # Change x-axis label
       y = "Count") +                          # Change y-axis label
  theme_classic() +                         # Set theme to classic
  theme(axis.title = element_text(size = 14), # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black") # Increase size of axis labels and make black
  )
# #remove unnecessary variables
remove(data)


# Mean imputation ---------------------------------------------------------
# Create new dataframe
df_MeanI <- df_NA
# Replace missing values in weight with new mean weight
df_MeanI[is.na(df_MeanI$weight), "weight"] <- mean(df_MeanI$weight, na.rm = TRUE)
# Calculate new BMI based on imputed data
df_MeanI$bmi <- df_MeanI$weight / df_MeanI$height^2

## Plots
# Add missingness indicator for plotting
df_MeanI$miss <- as.factor(ifelse(is.na(df_NA$weight), 1, 0))
data <- df_MeanI

# Scatterplot
f1_MeanI <- ggplot(data, aes(weight, height, colour = miss)) + # Use missingness indicator in different colour
        geom_point(shape=1, size=1.5, stroke=1.5) +            # Use hollow circles, and make bigger and thicker
        geom_smooth(method=lm,                             # Add linear regression line
                    se=FALSE,                             # Don't add shaded confidence region
                    color="black") +                      # Make regression line black
        labs(#title = "Mean imputation",                    # Add plot title
              x = "Weight (kg)",                            # Change x-axis label
              y = "Height (m)",                             # Change y-axis label
              colour = "Imputed") +                       # Change lagend title
        theme_classic() +                                   # Set theme to classic
        theme(axis.title = element_text(size = 14),         # Increase size of axis title
              axis.text = element_text(size = 12, colour = "black"), # Increase size of axis labels and make black
              legend.title = element_text(size = 12), # Increase size of legend title
              legend.text = element_text(size = 10), # Decrease size of legend labels
              legend.position = "top" # Place legend on the top of the figure
              ) +
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) # Change labels and colour of legend labels
# Barplot weight
f2_MeanI <- ggplot(data, aes(weight, colour = miss)) + # Use missingness indicator in different colour
        geom_histogram(fill="white", bins = 40) + # Use hollow circles, and make bigger and thicker
        labs(x = "Weight (kg)",                   # Change x-axis label
            y = "Count",                          # Change y-axis label
            colour = "Imputed") +               # Change lagend title
        theme_classic() +                         # Set theme to classic
        theme(axis.title = element_text(size = 14), # Increase size of axis title
              axis.text = element_text(size = 12, colour = "black"), # Increase size of axis labels and make black
              legend.title = element_text(size = 12), # Increase size of legend title
              legend.text = element_text(size = 10), # Decrease size of legend labels
              legend.position = "top" # Place legend on the top of the figure
        ) +
        scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) # Change labels and colour of legend labels
# Barplot BMI
f3_MeanI <- ggplot(data, aes(bmi, colour = miss)) + # Use missingness indicator in different colour
  geom_histogram(fill="white", bins = 40) + # Use hollow circles, and make bigger and thicker
  labs(x = "Body-mass index",                   # Change x-axis label
       y = "Count",                          # Change y-axis label
       colour = "Imputed") +               # Change lagend title
  theme_classic() +                         # Set theme to classic
  theme(axis.title = element_text(size = 14), # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black"), # Increase size of axis labels and make black
        legend.title = element_text(size = 12), # Increase size of legend title
        legend.text = element_text(size = 10), # Decrease size of legend labels
        legend.position = "top" # Place legend on the top of the figure
  ) +
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) # Change labels and colour of legend labels

# #remove unnecessary variables
remove(data)


# Regression imputation ---------------------------------------------------

# Use package 'mice'
# Create a default predictor matrix.
pred <- make.predictorMatrix(df_NA) # (rows are predicted by the columns.)
# Do not predict BMI
pred[c("height", "weight"), "bmi"] <- 0

## Get method matrix
meth <- make.method(df_NA, default = "norm.predict")
# Use passive imputation for BMI
meth["bmi"]<- "~I(weight / height^2)"

# Perform regression imputation with mice package
df_RI_imp <- mice(data = df_NA,
                  predictorMatrix = pred, # Custom predictor matrix
                  method = meth, # custom method matrix
                  m = 1, # Only impute one dataset
                  maxit = 10, # Number of iterations
                  seed = 123, # Seed
                  printFlag = FALSE) # Don't print all output 
# Store data
df_RI <- complete(df_RI_imp)

## Alternative way: use lm and predict functions
## Create new dataframe
#df_RI <- df_NA
## Calculate regression model based on complete cases
#fit <- lm(weight ~ height, data = df_CCA)
## Get the unused height values
#new <- data.frame(height = df_NA[which(is.na(df_NA$weight)), "height"])
## Predict new weight values from unused height values and put them in NA weight spots
#df_RI[which(is.na(df_RI$weight)), "weight"] <- predict(fit, newdata = new)
## Calculate new BMI
#df_RI$bmi <- df_RI$weight / df_RI$height^2

## Plots
# Add missingness indicator for plotting
df_RI$miss <- as.factor(ifelse(is.na(df_NA$weight), 1, 0))
data <- df_RI

# Scatterplot
f1_RI <- ggplot(data, aes(weight, height, colour = miss)) + # Use missingness indicator in different colour
  geom_point(shape=1, size=1.5, stroke=1.5) +            # Use hollow circles, and make bigger and thicker
  geom_smooth(method=lm,                             # Add linear regression line
              se=FALSE,                             # Don't add shaded confidence region
              color="black") +                      # Make regression line black
  labs(#title = "Mean imputation",                    # Add plot title
    x = "Weight (kg)",                            # Change x-axis label
    y = "Height (m)",                             # Change y-axis label
    colour = "Imputed") +                       # Change lagend title
  theme_classic() +                                   # Set theme to classic
  theme(axis.title = element_text(size = 14),         # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black"), # Increase size of axis labels and make black
        legend.title = element_text(size = 12), # Increase size of legend title
        legend.text = element_text(size = 10), # Decrease size of legend labels
        legend.position = "top" # Place legend on the top of the figure
  ) +
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) # Change labels and colour of legend labels
# Barplot weight
f2_RI <- ggplot(data, aes(weight, colour = miss)) + # Use missingness indicator in different colour
  geom_histogram(fill="white", bins = 40) + # Use hollow circles, and make bigger and thicker
  labs(x = "Weight (kg)",                   # Change x-axis label
       y = "Count",                          # Change y-axis label
       colour = "Imputed") +               # Change lagend title
  theme_classic() +                         # Set theme to classic
  theme(axis.title = element_text(size = 14), # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black"), # Increase size of axis labels and make black
        legend.title = element_text(size = 12), # Increase size of legend title
        legend.text = element_text(size = 10), # Decrease size of legend labels
        legend.position = "top" # Place legend on the top of the figure
  ) +
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) # Change labels and colour of legend labels
# Barplot BMI
f3_RI <- ggplot(data, aes(bmi, colour = miss)) + # Use missingness indicator in different colour
  geom_histogram(fill="white", bins = 40) + # Use hollow circles, and make bigger and thicker
  labs(x = "Body-mass index",                   # Change x-axis label
       y = "Count",                          # Change y-axis label
       colour = "Imputed") +               # Change lagend title
  theme_classic() +                         # Set theme to classic
  theme(axis.title = element_text(size = 14), # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black"), # Increase size of axis labels and make black
        legend.title = element_text(size = 12), # Increase size of legend title
        legend.text = element_text(size = 10), # Decrease size of legend labels
        legend.position = "top" # Place legend on the top of the figure
  ) +
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) # Change labels and colour of legend labels

# #remove unnecessary variables
remove(meth, pred, data)
##remove(fit, new)


# Stochastic regression imputation ---------------------------------------------------

# Use package 'mice'
# Create a default predictor matrix.
pred <- make.predictorMatrix(df_NA) # (rows are predicted by the columns.)
# Do not predict BMI
pred[c("height", "weight"), "bmi"] <- 0

## Get method matrix
meth <- make.method(df_NA, default = "norm.nob")
# Use passive imputation for BMI
meth["bmi"]<- "~I(weight / height^2)"

# Perform regression imputation with mice package
df_SRI_imp <- mice(data = df_NA, 
                   predictorMatrix = pred,
                   method = meth, # custom method matrix
                   m = 1, # Only impute one dataset
                   maxit = 10, # Number of iterations
                   seed = 123, # Seed
                   printFlag = FALSE) # Do not print all output
# Store data
df_SRI <- complete(df_SRI_imp)

## Plots
# Add missingness indicator for plotting
df_SRI$miss <- as.factor(ifelse(is.na(df_NA$weight), 1, 0))
data <- df_SRI
# Scatterplot
f1_SRI <- ggplot(data, aes(weight, height, colour = miss)) + # Use missingness indicator in different colour
  geom_point(shape=1, size=1.5, stroke=1.5) +            # Use hollow circles, and make bigger and thicker
  geom_smooth(method=lm,                             # Add linear regression line
              se=FALSE,                             # Don't add shaded confidence region
              color="black") +                      # Make regression line black
  labs(#title = "Mean imputation",                    # Add plot title
    x = "Weight (kg)",                            # Change x-axis label
    y = "Height (m)",                             # Change y-axis label
    colour = "Imputed") +                       # Change lagend title
  theme_classic() +                                   # Set theme to classic
  theme(axis.title = element_text(size = 14),         # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black"), # Increase size of axis labels and make black
        legend.title = element_text(size = 12), # Increase size of legend title
        legend.text = element_text(size = 10), # Decrease size of legend labels
        legend.position = "top" # Place legend on the top of the figure
  ) +
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) # Change labels and colour of legend labels
# Barplot weight
f2_SRI <- ggplot(data, aes(weight, colour = miss)) + # Use missingness indicator in different colour
  geom_histogram(fill="white", bins = 40) + # Use hollow circles, and make bigger and thicker
  labs(x = "Weight (kg)",                   # Change x-axis label
       y = "Count",                          # Change y-axis label
       colour = "Imputed") +               # Change lagend title
  theme_classic() +                         # Set theme to classic
  theme(axis.title = element_text(size = 14), # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black"), # Increase size of axis labels and make black
        legend.title = element_text(size = 12), # Increase size of legend title
        legend.text = element_text(size = 10), # Decrease size of legend labels
        legend.position = "top" # Place legend on the top of the figure
  ) +
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) # Change labels and colour of legend labels
# Barplot BMI
f3_SRI <- ggplot(data, aes(bmi, colour = miss)) + # Use missingness indicator in different colour
  geom_histogram(fill="white", bins = 40) + # Use hollow circles, and make bigger and thicker
  labs(x = "Body-mass index",                   # Change x-axis label
       y = "Count",                          # Change y-axis label
       colour = "Imputed") +               # Change lagend title
  theme_classic() +                         # Set theme to classic
  theme(axis.title = element_text(size = 14), # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black"), # Increase size of axis labels and make black
        legend.title = element_text(size = 12), # Increase size of legend title
        legend.text = element_text(size = 10), # Decrease size of legend labels
        legend.position = "top" # Place legend on the top of the figure
  ) +
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) # Change labels and colour of legend labels

# #remove unnecessary variables
remove(meth, pred, df_SRI_imp, data)


# Multiple imputation ---------------------------------------------------

# Create a default predictor matrix.
pred <- make.predictorMatrix(df_NA) # (rows are predicted by the columns.)
# Do not predict BMI
pred[c("height", "weight"), "bmi"] <- 0

## Get method matrix
meth <- make.method(df_NA)
# Use passive imputation for BMI
meth["bmi"]<- "~I(weight / height^2)"

# Perform imputation with mice package
df_MI_imp <- mice(data = df_NA, 
                  predictorMatrix = pred, # custom predictor matrix
                  method = meth, # custom method matrix
                  m = I(100-sum(is.na(df_NA))*100/nrow(df_NA)), # Impute as many datasets as percentage of missing data
                  maxit = 10, # Number of iterations
                  seed = 123, # Seed
                  print = FALSE) # Do not print all computational history
# Store data
df_MI <- complete(df_MI_imp)

## Check imputation quality
# Check multicolinearity
df_MI$loggedEvents # -> If NULL, no multicolinearity
# Plot the iterations and check convergence.
plot(df_MI_imp, layout = c(2, 2)) # The passive imputations plots will be straight lines, all other variables need to have overlapping squiggly lines that are close together on the right side.
# Boxplots
bwplot(df_MI_imp, cex = 0.5, layout = c(1, 3)) # See boxplots of the imputated data for each set and variable to compare with original data.
# Stripplots
stripplot(df_MI_imp, cex = 0.5, layout = c(1, 3))# See stripplots of the imputated data for each set and variable to compare with original data.
# Check out the imputed data
head(df_MI_imp$imp$weight) # Do values look plausible?

## Plots
# Add missingness indicator for plotting
df_MI$miss <- as.factor(ifelse(is.na(df_NA$weight), 1, 0))
data <- df_MI
# Scatterplot
jitter <- position_jitter(width = 0.05, height = 0.05, seed = 123) # Add jitter to point positions
f1_MI <- ggplot(data, aes(weight, height, colour = miss)) + # Use missingness indicator in different colour
  geom_point(shape=1, size=1.5, stroke=1.5, position = jitter) +            # Use hollow circles, and make bigger and thicker
  geom_smooth(method=lm,                             # Add linear regression line
              se=FALSE,                             # Don't add shaded confidence region
              color="black") +                      # Make regression line black
  labs(#title = "Mean imputation",                    # Add plot title
    x = "Weight (kg)",                            # Change x-axis label
    y = "Height (m)",                             # Change y-axis label
    colour = "Imputed") +                       # Change lagend title
  theme_classic() +                                   # Set theme to classic
  theme(axis.title = element_text(size = 14),         # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black"), # Increase size of axis labels and make black
        legend.title = element_text(size = 12), # Increase size of legend title
        legend.text = element_text(size = 10), # Decrease size of legend labels
        legend.position = "top" # Place legend on the top of the figure
  ) +
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) # Change labels and colour of legend labels
# Barplot weight
f2_MI <- ggplot(data, aes(weight, colour = miss)) + # Use missingness indicator in different colour
  geom_histogram(fill="white", bins = 40) + # Use hollow circles, and make bigger and thicker
  labs(x = "Weight (kg)",                   # Change x-axis label
       y = "Count",                          # Change y-axis label
       colour = "Imputed") +               # Change lagend title
  theme_classic() +                         # Set theme to classic
  theme(axis.title = element_text(size = 14), # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black"), # Increase size of axis labels and make black
        legend.title = element_text(size = 12), # Increase size of legend title
        legend.text = element_text(size = 10), # Decrease size of legend labels
        legend.position = "top" # Place legend on the top of the figure
  ) +
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) # Change labels and colour of legend labels
# Barplot BMI
f3_MI <- ggplot(data, aes(bmi, colour = miss)) + # Use missingness indicator in different colour
  geom_histogram(fill="white", bins = 40) + # Use hollow circles, and make bigger and thicker
  labs(x = "Body-mass index",                   # Change x-axis label
       y = "Count",                          # Change y-axis label
       colour = "Imputed") +               # Change lagend title
  theme_classic() +                         # Set theme to classic
  theme(axis.title = element_text(size = 14), # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black"), # Increase size of axis labels and make black
        legend.title = element_text(size = 12), # Increase size of legend title
        legend.text = element_text(size = 10), # Decrease size of legend labels
        legend.position = "top" # Place legend on the top of the figure
  ) +
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) # Change labels and colour of legend labels

# #remove unnecessary variables
remove(meth, pred, data, jitter)


# Ultimate comparison -----------------------------------------------------

## Difference in mean BMI
# Calculate means
df_mean <- data.frame(
  Original = c( I(mean(df$weight)), I(mean(df$bmi)) ) %>% round(2),
  Complete_Cases_Analysis = c( I(mean(df_CCA$weight)), I(mean(df_CCA$bmi)) ) %>% round(2), 
  Mean_Imputation = c( I(mean(df_MeanI$weight)), I(mean(df_MeanI$bmi)) ) %>% round(2), 
  Regression_Imputation = c( I(mean(df_RI$weight)), I(mean(df_RI$bmi)) ) %>% round(2), 
  Stochastic_Regression_Imputation = c( I(mean(df_SRI$weight)), I(mean(df_SRI$bmi)) ) %>% round(2),
  Multiple_Imputation = c( I(with(df_MI_imp, expr=c(mean(weight) ) ) %>% withPool_MI()) %>% round(2), 
                           I(with(df_MI_imp, expr=c(mean(bmi) ) ) %>% withPool_MI()) ) %>% round(2),
  row.names = c("Mean_Weight", "Mean_BMI")
  )

## Calculate difference between original mean and imputed mean
# Weight
for (i in 1:ncol(df_mean)) {
  df_mean["Residual_Weight", i] <- df_mean[1, "Original"] - df_mean[1, i]
}
# BMI
for (i in 1:ncol(df_mean)) {
  df_mean["Residual_BMI", i] <- df_mean[2, "Original"] - df_mean[2, i]
}
# Calculate winner
df_mean["Mean_Weight", "Most_Accurate"] <- colnames(df_mean["Residual_Weight", -1])[which.min(abs(df_mean["Residual_Weight", -1]))]
df_mean["Mean_BMI", "Most_Accurate"] <- colnames(df_mean[, -c(which(colnames(df_mean) == "Original"), 
                                                           which(colnames(df_mean) == "Most_Accurate"))])[which.min(abs(df_mean["Residual_BMI", -c(which(colnames(df_mean) == "Original"), 
                                                                                                                                            which(colnames(df_mean) == "Most_Accurate"))]))]
# Remove difference calculation
df_mean <- df_mean[-c(3:nrow(df_mean)), ]

## Difference in standard error of BMI
# Calculate SDs
df_sd <- data.frame(
  Original = c( I(sd(df$weight)), I(sd(df$bmi)) ) %>% round(2),
  Complete_Cases_Analysis = c( I(sd(df_CCA$weight)), I(sd(df_CCA$bmi)) ) %>% round(2), 
  Mean_Imputation = c( I(sd(df_MeanI$weight)), I(sd(df_MeanI$bmi)) ) %>% round(2), 
  Regression_Imputation = c( I(sd(df_RI$weight)), I(sd(df_RI$bmi)) ) %>% round(2), 
  Stochastic_Regression_Imputation = c( I(sd(df_SRI$weight)), I(sd(df_SRI$bmi)) ) %>% round(2),
  Multiple_Imputation = c( I(with(df_MI_imp, expr=c(sd(weight) ) ) %>% withPool_MI()) %>% round(2), 
                           I(with(df_MI_imp, expr=c(sd(bmi) ) ) %>% withPool_MI()) ) %>% round(2),
  row.names = c("sd_Weight", "sd_BMI")
)

## Calculate difference between original sd and imputed sd
# Weight
for (i in 1:ncol(df_sd)) {
  df_sd["Residual_Weight", i] <- df_sd[1, "Original"] - df_sd[1, i]
}
# BMI
for (i in 1:ncol(df_sd)) {
  df_sd["Residual_BMI", i] <- df_sd[2, "Original"] - df_sd[2, i]
}

# Calculate winner
df_sd["sd_Weight", "Most_Accurate"] <- colnames(df_sd["Residual_Weight", -1])[which.min(abs(df_sd["Residual_Weight", -1]))]
df_sd["sd_BMI", "Most_Accurate"] <- colnames(df_sd[, -c(which(colnames(df_sd) == "Original"), 
                                                           which(colnames(df_sd) == "Most_Accurate"))])[which.min(abs(df_sd["Residual_BMI", -c(which(colnames(df_sd) == "Original"), 
                                                                                                                                            which(colnames(df_sd) == "Most_Accurate"))]))]
# Remove difference calculation
df_sd <- df_sd[-c(3:nrow(df_sd)), ]

## Difference in correlations between height and weight
# Calculate the correlation for multiple imputation
corr <- micombine.cor(mi.res = df_MI_imp, variables = c("weight", "height")) %>% attr("r_matrix")
corr <- corr[lower.tri(corr)] %>% round(2)
# Calculate other correlations
df_cor <- data.frame(
  Original = cor(df$weight, df$height) %>% round(2),
  Complete_Cases_Analysis = cor(df_CCA$weight, df_CCA$height) %>% round(2),
  Mean_Imputation = cor(df_MeanI$weight, df_MeanI$height) %>% round(2), 
  Regression_Imputation = cor(df_RI$weight, df_RI$height) %>% round(2), 
  Stochastic_Regression_Imputation = cor(df_SRI$weight, df_SRI$height) %>% round(2),
  Multiple_Imputation = corr,
  row.names = c("Correlation_Height_Weight")
)

# Calculate difference between original correlation and imputed correlation
for (i in 1:ncol(df_cor)) {
  df_cor["Residual_Correlation", i] <- df_cor[1, "Original"] - df_cor[1, i]
}

# Calculate Winner
df_cor["Correlation_Height_Weight", "Most_Accurate"] <- colnames(df_cor["Residual_Correlation", -1])[which.min(abs(df_cor["Residual_Correlation", -1]))]

# Remove difference calculation
df_cor <- df_cor[-c(2:nrow(df_cor)), ]

# #remove unnecessary variables
remove(corr, i)
