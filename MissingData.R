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

# Remove unnecessary variables
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
df_NA[ which(is.na(df_NA$weight)), "bmi"] <- NA

# Remove unnecessary variables
remove(i)

# Complete-case analysis --------------------------------------------------
# Create new dataframe with only complete observations
df_CCA <- na.omit(df_NA)

## Plots
data <- df_CCA

# Scatterplot
corr_orig <- cor(df$weight, df$height) %>% round(2) # Calculate original correlation between height and weight
grob_orig <- grobTree(textGrob(bquote("(Orig." ~ italic("r") ~ "=" ~ italic(.(corr_orig)) ~ ")"), hjust=0, x=0.05,  y=0.90, # Put into grob
                          gp=gpar(col="black", fontsize=9)))
corr <- cor(data$weight, data$height) %>% round(2) # Calculate correlation between height and weight
grob <- grobTree(textGrob(bquote(italic("r") ~ "=" ~ italic(.(corr))), hjust=0, x=0.05,  y=0.95, # Put into grob
                          gp=gpar(col="black", fontsize=9)))
f1 <- ggplot(data, aes(weight, height)) + 
  geom_point(shape=1, size=1.5, stroke=1.5, colour = "black") + # Use hollow circles, and make bigger and thicker
  geom_smooth(method=lm,                            # Add linear regression line
              se=FALSE,                             # Don't add shaded confidence region
              color="black") +                      # Make regression line black
  labs(x = "Weight (kg)",                           # Change x-axis label
       y = "Height (m)") +                          # Change y-axis label
  theme_classic() +                                 # Set theme to classic
  theme(axis.title = element_text(size = 14),       # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black") # Increase size of axis labels and make black
  ) +
  annotation_custom(grob) +                          # Add correlation text
  annotation_custom(grob_orig)                       # Add original correlation text
# Barplot
mean_weight_orig <- mean(df$weight) %>% round(1) # Calculate mean weight
sd_weight_orig <- sd(df$weight) %>% round(1) # Calculate sd of weight
grob_orig <- grobTree(textGrob(bquote("(Orig. M =" ~ .(mean_weight_orig) ~ "(SD =" ~ .(sd_weight_orig) ~ "))"), hjust = 0, x=0.55,  y=0.90, # Put into grob
                               gp=gpar(col="black", fontsize=9)))
mean_weight <- mean(data$weight) %>% round(1) # Calculate mean weight
sd_weight <- sd(data$weight) %>% round(1) # Calculate sd of weight
grob <- grobTree(textGrob(bquote("M =" ~ .(mean_weight) ~ "(SD =" ~ .(sd_weight) ~ ")"), hjust = 0, x=0.55,  y=0.95, # Put into grob
                          gp=gpar(col="black", fontsize=9)))
f2 <- ggplot(data, aes(weight)) + 
  geom_histogram(fill="white", bins = 40, colour = "black") + # Use hollow circles, and make bigger and thicker
  labs(x = "Weight (kg)",                   # Change x-axis label
       y = "Count") +                          # Change y-axis label
  theme_classic() +                         # Set theme to classic
  theme(axis.title = element_text(size = 14), # Increase size of axis title
        axis.text = element_text(size = 12, colour = "black") # Increase size of axis labels and make black
  ) +
  annotation_custom(grob) +                          # Add mean text
  annotation_custom(grob_orig)                       # Add original mean text
## Plot all subfigures simultaneously
tiff(filename = "CompleteCases.tiff",
     width = 2500,
     height = 1000,
     family = "serif",
     res = 300)
# Plot side-by-side
grid.arrange(f1, f2, nrow = 1,
             top = textGrob("Complete cases", gp = gpar(fontsize = 20))) # Add title
# Print figure
dev.off()

# Remove unnecessary variables
remove(f1, f2, data, grob, grob_orig, corr, corr_orig, mean_weight, mean_weight_orig, sd_weight, sd_weight_orig)


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
corr_orig <- cor(df$weight, df$height) %>% round(2) # Calculate original correlation between height and weight
grob_orig <- grobTree(textGrob(bquote("(Orig." ~ italic("r") ~ "=" ~ italic(.(corr_orig)) ~ ")"), hjust=0, x=0.05,  y=0.90, # Put into grob
                               gp=gpar(col="black", fontsize=9)))
corr <- cor(data$weight, data$height) %>% round(2) # Calculate correlation between height and weight
grob <- grobTree(textGrob(bquote(italic("r") ~ "=" ~ italic(.(corr))), hjust=0, x=0.05,  y=0.95, # Put into grob
                          gp=gpar(col="black", fontsize=9)))
f1 <- ggplot(data, aes(weight, height, colour = miss)) + # Use missingness indicator in different colour
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
        scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) + # Change labels and colour of legend labels
  annotation_custom(grob) +                          # Add correlation text
  annotation_custom(grob_orig)                       # Add original correlation text
# Barplot
mean_weight_orig <- mean(df$weight) %>% round(1) # Calculate mean weight
sd_weight_orig <- sd(df$weight) %>% round(1) # Calculate sd of weight
grob_orig <- grobTree(textGrob(bquote("(Orig. M =" ~ .(mean_weight_orig) ~ "(SD =" ~ .(sd_weight_orig) ~ "))"), hjust = 0, x=0.55,  y=0.90, # Put into grob
                               gp=gpar(col="black", fontsize=9)))
mean_weight <- mean(data$weight) %>% round(1) # Calculate mean weight
sd_weight <- sd(data$weight) %>% round(1) # Calculate sd of weight
grob <- grobTree(textGrob(bquote("M =" ~ .(mean_weight) ~ "(SD =" ~ .(sd_weight) ~ ")"), hjust = 0, x=0.55,  y=0.95, # Put into grob
                          gp=gpar(col="black", fontsize=9)))
f2 <- ggplot(data, aes(weight, colour = miss)) + # Use missingness indicator in different colour
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
        scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) + # Change labels and colour of legend labels
  annotation_custom(grob) +                          # Add mean text
  annotation_custom(grob_orig)                       # Add original mean text
## Plot all subfigures simultaneously
tiff(filename = "MeanImputation.tiff",
     width = 2500,
     height = 1000,
     family = "serif",
     res = 300)
# Plot side-by-side
grid.arrange(f1, f2, nrow = 1,
             top = textGrob("Mean imputation", gp = gpar(fontsize = 20))) # Add title
# Print figure
dev.off()

# Remove unnecessary variables
remove(f1, f2, data, grob, grob_orig, corr, corr_orig, mean_weight, mean_weight_orig, sd_weight, sd_weight_orig)


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
corr_orig <- cor(df$weight, df$height) %>% round(2) # Calculate original correlation between height and weight
grob_orig <- grobTree(textGrob(bquote("(Orig." ~ italic("r") ~ "=" ~ italic(.(corr_orig)) ~ ")"), hjust=0, x=0.05,  y=0.90, # Put into grob
                               gp=gpar(col="black", fontsize=9)))
corr <- cor(data$weight, data$height) %>% round(2) # Calculate correlation between height and weight
grob <- grobTree(textGrob(bquote(italic("r") ~ "=" ~ italic(.(corr))), hjust=0, x=0.05,  y=0.95, # Put into grob
                          gp=gpar(col="black", fontsize=9)))
f1 <- ggplot(data, aes(weight, height, colour = miss)) + # Use missingness indicator in different colour
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
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) + # Change labels and colour of legend labels
  annotation_custom(grob) +                          # Add correlation text
  annotation_custom(grob_orig)                       # Add original correlation text
# Barplot
mean_weight_orig <- mean(df$weight) %>% round(1) # Calculate mean weight
sd_weight_orig <- sd(df$weight) %>% round(1) # Calculate sd of weight
grob_orig <- grobTree(textGrob(bquote("(Orig. M =" ~ .(mean_weight_orig) ~ "(SD =" ~ .(sd_weight_orig) ~ "))"), hjust = 0, x=0.55,  y=0.90, # Put into grob
                               gp=gpar(col="black", fontsize=9)))
mean_weight <- mean(data$weight) %>% round(1) # Calculate mean weight
sd_weight <- sd(data$weight) %>% round(1) # Calculate sd of weight
grob <- grobTree(textGrob(bquote("M =" ~ .(mean_weight) ~ "(SD =" ~ .(sd_weight) ~ ")"), hjust = 0, x=0.55,  y=0.95, # Put into grob
                          gp=gpar(col="black", fontsize=9)))
f2 <- ggplot(data, aes(weight, colour = miss)) + # Use missingness indicator in different colour
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
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) +# Change labels and colour of legend labels
  annotation_custom(grob) +                          # Add mean text
  annotation_custom(grob_orig)                       # Add original mean text
## Plot all subfigures simultaneously
tiff(filename = "RegressionImputation.tiff",
     width = 2500,
     height = 1000,
     family = "serif",
     res = 300)
# Plot side-by-side
grid.arrange(f1, f2, nrow = 1,
             top = textGrob("Regression imputation", gp = gpar(fontsize = 20))) # Add title
# Print figure
dev.off()

# Remove unnecessary variables
remove(meth, pred, f1, f2, data, df_RI_imp, grob, grob_orig, corr, corr_orig, mean_weight, mean_weight_orig, sd_weight, sd_weight_orig)
#remove(fit, new)


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
corr_orig <- cor(df$weight, df$height) %>% round(2) # Calculate original correlation between height and weight
grob_orig <- grobTree(textGrob(bquote("(Orig." ~ italic("r") ~ "=" ~ italic(.(corr_orig)) ~ ")"), hjust=0, x=0.05,  y=0.90, # Put into grob
                               gp=gpar(col="black", fontsize=9)))
corr <- cor(data$weight, data$height) %>% round(2) # Calculate correlation between height and weight
grob <- grobTree(textGrob(bquote(italic("r") ~ "=" ~ italic(.(corr))), hjust=0, x=0.05,  y=0.95, # Put into grob
                          gp=gpar(col="black", fontsize=9)))
f1 <- ggplot(data, aes(weight, height, colour = miss)) + # Use missingness indicator in different colour
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
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) + # Change labels and colour of legend labels
  annotation_custom(grob) +                          # Add correlation text
  annotation_custom(grob_orig)                       # Add original correlation text
# Barplot
mean_weight_orig <- mean(df$weight) %>% round(1) # Calculate mean weight
sd_weight_orig <- sd(df$weight) %>% round(1) # Calculate mean weight
grob_orig <- grobTree(textGrob(bquote("(Orig. M =" ~ .(mean_weight_orig) ~ "(SD =" ~ .(sd_weight_orig) ~ "))"), hjust = 0, x=0.55,  y=0.90, # Put into grob
                               gp=gpar(col="black", fontsize=9)))
mean_weight <- mean(data$weight) %>% round(1) # Calculate mean weight
sd_weight <- sd(data$weight) %>% round(1) # Calculate mean weight
grob <- grobTree(textGrob(bquote("M =" ~ .(mean_weight) ~ "(SD =" ~ .(sd_weight) ~ ")"), hjust = 0, x=0.55,  y=0.95, # Put into grob
                          gp=gpar(col="black", fontsize=9)))
f2 <- ggplot(data, aes(weight, colour = miss)) + # Use missingness indicator in different colour
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
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) + # Change labels and colour of legend labels
  annotation_custom(grob) +                          # Add mean text
  annotation_custom(grob_orig)                       # Add original mean text
## Plot all subfigures simultaneously
tiff(filename = "StochasticRegressionImputation.tiff",
     width = 2500,
     height = 1000,
     family = "serif",
     res = 300)
# Plot side-by-side
grid.arrange(f1, f2, nrow = 1,
             top = textGrob("Stochastic regression imputation", gp = gpar(fontsize = 20))) # Add title
# Print figure
dev.off()

# Remove unnecessary variables
remove(meth, pred, df_SRI_imp, f1, f2, data, grob, grob_orig, corr, corr_orig, mean_weight, mean_weight_orig, sd_weight, sd_weight_orig)


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
corr_orig <- cor(df$weight, df$height) %>% round(2) # Calculate original correlation between height and weight
grob_orig <- grobTree(textGrob(bquote("(Orig." ~ italic("r") ~ "=" ~ italic(.(corr_orig)) ~ ")"), hjust=0, x=0.05,  y=0.90, # Put into grob
                               gp=gpar(col="black", fontsize=9)))
corr <- corr <- micombine.cor(mi.res = df_MI_imp, variables = c("weight", "height")) %>% attr("r_matrix") # Calculate correlation between height and weight
corr <- corr[lower.tri(corr)] %>% round(2)
grob <- grobTree(textGrob(bquote(italic("r") ~ "=" ~ italic(.(corr))), hjust=0, x=0.05,  y=0.95, # Put into grob
                          gp=gpar(col="black", fontsize=9)))
jitter <- position_jitter(width = 0.05, height = 0.05, seed = 123) # Add jitter to point positions
f1 <- ggplot(data, aes(weight, height, colour = miss)) + # Use missingness indicator in different colour
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
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) + # Change labels and colour of legend labels
  annotation_custom(grob) +                          # Add correlation text
  annotation_custom(grob_orig)                       # Add original correlation text
# Barplot
mean_weight_orig <- mean(df$weight) %>% round(1) # Calculate mean weight
sd_weight_orig <- sd(df$weight) %>% round(1) # Calculate mean weight
grob_orig <- grobTree(textGrob(bquote("(Orig. M =" ~ .(mean_weight_orig) ~ "(SD =" ~ .(sd_weight_orig) ~ "))"), hjust = 0, x=0.55,  y=0.90, # Put into grob
                               gp=gpar(col="black", fontsize=9)))
mean_weight <- with(df_MI_imp, expr=c(mean(weight) ) ) %>% withPool_MI() %>% round(1) # Calculate mean weight
sd_weight <- with(df_MI_imp, expr=c(sd(weight) ) ) %>% withPool_MI() %>% round(1) # Calculate mean weight
grob <- grobTree(textGrob(bquote("M =" ~ .(mean_weight) ~ "(SD =" ~ .(sd_weight) ~ ")"), hjust = 0, x=0.55,  y=0.95, # Put into grob
                          gp=gpar(col="black", fontsize=9)))
f2 <- ggplot(data, aes(weight, colour = miss)) + # Use missingness indicator in different colour
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
  scale_color_manual(labels = c("No", "Yes"), values = c("black", "darkred")) + # Change labels and colour of legend labels
  annotation_custom(grob) +                          # Add mean text
  annotation_custom(grob_orig)                       # Add original mean text
## Plot all subfigures simultaneously
tiff(filename = "MultipleImputation.tiff",
     width = 2500,
     height = 1000,
     family = "serif",
     res = 300)
# Plot side-by-side
grid.arrange(f1, f2, nrow = 1,
             top = textGrob("Multiple imputation", gp = gpar(fontsize = 20))) # Add title
# Print figure
dev.off()

# Remove unnecessary variables
remove(meth, pred, f1, f2, data, jitter, grob, grob_orig, corr, corr_orig, mean_weight, mean_weight_orig, sd_weight, sd_weight_orig)


# Ultimate comparison -----------------------------------------------------

## Difference in mean BMI
# Calculate means
df_mean <- c(mean(df$bmi),
             mean(df_CCA$bmi), 
             mean(df_RI$bmi), 
             mean(df_SRI$bmi),
             with(df_MI_imp, expr=c(mean(bmi) ) ) %>% withPool_MI()
             )
# Put in dataframe
df_mean <- data.frame(M = df_mean)
# Set rownames
rownames(df_mean) <- c("orig", "CCA", "RI", "SRI", "MI")
# Calculate difference between original mean and imputed mean
for (i in 1:nrow(df_mean)) {
  df_mean[i, "diff"] <- df_mean["orig", 1] - df_mean[i, 1]
}
## Difference in standard error of BMI
# Calculate SDs
df_sd <- c(sd(df$bmi),
             sd(df_CCA$bmi), 
             sd(df_RI$bmi), 
             sd(df_SRI$bmi),
             with(df_MI_imp, expr=c(sd(bmi) ) ) %>% withPool_MI()
)
# Put in dataframe
df_sd <- data.frame(sd = df_sd)
# Set rownames
rownames(df_sd) <- c("orig", "CCA", "RI", "SRI", "MI")
# Calculate difference between original sd and imputed sd
for (i in 1:nrow(df_sd)) {
  df_sd[i, "diff"] <- df_sd["orig", 1] - df_sd[i, 1]
}

## Difference in correlations between height and weight
# Calculate the correlation for multiple imputation
corr <- micombine.cor(mi.res = df_MI_imp, variables = c("weight", "height")) %>% attr("r_matrix")
corr <- corr[lower.tri(corr)] %>% round(2)
# Calculate other correlations
df_cor <- c(cor(df$weight, df$height),
           cor(df_CCA$weight, df_CCA$height), 
           cor(df_RI$weight, df_RI$height), 
           cor(df_SRI$weight, df_SRI$height),
           corr
)
# Put in dataframe
df_cor <- data.frame(cor = df_cor)
# Set rownames
rownames(df_cor) <- c("orig", "CCA", "RI", "SRI", "MI")
# Calculate difference between original correlation and imputed correlation
for (i in 1:nrow(df_cor)) {
  df_cor[i, "diff"] <- df_cor["orig", 1] - df_cor[i, 1]
}

## Announce winners
c("The most accurate mean is: " = rownames(df_mean[-1, ])[which.min(abs(df_mean[-1, "diff"]))],
  "The most accurate SD is: " = rownames(df_sd[-1, ])[which.min(abs(df_sd[-1, "diff"]))],
  "The most accurate correlation is: " = rownames(df_cor[-1, ])[which.min(abs(df_cor[-1, "diff"]))]
  )

# Remove unnecessary variables
remove(corr, df_mean, df_sd, df_cor, i)
