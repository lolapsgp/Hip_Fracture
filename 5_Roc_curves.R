#qPCR results
library(readxl)
metadata <- read_excel("~/Documents/Doctorado/Estudios/IATA/Rotura_de_cadera/MetadatosClinMatriz_con_sexo.xlsx")
qPCR_results <- read_excel("~/Documents/Doctorado/Estudios/IATA/Rotura_de_cadera/qPCR/qPCRs_Pgingi_16S_Porphy2for-Porphyrev.xlsx", 
                           sheet = "Hoja2")
metadata<-data.frame(metadata)
qPCR_results<-data.frame(qPCR_results)

qPCR_results <- qPCR_results[order(qPCR_results$FRC), ]
qPCR_results <- qPCR_results[-c(1,2),]

merged_data <- merge(metadata, qPCR_results, by = "FRC")
#Correlations
Abundances <- read_excel("~/Documents/Doctorado/Estudios/IATA/Rotura_de_cadera/qPCR/qPCRs_Pgingi_16S_Porphy2for-Porphyrev.xlsx", 
                         sheet = "Abundancias_porphyromonas_seq")
Abundances<-(data.frame(Abundances))
Abundances <- Abundances[rowSums(is.na(Abundances)) != ncol(Abundances), ]

rownames(Abundances)<-Abundances$`...1`
Abundances$`...1`<-NULL
Abundances<-data.frame(t(Abundances))
Abundances$FRC<-rownames(Abundances)
Abundances<-Abundances[,c("Total", "FRC")]
Abundances <- Abundances[order(Abundances$FRC), ]
merged_data <- merge(merged_data, Abundances, by = "FRC")

merged_data$Group <-merged_data$Fractura.Caso.Control.Control
merged_data$Group <- ifelse(merged_data$Group == "Fractura", 1, 0)

######################## Copies qPCR + Calcio_total ########################
merged_data2<-data.frame(merged_data[,c("Group", "Copies_num", "Calcio_total")])
rownames(merged_data2)<-merged_data$FRC

# Logistic regression
model <- glm(Group ~ Copies_num + Calcio_total, data = merged_data, family = binomial)

# Predict probabilities
merged_data2$prob <- predict(model, type = "response")

# ROC curve
library(performance)
roc_obj <- performance_roc(model)
plot(roc_obj)

auc <- sum(diff(roc_obj$Specificity) * head(roc_obj$Sensitivity, -1))
print(auc)

# Get optimal threshold and convert (best specificity and sensitivity)
library(performance)
library(dplyr)

# Create thresholds manually
thresholds <- sort(unique(merged_data2$prob))

# Initialize empty list to store ROC values
roc_data <- data.frame(
  Threshold = thresholds,
  Sensitivity = numeric(length(thresholds)),
  Specificity = numeric(length(thresholds))
)

for (i in seq_along(thresholds)) {
  t <- thresholds[i]
  predicted <- ifelse(merged_data2$prob >= t, 1, 0)
  
  TP <- sum(predicted == 1 & merged_data$Group == 1)
  TN <- sum(predicted == 0 & merged_data$Group == 0)
  FP <- sum(predicted == 1 & merged_data$Group == 0)
  FN <- sum(predicted == 0 & merged_data$Group == 1)
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  
  roc_data$Sensitivity[i] <- sensitivity
  roc_data$Specificity[i] <- specificity
}

# Compute Youden Index
roc_data$YoudenIndex <- roc_data$Sensitivity + roc_data$Specificity - 1

# Find best threshold
best <- roc_data[which.max(roc_data$YoudenIndex), ]
print(best)

set.seed(123)
n_boot <- 1000
auc_values <- numeric(n_boot)

for (i in 1:n_boot) {
  boot_idx <- sample(nrow(merged_data), replace = TRUE)
  boot_model <- glm(Group ~ Total + Proteína_C_reactiva, data = merged_data[boot_idx, ], family = binomial)
  boot_roc <- performance::performance_roc(boot_model)
  auc_values[i] <- sum(diff(boot_roc$Specificity) * head(boot_roc$Sensitivity, -1))
}

# 95% CI
ci <- quantile(auc_values, probs = c(0.025, 0.975))
cat("95% CI for AUC:", round(ci[1], 3), "-", round(ci[2], 3), "\n")

#VIF (Variance Inflation Factor)
library(car)
vif(model)

########################LOOCV model copies + Calcio_Total #####################
# Prepare data
df <- merged_data[, c("Group", "Copies_num", "Calcio_total")]
df <- na.omit(df)

n <- nrow(df)
pred_probs <- numeric(n)

# LOOCV loop
for (i in 1:n) {
  train_data <- df[-i, ]
  test_data  <- df[i, , drop = FALSE]
  
  # Fit model on training set
  model_loocv <- glm(Group ~ Copies_num + Calcio_total,
                     data = train_data,
                     family = binomial)
  
  # Predict on left-out sample
  pred_probs[i] <- predict(model_loocv, newdata = test_data, type = "response")
}

# Store predictions
df$pred_loocv <- pred_probs

library(performance)

# Create a fake model using LOOCV predictions
model_loocv_fake <- glm(Group ~ pred_loocv, data = df, family = binomial)

roc_obj_loocv <- performance_roc(model_loocv_fake)
plot(roc_obj_loocv)

auc_loocv <- sum(diff(roc_obj_loocv$Specificity) *
                   head(roc_obj_loocv$Sensitivity, -1))

print(auc_loocv)

######################## Abundances 16S + Calcio_total ########################
merged_data2<-data.frame(merged_data[,c("Group", "Total", "Calcio_total")])
rownames(merged_data2)<-merged_data$FRC

# Logistic regression
model <- glm(Group ~ Total + Calcio_total, data = merged_data, family = binomial)

# Predict probabilities
merged_data2$prob <- predict(model, type = "response")

# ROC curve
library(performance)
roc_obj <- performance_roc(model)
plot(roc_obj)

auc <- sum(diff(roc_obj$Specificity) * head(roc_obj$Sensitivity, -1))
print(auc)

# Get optimal threshold and convert (best specificity and sensitivity)
library(performance)
library(dplyr)

# Create thresholds manually
thresholds <- sort(unique(merged_data2$prob))

# Initialize empty list to store ROC values
roc_data <- data.frame(
  Threshold = thresholds,
  Sensitivity = numeric(length(thresholds)),
  Specificity = numeric(length(thresholds))
)

for (i in seq_along(thresholds)) {
  t <- thresholds[i]
  predicted <- ifelse(merged_data2$prob >= t, 1, 0)
  
  TP <- sum(predicted == 1 & merged_data$Group == 1)
  TN <- sum(predicted == 0 & merged_data$Group == 0)
  FP <- sum(predicted == 1 & merged_data$Group == 0)
  FN <- sum(predicted == 0 & merged_data$Group == 1)
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  
  roc_data$Sensitivity[i] <- sensitivity
  roc_data$Specificity[i] <- specificity
}

# Compute Youden Index
roc_data$YoudenIndex <- roc_data$Sensitivity + roc_data$Specificity - 1

# Find best threshold
best <- roc_data[which.max(roc_data$YoudenIndex), ]
print(best)

set.seed(123)
n_boot <- 1000
auc_values <- numeric(n_boot)

for (i in 1:n_boot) {
  boot_idx <- sample(nrow(merged_data), replace = TRUE)
  boot_model <- glm(Group ~ Total + Proteína_C_reactiva, data = merged_data[boot_idx, ], family = binomial)
  boot_roc <- performance::performance_roc(boot_model)
  auc_values[i] <- sum(diff(boot_roc$Specificity) * head(boot_roc$Sensitivity, -1))
}

# 95% CI
ci <- quantile(auc_values, probs = c(0.025, 0.975))
cat("95% CI for AUC:", round(ci[1], 3), "-", round(ci[2], 3), "\n")

#VIF (Variance Inflation Factor)
library(car)
vif(model)


######################## Calcio_total ########################
merged_data2<-data.frame(merged_data[,c("Group", "Calcio_total")])
rownames(merged_data2)<-merged_data$FRC

# Logistic regression
model <- glm(Group ~ Calcio_total, data = merged_data, family = binomial)

# Predict probabilities
merged_data2$prob <- predict(model, type = "response")

# ROC curve
library(performance)
roc_obj <- performance_roc(model)
plot(roc_obj)

auc <- sum(diff(roc_obj$Specificity) * head(roc_obj$Sensitivity, -1))
print(auc)

# Get optimal threshold and convert (best specificity and sensitivity)
library(performance)
library(dplyr)

# Create thresholds manually
thresholds <- sort(unique(merged_data2$prob))

# Initialize empty list to store ROC values
roc_data <- data.frame(
  Threshold = thresholds,
  Sensitivity = numeric(length(thresholds)),
  Specificity = numeric(length(thresholds))
)

for (i in seq_along(thresholds)) {
  t <- thresholds[i]
  predicted <- ifelse(merged_data2$prob >= t, 1, 0)
  
  TP <- sum(predicted == 1 & merged_data$Group == 1)
  TN <- sum(predicted == 0 & merged_data$Group == 0)
  FP <- sum(predicted == 1 & merged_data$Group == 0)
  FN <- sum(predicted == 0 & merged_data$Group == 1)
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  
  roc_data$Sensitivity[i] <- sensitivity
  roc_data$Specificity[i] <- specificity
}

# Compute Youden Index
roc_data$YoudenIndex <- roc_data$Sensitivity + roc_data$Specificity - 1

# Find best threshold
best <- roc_data[which.max(roc_data$YoudenIndex), ]
print(best)

set.seed(123)
n_boot <- 1000
auc_values <- numeric(n_boot)

for (i in 1:n_boot) {
  boot_idx <- sample(nrow(merged_data), replace = TRUE)
  boot_model <- glm(Group ~ Total + Proteína_C_reactiva, data = merged_data[boot_idx, ], family = binomial)
  boot_roc <- performance::performance_roc(boot_model)
  auc_values[i] <- sum(diff(boot_roc$Specificity) * head(boot_roc$Sensitivity, -1))
}

# 95% CI
ci <- quantile(auc_values, probs = c(0.025, 0.975))
cat("95% CI for AUC:", round(ci[1], 3), "-", round(ci[2], 3), "\n")

#VIF (Variance Inflation Factor)
library(car)
vif(model)


######################## Copies qPCR  ########################
merged_data2<-data.frame(merged_data[,c("Group", "Copies_num", "Calcio_total")])
rownames(merged_data2)<-merged_data$FRC

# Logistic regression
model <- glm(Group ~ Copies_num , data = merged_data, family = binomial)

# Predict probabilities
merged_data2$prob <- predict(model, type = "response")

# ROC curve
library(performance)
roc_obj <- performance_roc(model)
plot(roc_obj)

auc <- sum(diff(roc_obj$Specificity) * head(roc_obj$Sensitivity, -1))
print(auc)

# Get optimal threshold and convert (best specificity and sensitivity)
library(performance)
library(dplyr)

# Create thresholds manually
thresholds <- sort(unique(merged_data2$prob))

# Initialize empty list to store ROC values
roc_data <- data.frame(
  Threshold = thresholds,
  Sensitivity = numeric(length(thresholds)),
  Specificity = numeric(length(thresholds))
)

for (i in seq_along(thresholds)) {
  t <- thresholds[i]
  predicted <- ifelse(merged_data2$prob >= t, 1, 0)
  
  TP <- sum(predicted == 1 & merged_data$Group == 1)
  TN <- sum(predicted == 0 & merged_data$Group == 0)
  FP <- sum(predicted == 1 & merged_data$Group == 0)
  FN <- sum(predicted == 0 & merged_data$Group == 1)
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  
  roc_data$Sensitivity[i] <- sensitivity
  roc_data$Specificity[i] <- specificity
}

# Compute Youden Index
roc_data$YoudenIndex <- roc_data$Sensitivity + roc_data$Specificity - 1

# Find best threshold
best <- roc_data[which.max(roc_data$YoudenIndex), ]
print(best)

set.seed(123)
n_boot <- 1000
auc_values <- numeric(n_boot)

for (i in 1:n_boot) {
  boot_idx <- sample(nrow(merged_data), replace = TRUE)
  boot_model <- glm(Group ~ Total + Proteína_C_reactiva, data = merged_data[boot_idx, ], family = binomial)
  boot_roc <- performance::performance_roc(boot_model)
  auc_values[i] <- sum(diff(boot_roc$Specificity) * head(boot_roc$Sensitivity, -1))
}

# 95% CI
ci <- quantile(auc_values, probs = c(0.025, 0.975))
cat("95% CI for AUC:", round(ci[1], 3), "-", round(ci[2], 3), "\n")



######################## Abundances 16S  ########################
merged_data2<-data.frame(merged_data[,c("Group", "Total", "Calcio_total")])
rownames(merged_data2)<-merged_data$FRC

# Logistic regression
model <- glm(Group ~ Total , data = merged_data, family = binomial)

# Predict probabilities
merged_data2$prob <- predict(model, type = "response")

# ROC curve
library(performance)
roc_obj <- performance_roc(model)
plot(roc_obj)

auc <- sum(diff(roc_obj$Specificity) * head(roc_obj$Sensitivity, -1))
print(auc)

# Get optimal threshold and convert (best specificity and sensitivity)
library(performance)
library(dplyr)

# Create thresholds manually
thresholds <- sort(unique(merged_data2$prob))

# Initialize empty list to store ROC values
roc_data <- data.frame(
  Threshold = thresholds,
  Sensitivity = numeric(length(thresholds)),
  Specificity = numeric(length(thresholds))
)

for (i in seq_along(thresholds)) {
  t <- thresholds[i]
  predicted <- ifelse(merged_data2$prob >= t, 1, 0)
  
  TP <- sum(predicted == 1 & merged_data$Group == 1)
  TN <- sum(predicted == 0 & merged_data$Group == 0)
  FP <- sum(predicted == 1 & merged_data$Group == 0)
  FN <- sum(predicted == 0 & merged_data$Group == 1)
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  
  roc_data$Sensitivity[i] <- sensitivity
  roc_data$Specificity[i] <- specificity
}

# Compute Youden Index
roc_data$YoudenIndex <- roc_data$Sensitivity + roc_data$Specificity - 1

# Find best threshold
best <- roc_data[which.max(roc_data$YoudenIndex), ]
print(best)

set.seed(123)
n_boot <- 1000
auc_values <- numeric(n_boot)

for (i in 1:n_boot) {
  boot_idx <- sample(nrow(merged_data), replace = TRUE)
  boot_model <- glm(Group ~ Total + Proteína_C_reactiva, data = merged_data[boot_idx, ], family = binomial)
  boot_roc <- performance::performance_roc(boot_model)
  auc_values[i] <- sum(diff(boot_roc$Specificity) * head(boot_roc$Sensitivity, -1))
}

# 95% CI
ci <- quantile(auc_values, probs = c(0.025, 0.975))
cat("95% CI for AUC:", round(ci[1], 3), "-", round(ci[2], 3), "\n")


