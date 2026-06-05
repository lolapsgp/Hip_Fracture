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
library(ggpubr)
library(ggplot2)
my_comparisons <- list(c("Control", "Fracture"))
merged_data$Group<- merged_data$Fractura.Caso.Control.Control
# Change "Fractura" to "Fracture" in the Group column
merged_data$Group <- ifelse(merged_data$Group == "Fractura", "Fracture", merged_data$Group)

plot1<-ggplot(merged_data, aes(x = Group, y = Copies_num)) +
  geom_violin(aes(color = Group), trim = FALSE) +  
  geom_boxplot(width = 0.2) +
  theme_classic() + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  stat_compare_means(label= "p.signif", comparisons = my_comparisons)

ggsave("plot1.pdf", plot = plot1, device = "pdf")

plot_sacaled<-ggplot(merged_data, aes(x = Group, y = Copies_num)) +
  geom_violin(aes(color = Group), trim = FALSE) +  
  geom_boxplot(width = 0.2) +
  theme_classic() + 
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90)) +
  stat_compare_means(label= "p.signif", comparisons = my_comparisons) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  labs(y = expression(paste("Copies ", italic("Porphyromonas"), " (scale y log)")))


ggsave("plot_sacaled.pdf", plot = plot_sacaled, device = "pdf")


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


ggplot(merged_data, aes(x = Copies_num, y = Total)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Correlation between Copies_num and Total",
       x = expression(paste("Number of copies qPCR ", italic("Porphyromonas"))),
       y = expression(paste("Total counts ", italic("Porphyromonas")))) +
  theme_minimal()

cor_test <- cor.test(merged_data$Copies_num, merged_data$Total)
print(cor_test)


#-----------variable correlations including qPCR---------
library(readxl)
library(corrplot)
library(GGally)
metadata <- read_excel("~/Documents/Doctorado/Estudios/IATA/Rotura_de_cadera/MetadatosClinMatriz_con_sexo.xlsx")
qPCR_results <- read_excel("~/Documents/Doctorado/Estudios/IATA/Rotura_de_cadera/qPCR/qPCRs_Pgingi_16S_Porphy2for-Porphyrev.xlsx", 
                           sheet = "Hoja2")
metadata<-data.frame(metadata)
qPCR_results<-data.frame(qPCR_results)

qPCR_results <- qPCR_results[order(qPCR_results$FRC), ]
qPCR_results <- qPCR_results[-c(1,2),]

merged_data <- merge(metadata, qPCR_results, by = "FRC")

df <- data.frame(merged_data)
df <- df %>% mutate(Fractura.Caso.Control.Control = ifelse(Fractura.Caso.Control.Control=="Fractura",1,0))
df <- df %>% mutate(Sexo = ifelse(Sexo=="Varon",1,0))
rownames(df)<-df$FRC
df$FRC<-NULL

M = cor(df)
testRes = cor.mtest(df, conf.level = 0.95)
corrplot(M, 
         p.mat = testRes$p, 
         method = "color",
         diag = FALSE, 
         type = 'lower', 
         sig.level = c(0.001, 0.01, 0.05), 
         pch.cex = 0.5, 
         tl.cex = 0.5,
         insig = 'label_sig', 
         pch.col = 'grey20', 
         order = 'hclust', 
         tl.col = 'black',
         tl.srt = 45, 
         col = colorRampPalette(c("dodgerblue3", "white", "red4"))(200))

#More detailed plot with significant correlations
merged_data <- merged_data %>%
  mutate(Group = case_when(
    Fractura.Caso.Control.Control == "Control" ~ "Control",
    Fractura.Caso.Control.Control == "Fractura" ~ "Fracture",
    TRUE ~ NA_character_ # Default case
  ))

newmetadata <- merged_data %>% mutate(Group = factor(Group, 
                                                     levels=c("Control", "Fracture"))) %>% 
  subset(select=c(Copies_num, Edad, PTH, Filtrado_Glomerular, Group)) %>% data.frame()

newmetadata$Age <- newmetadata$Edad
newmetadata$Edad<- NULL
newmetadata$`Glomerular Filtr.`<- newmetadata$Filtrado_Glomerular
newmetadata$Filtrado_Glomerular<- NULL
newmetadata<-newmetadata[,c("Copies_num", "PTH", "Glomerular Filtr.", "Age", "Group")]


newmetadata %>% 
  ggpairs(., columns = 1:5,
          title = "Significant correlations by group", 
          mapping = ggplot2::aes(colour=Group), 
          lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.1), discrete = "points"), 
          diag = list(discrete="barDiag", continuous = wrap("densityDiag", alpha=0.5 )), 
          upper = list(combo = wrap("box_no_facet", alpha=0.5),continuous = wrap("cor", size=4))) + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12, hjust = 1, vjust = 0.999),
        axis.text.y = element_text(size = 12, hjust = 1, vjust = 0.35),
        plot.title.position = "plot",
        plot.title = element_text(size = 12, hjust = 0),
        strip.text = element_text(size = 12, face = "bold")
  )



#Other formats
## leave blank on non-significant coefficient
## add significant correlation coefficients
corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=FALSE, tl.col = 'black')
## add significant level stars
corrplot(M, p.mat = testRes$p, diag = FALSE, type = 'upper',method = 'circle',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE', tl.col = 'black')

######################## Testing agreement #########################
agreement_data <- merged_data[, c("FRC", "Copies_num", "Total", "Group")]
agreement_data$qPCR<-agreement_data$Copies_num
agreement_data$MiSeq<-agreement_data$Total
agreement_data$Total<-NULL
agreement_data$Copies_num<-NULL

# Apply log transformation to both qPCR and MiSeq (add 1 to avoid log(0) errors)
agreement_data <- agreement_data %>%
  mutate(
    log_qPCR = log(qPCR + 1),
    log_MiSeq = log(MiSeq + 1)
  )

# Calculate the mean and difference for the log-transformed data
agreement_data <- agreement_data %>%
  mutate(
    Mean_log = (log_qPCR + log_MiSeq) / 2,
    Difference_log = log_qPCR - log_MiSeq
  )

# Calculate the mean difference and standard deviation of differences for the log-transformed data
mean_diff_log <- mean(agreement_data$Difference_log)
sd_diff_log <- sd(agreement_data$Difference_log)

# Calculate the limits of agreement for the log-transformed data
loa_lower_log <- mean_diff_log - 1.96 * sd_diff_log
loa_upper_log <- mean_diff_log + 1.96 * sd_diff_log

# Print the limits of agreement for the log-transformed data
cat("Lower limit of agreement (log scale):", loa_lower_log, "\n")
cat("Upper limit of agreement (log scale):", loa_upper_log, "\n")

# Create the Bland-Altman plot for the log-transformed data
ggplot(agreement_data, aes(x = Mean_log, y = Difference_log)) +
  geom_point(aes(color = Group), size = 3) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  # Line at 0 difference
  geom_hline(yintercept = loa_lower_log, color = "red", linetype = "dashed") +  # Lower limit of agreement
  geom_hline(yintercept = loa_upper_log, color = "red", linetype = "dashed") +  # Upper limit of agreement
  labs(
    title = "Bland-Altman Plot: qPCR vs MiSeq (Log Scale)",
    x = "Mean of log(qPCR and MiSeq)",
    y = "Difference (log(qPCR) - log(MiSeq))"
  ) +
  theme_minimal()

t.test(agreement_data$Difference_log)


#P-value differences wilcoxon
result <- wilcox.test(agreement_data$qPCR)
p_value <- result$p.value
adjusted_p <- p.adjust(p_value, method = "BH")
adjusted_p
