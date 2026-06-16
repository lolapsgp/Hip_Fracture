# 🧬 Gut Microbiome and Hip Fragility Analysis


This repository contains the full analysis workflow supporting the study:

> **Intestinal *Porphyromonas*, together with serum calcium, predict fragility hip fracture in elderly patients.**

---

## 📌 Project Overview

Osteoporosis is a chronic skeletal disorder characterized by reduced bone mass and increased fracture risk. Traditional clinical predictors such as Bone Mineral Density (BMD) have limited predictive power.

This project investigates the **gut microbiome as a biomarker of fragility hip fractures** in elderly individuals using 16S rRNA sequencing data and advanced statistical modeling.

---

## 🗂️ Repository Structure

├── 1_Phyloseq_object.Rmd

├── 2_Alpha_and_betaDiv.Rmd

├── 2_2_AlphaDiv_gender.R

├── 3_MetadeconfoundR.Rmd

├── 4_Picrust_to_GMMs_and_MetadeconfoundR.Rmd

├── 5_Roc_curves.R

├── Supplementary_plots.R


## 📊 File Description

### 🧱 Data Processing
- **1_Phyloseq_object.Rmd**  
  Creation and preprocessing of the phyloseq object from microbiome sequencing data.

### 📈 Diversity Analysis
- **2_Alpha_and_betaDiv.Rmd**  
  Alpha diversity (within-sample) and beta diversity (between-sample) analyses.

- **2_2_AlphaDiv_gender.R**  
  Alpha diversity analysis stratified by gender.

### ⚖️ Confounder Analysis
- **3_MetadeconfoundR.Rmd**  
  Identification and adjustment of confounders using the *MetadeconfoundR* framework.

### 🧪 Functional Profiling
- **4_Picrust_to_GMMs_and_MetadeconfoundR.Rmd**  
  Functional prediction using PICRUSt and integration with statistical models.

### 📉 Predictive Modeling
- **5_Roc_curves.R**  
  ROC curve analysis and evaluation of predictive models using LOOCV.

### 📎 Supplementary Analyses
- **Supplementary_plots.R**  
  Additional figures and visualizations for supporting results.

---

## 🧪 Methods Overview

- **Microbiome profiling:** 16S rRNA gene sequencing  
- **Statistical analysis:**  
  - Alpha and beta diversity metrics  
  - Differential abundance testing  
  - Generalized Linear Models (GLMs)  
- **Confounder correction:** *MetadeconfoundR*  
- **Functional inference:** PICRUSt  
- **Model validation:** Leave-One-Out Cross-Validation (LOOCV)  
- **Performance evaluation:** ROC curves and AUC  

---



📬 Contact
Lola Giner Pérez
For questions or collaboration, feel free to reach out.
