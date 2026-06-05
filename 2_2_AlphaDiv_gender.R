ps_sexo <- readRDS("~/OneDrive - Universitat de Valencia/Doctorado/Estudios/IATA/Rotura_de_cadera/output/ps_sexo.Rds")
metadata<-data.frame(sample_data(ps_sexo))
metadata$Gender <- metadata$Sexo
metadata <- metadata %>%
  mutate(Gender = case_when(
    Gender == "Mujer" ~ "Female",
    Gender == "Varon" ~ "Male",
    TRUE ~ NA_character_ # Default case
  ))
metadata$Group <- metadata$Fractura.Caso.Control.Control
metadata <- metadata %>%
  mutate(Group = case_when(
    Group == "Control" ~ "Control",
    Group == "Fractura" ~ "Hip Fracture",
    TRUE ~ NA_character_ # Default case
  ))
ps<-ps_sexo
sample_data(ps)<-metadata

newSTorder = c("Female", "Male")
my_comparisons <- list(c("Female", "Male"))
library(ggsignif)
library(ggpubr)


#Shannon
pShannon <- plot_richness(ps, x="Gender", measures=c("Shannon")) +
  geom_violin(aes(color = Gender), trim = FALSE) +  
  geom_boxplot(width = 0.2) +
  theme_classic() + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90), axis.text.x = element_text(size = 12)) +
  stat_compare_means(label= "p.signif", comparisons = my_comparisons, method = "t.test")+
  ylab("Shannon")
#Re-arranging the order of the figure, controls first
pShannon$data$Gender <- as.character(pShannon$data$Gender)
pShannon$data$Gender <- factor(pShannon$data$Gender, levels=newSTorder)
pShannon$layers<-pShannon$layers[-1]
pShannon 


#Chao1
pChao <- plot_richness(ps, x="Gender", measures=c("Chao1")) +
  geom_violin(aes(color = Gender), trim = FALSE) +  
  geom_boxplot(width = 0.2) +
  theme_classic() + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90), axis.text.x = element_text(size = 12)) +
  stat_compare_means(label= "p.signif", comparisons = my_comparisons, method = "t.test")+
  ylab("Chao1")
#Re-arranging the order of the figure, controls first
pChao$data$Gender <- as.character(pChao$data$Gender)
pChao$data$Gender <- factor(pChao$data$Gender, levels=newSTorder)
pChao$layers<-pChao$layers[-1]
pChao 

#ACE
pACE <-  plot_richness(ps, x="Gender", measures=c("ACE")) +
  geom_violin(aes(color = Gender), trim = FALSE) +  
  geom_boxplot(width = 0.2) +
  theme_classic() + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90), axis.text.x = element_text(size = 12)) +
  stat_compare_means(label= "p.signif", comparisons = my_comparisons, method = "t.test")+
  ylab("ACE")
pACE$data$Gender <- as.character(pACE$data$Gender)
pACE$data$Gender <- factor(pACE$data$Gender, levels=newSTorder)
pACE$layers<-pACE$layers[-1]
pACE 
ggarrange(pChao,pACE,pShannon, common.legend = TRUE, legend="bottom",nrow = 1)

#Beta diversity
ps<-ps_sexo
metadata$Group <- metadata$Gender
sample_data(ps)<-metadata
#Beta Diversity
library(tibble)
library(ggplot2)
library(forcats)
#No transformations Bray distance
dist = phyloseq::distance(ps, method="bray")
ordination = ordinate(ps, method="PCoA", distance=dist)
plot_ordination(ps, ordination, color="Group") + 
  theme_classic() +
  theme(strip.background = element_blank()) + stat_ellipse(aes(group = Group), linetype = 2)

#Relative abundance Bray distance
ps_comp<- microbiome::transform(ps, "compositional")
dist = phyloseq::distance(ps_comp, method="bray")
ordination = ordinate(ps_comp, method="PCoA", distance=dist)
plot_ordination(ps_comp, ordination, color="Group") + 
  theme_classic() +
  theme(strip.background = element_blank()) + stat_ellipse(aes(group = Group), linetype = 2)

#Cntered-log-ratio transformation euclidean distance
ps_clr<- microbiome::transform(ps, "clr")

dist = phyloseq::distance(ps_clr, method="euclidean")
ordination = ordinate(ps_clr, method="PCoA", distance=dist)
plot_ordination(ps_clr, ordination, color="Group") + 
  theme_classic() +
  theme(strip.background = element_blank()) + stat_ellipse(aes(group = Group), linetype = 2)

#PERMANOVA
ps<-ps_sexo
pseq.rel <- microbiome::transform(ps, "compositional")
otu <- data.frame(otu_table(pseq.rel))
meta <- data.frame(sample_data(pseq.rel))
permanova <- adonis(t(otu) ~   Group + Sexo, 
                    data = meta, permutations=999, method = "bray", by = "margin")
library(repmod)
permanova_table<-data.frame(permanova$aov.tab)
make_csv_table(permanova_table, "~/PERMANOVA_table", info = "")


