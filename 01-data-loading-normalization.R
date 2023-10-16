#Load libraries
library(tidyverse)
library(ggfortify)
library(ggridges)
library(limma)
library(vsn)

#Load dataset
dat <- read.delim("proteinGroups.txt")

#Rename samples
colnames <- colnames(dat)
label <- c("HD013","P022nm","HD014","P016nm","HD021","P018nm","HD031","P017","HD004","P015", "P006","P016","P023","P018","P019","P022")
collabs <- colnames(dat)
collabs[43:58] <- label
colnames(dat) <- collabs


#filter reverse, contaminate and unique peptides
dat <- dat %>% filter(Only.identified.by.site == "")
dat <- dat %>% filter(Reverse == "")
dat <- dat %>% filter(Unique.peptides > 0)
dat <- dat %>% filter(Potential.contaminant == "")



#log2 appropiate cols
dat1 <- dat
dat1[,43:58] <- lapply(dat1[43:58], log2)


#remove -Inf
dat2 <- dat1[apply(dat1[,c(43:58)],1, function(x) !any(x=="-Inf")),]

#VSN norm
dat_vsn <- dat2
dat_vsn[,43:58] <- as.data.frame(normalizeVSN(dat_vsn[,43:58]))


#pivot longer
dat3 <- pivot_longer(data = dat2, cols = 43:58, names_to = "sample", values_to = "intensity")
#relevel
dat3 <- dat3 %>%
  mutate(sample = fct_relevel(sample, levels = "P022nm","P018nm","P016nm","HD013","HD014","HD021","HD031","P017","HD004","P015", "P006","P016","P023","P018","P019","P022"))

dat_vsn_long <- pivot_longer(data = dat_vsn, cols = 43:58, names_to = "sample", values_to = "intensity")
#relevel
dat_vsn_long <- dat_vsn_long %>%
  mutate(sample = fct_relevel(sample, levels = "P022nm","P018nm","P016nm","HD013","HD014","HD021","HD031","P017","HD004","P015", "P006","P016","P023","P018","P019","P022"))

#plot profile plots before and after normalization
profilplots <- ggplot(dat3, aes(x = sample, y = intensity, group = Gene.names))
profilplots+geom_line()+geom_boxplot(aes(group = sample))

profilplots_norm <- ggplot(dat_vsn_long, aes(x = sample, y = intensity, group = Gene.names))
profilplots_norm+geom_line()+geom_boxplot(aes(group = sample))

#plot density plots
p <- ggplot(dat3, aes(x = intensity, color = sample))
P_vsn <- ggplot(dat_vsn_long, aes(x = intensity, color = sample))
p+geom_density()
P_vsn+geom_density()


#plot scale ridgeline density plots
ridge <- ggplot(dat3, aes(x = intensity, y = sample, fill = stat(x)))
ridge + geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "intensity", option = "C")

ridge_norm <- ggplot(dat_vsn_long, aes(x = intensity, y = sample, fill = stat(x)))
ridge_norm + geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "intensity", option = "C")

#plot PCA analysis before and after normalization
transpose <- dat2[43:58]
rownames(transpose) <- make.unique(dat2$Gene.names)
transpos <- t(transpose)
pca_res <- prcomp(transpos, scale. = TRUE)
autoplot(pca_res)

transpose_norm <- dat_vsn[43:58]
rownames(transpose_norm) <- make.unique(dat2$Gene.names)
transpos_norm <- t(transpose_norm)
pca_res_norm <- prcomp(transpos_norm, scale. = TRUE)
autoplot(pca_res_norm, label = T, label.size = 3)
autoplot(pca_res_norm)+geom_label(aes(label = rownames(transpos_norm)), nudge_x = .03, size = 1.5)
