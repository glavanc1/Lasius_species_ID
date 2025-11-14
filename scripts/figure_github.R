# Lasius Paper Figures
# 2025
# R project contains : /scripts /data /figures
# PACKAGES####
library(vcfR)
library(ape)
library(tidyverse)
library(MASS)
library(adegenet)
library(ggpubr)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
library(reshape2)
library(plotly)
library(ggtree)
library(TDbook)
library(aplot)


# Common variables ####
vcf.contam <- read.vcfR("data/conta_Lasius_1snp_DP8_meanDP200_mac2_miss75_miss_ind_0.5.vcf.gz")

compmap.vcf <- read.vcfR("data/Lasius_1snp_DP8_meanDP200_mac2_miss75_miss_ind_0.5.vcf.gz")

ADR.vcf0525 <- read.vcfR("data/Lasius_1snp_DP8_meanDP200_mac2_miss75_ADR_filtered_discard0.05_correct0.25.vcf.gz")

haploid.vcf <- read.vcfR("data/haploidLasius_1snp_DP8_meanDP200_mac2_miss75_miss_ind_0.5.vcf.gz")


vcf.mds <- function(vcf) {
  gl <- vcfR2genlight(vcf)
  d <- dist(gl)
  m <- as.matrix(d)
  rm(Lasi_)
  rm(`Lasi Chtono`)
  rm(`Lasi jaune`)
  rm(Lasi.Chtono)
  rm(Lasi.jaune)
  rm(Lasi_alie)
  rm(`Lasi_alie gr`)
  rm(Lasi_alie.gr)
  rm(Lasi_brun)
  rm(Lasi_emar)
  rm(Lasi_flav)
  rm(Lasi_fuli)
  rm(Lasi_mixt)
  rm(Lasi_myop)
  rm(Lasi_nige)
  rm(Lasi_umbr)
  rm(Lasi_nige.plat)
  rm(`Lasi_nige/plat`)
  rm(Lasi_para)
  rm(Lasi_plat)
  rm(Lasi_psam)
  ants <- read.csv("data/antmetadata.csv")
  colnames(ants)[1] <- "CATALOGUENUMBER"

  id.list <- read.table("data/all_individuals.txt")
  lasius <- subset(ants, grepl("Las", SPECIESSUMMARY))

  id.las <- lasius[lasius[, "CATALOGUENUMBER"] %in% unlist(id.list$V1), ]
  id.las <- as.data.frame(cbind(id.las$CATALOGUENUMBER, id.las$SPECIESSUMMARY))

  ir <- split(id.las, id.las$V2)

  for (name in (names(ir))) {
    list <- (ir[[name]]$V1)
    assign(name, list)
  }
  Lasi_plat <- c(Lasi_plat, "4457") # mislabels
  Lasi_flav <- Lasi_flav[Lasi_flav != "9990946"]

  MDS_all <- isoMDS(d = m, k = 2, maxit = 100)

  MDS.p <- as.data.frame(MDS_all$points)
  MDS.col <- MDS.p %>%
    mutate(oldcolors = case_when(
      try((rownames(MDS.p) %in% Lasi) ~ "grey40"),
      try((rownames(MDS.p) %in% `Lasi Chtono`) ~ "grey40"),
      try((rownames(MDS.p) %in% `Lasi jaune`) ~ "grey40"),
      try((rownames(MDS.p) %in% Lasi_alie) ~ "#27AAE1"),
      try((rownames(MDS.p) %in% `Lasi_alie gr`) ~ "grey40"),
      try((rownames(MDS.p) %in% Lasi_brun) ~ "#006838"),
      try((rownames(MDS.p) %in% Lasi_emar) ~ "#8DC63F"),
      try((rownames(MDS.p) %in% Lasi_flav) ~ "#F9ED32"),
      try((rownames(MDS.p) %in% Lasi_fuli) ~ "#1B75BC"),
      try((rownames(MDS.p) %in% Lasi_mixt) ~ "#2B3990"),
      try((rownames(MDS.p) %in% Lasi_myop) ~ "#92278F"),
      try((rownames(MDS.p) %in% Lasi_nige) ~ "#EC008C"),
      try((rownames(MDS.p) %in% `Lasi_nige/plat`) ~ "grey40"),
      try((rownames(MDS.p) %in% Lasi_para) ~ "#F15A29"),
      try((rownames(MDS.p) %in% Lasi_plat) ~ "#A97C50"),
      try((rownames(MDS.p) %in% Lasi_psam) ~ "#00A79D"),
      try((rownames(MDS.p) %in% Lasi_umbr) ~ "#C2B59B")
    ))

  MDS.col <- MDS.col %>%
    mutate(Species = case_when(
      try((rownames(MDS.p) %in% Lasi) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% `Lasi Chtono`) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% `Lasi jaune`) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% Lasi_alie) ~ "L.alienus"),
      try((rownames(MDS.p) %in% `Lasi_alie gr`) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% Lasi_brun) ~ "L.brunneus"),
      try((rownames(MDS.p) %in% Lasi_emar) ~ "L.emarginatus"),
      try((rownames(MDS.p) %in% Lasi_flav) ~ "L.flavus"),
      try((rownames(MDS.p) %in% Lasi_fuli) ~ "L.fuliginosus"),
      try((rownames(MDS.p) %in% Lasi_mixt) ~ "L.mixtus"),
      try((rownames(MDS.p) %in% Lasi_myop) ~ "L.myops"),
      try((rownames(MDS.p) %in% Lasi_nige) ~ "L.niger"),
      try((rownames(MDS.p) %in% `Lasi_nige/plat`) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% Lasi_para) ~ "L.paralienus"),
      try((rownames(MDS.p) %in% Lasi_plat) ~ "L.platythorax"),
      try((rownames(MDS.p) %in% Lasi_psam) ~ "L.neglectus"),
      try((rownames(MDS.p) %in% Lasi_umbr) ~ "L.umbratus")
    ))


  # corrected ids####
  rm(Lasi_)
  rm(`Lasi Chtono`)
  rm(`Lasi jaune`)
  rm(Lasi.Chtono)
  rm(Lasi.jaune)
  rm(Lasi_alie)
  rm(`Lasi_alie gr`)
  rm(Lasi_alie.gr)
  rm(Lasi_brun)
  rm(Lasi_emar)
  rm(Lasi_flav)
  rm(Lasi_fuli)
  rm(Lasi_mixt)
  rm(Lasi_myop)
  rm(Lasi_nige)
  rm(Lasi_umbr)
  rm(Lasi_nige.plat)
  rm(`Lasi_nige/plat`)
  rm(Lasi_para)
  rm(Lasi_plat)
  rm(Lasi_psam)

  correctedID <- read.csv("data/correctID1.csv", header = T)
  for (name in (colnames(correctedID))) {
    list <- na.exclude(correctedID[[name]])
    assign(name, list)
  }

  MDS.col <- MDS.col %>%
    mutate(New.Species = case_when(
      try((rownames(MDS.p) %in% Lasi) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% Lasi_alie) ~ "L.alienus"),
      try((rownames(MDS.p) %in% Lasi_brun) ~ "L.brunneus"),
      try((rownames(MDS.p) %in% Lasi_emar) ~ "L.emarginatus"),
      try((rownames(MDS.p) %in% Lasi_flav) ~ "L.flavus"),
      try((rownames(MDS.p) %in% Lasi_fuli) ~ "L.fuliginosus"),
      try((rownames(MDS.p) %in% Lasi_mixt) ~ "L.mixtus"),
      try((rownames(MDS.p) %in% Lasi_myop) ~ "L.myops"),
      try((rownames(MDS.p) %in% Lasi_nige) ~ "L.niger"),
      try((rownames(MDS.p) %in% Lasi_para) ~ "L.paralienus"),
      try((rownames(MDS.p) %in% Lasi_plat) ~ "L.platythorax"),
      try((rownames(MDS.p) %in% Lasi_psam) ~ "L.neglectus"),
      try((rownames(MDS.p) %in% Lasi_umbr) ~ "L.umbratus")
    ))
  MDS.col <- MDS.col %>%
    mutate(colors = case_when(
      try((rownames(MDS.p) %in% Lasi_alie) ~ "#27AAE1"),
      try((rownames(MDS.p) %in% Lasi_brun) ~ "#006838"),
      try((rownames(MDS.p) %in% Lasi_emar) ~ "#8DC63F"),
      try((rownames(MDS.p) %in% Lasi_flav) ~ "#F9ED32"),
      try((rownames(MDS.p) %in% Lasi_fuli) ~ "#1B75BC"),
      try((rownames(MDS.p) %in% Lasi_mixt) ~ "#2B3990"),
      try((rownames(MDS.p) %in% Lasi_myop) ~ "#92278F"),
      try((rownames(MDS.p) %in% Lasi_nige) ~ "#EC008C"),
      try((rownames(MDS.p) %in% Lasi_para) ~ "#F15A29"),
      try((rownames(MDS.p) %in% Lasi_plat) ~ "#A97C50"),
      try((rownames(MDS.p) %in% Lasi_psam) ~ "#00A79D"),
      try((rownames(MDS.p) %in% Lasi_umbr) ~ "#C2B59B")
    ))

  MDS.col <- MDS.col %>%
    mutate(Spec.change = ifelse(MDS.col$Species == MDS.col$New.Species, "Same as Morphological Identification", "Different from Morphological Identification"))

  MDS.col$idnum <- rownames(MDS.col)

  MDS.col <- MDS.col %>%
    mutate(alph = case_when(
      try((Spec.change == "Same as Morphological Identification") ~ 0.2),
      try((Spec.change == "Different from Morphological Identification") ~ 1)
    ))

  return(MDS.col)
}

MDS.col.contam <- vcf.mds(vcf.contam)
MDS.col.adr <- vcf.mds(ADR.vcf0525)
MDS.col.comp <- vcf.mds(compmap.vcf)
MDS.col.hap <- vcf.mds(haploid.vcf)
# FIGURE 2. ADR and Het missing comparison ####

vcf.ad.plot <- function(vcf) {
  dp <- extract.gt(vcf, "DP", as.numeric = T)

  ad.ref <- extract.gt(vcf, "AD", as.numeric = T)

  gt <- extract.gt(vcf, "GT")
  ad.ref.het <- ad.ref
  ad.ref.het[which(gt %in% c("0/0", "1/1"))] <- NA

  ad.alt.het <- dp - ad.ref
  ad.alt.het[which(gt %in% c("0/0", "1/1"))] <- NA

  get.ad.ratio <- function(ad.ref, ad.alt) {
    return(pmax(ad.ref, ad.alt) / psum(ad.ref, ad.alt))
  } # psum comes from the package 'rccmisc'
  psum <- function(..., na.rm = FALSE) {
    rowSums(do.call(cbind, list(...)), na.rm = na.rm)
  }

  all.lines <- data.frame()

  for (id in colnames(ad.ref)) {
    ad.ratio <- mapply(get.ad.ratio, ad.ref.het[, which(colnames(ad.ref) == id)], ad.alt.het[, which(colnames(ad.ref) == id)])
    if (length(ad.ratio[which(ad.ratio != 1)]) > 1) {
      dens <- density(ad.ratio[which(ad.ratio != 1)], na.rm = T)
      dens.data <- as.data.frame(cbind(dens[["x"]], dens[["y"]]))
      all.lines <- rbind(all.lines, cbind(id, dens.data))
    }
  }
  ad.ggplot <- ggplot() +
    scale_x_continuous(name = "Allelic depth ratio (heterozygous alleles only)", limits = c(0.45, 1)) +
    scale_y_continuous(name = "Frequency", limits = c(0, 45), breaks = NULL) +
    geom_line(data = all.lines, mapping = aes(x = V1, y = V2, group = id), col = "#000000") +
    theme_bw(base_size = 20) +
    theme(axis.ticks = element_blank())

  return(ad.ggplot)
}

ggcontam <- vcf.ad.plot(vcf.contam) + ggtitle("Original data")
ggADR0525 <- vcf.ad.plot(ADR.vcf0525) + ggtitle("Competitive mapping + ADR")
ggredo <- vcf.ad.plot(compmap.vcf) + ggtitle("Competitive mapping")

vcf.het.plot <- function(vcf) {
  dp <- extract.gt(vcf, "DP", as.numeric = T)

  ad.ref <- extract.gt(vcf, "AD", as.numeric = T)

  gt <- extract.gt(vcf, "GT")
  gt.d <- as.data.frame(gt)

  missing <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(missing) <- c("sampleid", "prop.missing", "prop.heterozygous")
  for (i in colnames(gt.d)) {
    sample <- gt.d[, i]
    het.tab <- as.data.frame(table(sample, useNA = "no"))
    rownames(het.tab) <- het.tab[, 1]
    prop.miss <- (sum(is.na(sample))) / length(sample)
    prop.het <- ((het.tab["0/1", "Freq"]) / sum(het.tab$Freq))

    row <- cbind(i, as.numeric(prop.miss), as.numeric(prop.het))
    missing <- rbind(missing, row)
  }
  colnames(missing) <- c("id", "prop.miss", "prop.het")
  missing[is.na(missing)] <- 0


  hetplot <- ggplot(missing, aes(x = as.numeric(prop.miss), y = as.numeric(prop.het))) +
    geom_point(size = 2, col = "#000000") +
    labs(x = "Proportion of missing genotypes", y = "Proportion of heterozygous loci") +
    ylim(0, 0.17) +
    xlim(0, 0.6) +
    theme_bw(base_size = 20)
  return(hetplot)
}

hetcon <- vcf.het.plot(vcf.contam) + ggtitle("Original data")
hetADR0525 <- vcf.het.plot(ADR.vcf0525) + ggtitle("Competitive mapping + ADR")
hetredo <- vcf.het.plot(compmap.vcf) + ggtitle("Competitive mapping")

setEPS()
postscript("./figures/ADR_Het_compare.eps", width = 25, height = 18)

ggarrange(ggcontam, ggredo, ggADR0525,
  hetcon, hetredo, hetADR0525,
  labels = c("A", "C", "E", "B", "D", "F"), font.label = list(size = 20),
  ncol = 3, nrow = 2
)
dev.off()

# FIGURE 3. MDS comparison ####
vcf.mds.plot <- function(vcf) {
  gl <- vcfR2genlight(vcf)
  d <- dist(gl)
  m <- as.matrix(d)
  rm(Lasi_)
  rm(`Lasi Chtono`)
  rm(`Lasi jaune`)
  rm(Lasi.Chtono)
  rm(Lasi.jaune)
  rm(Lasi_alie)
  rm(`Lasi_alie gr`)
  rm(Lasi_alie.gr)
  rm(Lasi_brun)
  rm(Lasi_emar)
  rm(Lasi_flav)
  rm(Lasi_fuli)
  rm(Lasi_mixt)
  rm(Lasi_myop)
  rm(Lasi_nige)
  rm(Lasi_umbr)
  rm(Lasi_nige.plat)
  rm(`Lasi_nige/plat`)
  rm(Lasi_para)
  rm(Lasi_plat)
  rm(Lasi_psam)
  ants <- read.csv("data/antmetadata.csv")
  colnames(ants)[1] <- "CATALOGUENUMBER"

  id.list <- read.table("data/all_individuals.txt")
  lasius <- subset(ants, grepl("Las", SPECIESSUMMARY))

  id.las <- lasius[lasius[, "CATALOGUENUMBER"] %in% unlist(id.list$V1), ]
  id.las <- as.data.frame(cbind(id.las$CATALOGUENUMBER, id.las$SPECIESSUMMARY))

  ir <- split(id.las, id.las$V2)

  for (name in (names(ir))) {
    list <- (ir[[name]]$V1)
    assign(name, list)
  }
  Lasi_plat <- c(Lasi_plat, "4457")
  Lasi_flav <- Lasi_flav[Lasi_flav != "9990946"]

  MDS_all <- isoMDS(d = d, k = 2, maxit = 100)

  MDS.p <- as.data.frame(MDS_all$points)
  MDS.col <- MDS.p %>%
    mutate(colors = case_when(
      try((rownames(MDS.p) %in% Lasi) ~ "grey40"),
      try((rownames(MDS.p) %in% `Lasi Chtono`) ~ "grey40"),
      try((rownames(MDS.p) %in% `Lasi jaune`) ~ "grey40"),
      try((rownames(MDS.p) %in% Lasi_alie) ~ "#27AAE1"),
      try((rownames(MDS.p) %in% `Lasi_alie gr`) ~ "grey40"),
      try((rownames(MDS.p) %in% Lasi_brun) ~ "#006838"),
      try((rownames(MDS.p) %in% Lasi_emar) ~ "#8DC63F"),
      try((rownames(MDS.p) %in% Lasi_flav) ~ "#F9ED32"),
      try((rownames(MDS.p) %in% Lasi_fuli) ~ "#1B75BC"),
      try((rownames(MDS.p) %in% Lasi_mixt) ~ "#2B3990"),
      try((rownames(MDS.p) %in% Lasi_myop) ~ "#92278F"),
      try((rownames(MDS.p) %in% Lasi_nige) ~ "#EC008C"),
      try((rownames(MDS.p) %in% `Lasi_nige/plat`) ~ "grey40"),
      try((rownames(MDS.p) %in% Lasi_para) ~ "#F15A29"),
      try((rownames(MDS.p) %in% Lasi_plat) ~ "#A97C50"),
      try((rownames(MDS.p) %in% Lasi_psam) ~ "#00A79D"),
      try((rownames(MDS.p) %in% Lasi_umbr) ~ "#C2B59B")
    ))

  MDS.col <- MDS.col %>%
    mutate(Species = case_when(
      try((rownames(MDS.p) %in% Lasi) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% `Lasi Chtono`) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% `Lasi jaune`) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% Lasi_alie) ~ "L.alienus"),
      try((rownames(MDS.p) %in% `Lasi_alie gr`) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% Lasi_brun) ~ "L.brunneus"),
      try((rownames(MDS.p) %in% Lasi_emar) ~ "L.emarginatus"),
      try((rownames(MDS.p) %in% Lasi_flav) ~ "L.flavus"),
      try((rownames(MDS.p) %in% Lasi_fuli) ~ "L.fuliginosus"),
      try((rownames(MDS.p) %in% Lasi_mixt) ~ "L.mixtus"),
      try((rownames(MDS.p) %in% Lasi_myop) ~ "L.myops"),
      try((rownames(MDS.p) %in% Lasi_nige) ~ "L.niger"),
      try((rownames(MDS.p) %in% `Lasi_nige/plat`) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% Lasi_para) ~ "L.paralienus"),
      try((rownames(MDS.p) %in% Lasi_plat) ~ "L.platythorax"),
      try((rownames(MDS.p) %in% Lasi_psam) ~ "L.neglectus"),
      try((rownames(MDS.p) %in% Lasi_umbr) ~ "L.umbratus")
    ))

  color.listuncor <- c(
    "#27AAE1", "#006838", "#8DC63F", "#F9ED32", "#1B75BC", "#2B3990",
    "#92278F", "#EC008C", "#F15A29", "#A97C50", "#00A79D", "blue", "gray40"
  )

  uncorr.mds <- ggplot(data = MDS.col, mapping = aes(x = V1, y = V2, col = Species)) +
    geom_point(cex = 4, alpha = 1, aes(col = Species)) +
    labs(x = "Dimension 1", y = "Dimension 2", col = "Species") +
    theme_bw(base_size = 20) +
    scale_color_manual(values = color.listuncor) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))

  return(uncorr.mds)
}
mds.contam <- vcf.mds.plot(vcf.contam) + ggtitle("Contaminated")
mds.ADR0525 <- vcf.mds.plot(ADR.vcf0525) + ggtitle("ADR 0525")
mds.haploid <- vcf.mds.plot(haploid.vcf) + ggtitle("Haploid")
mds.redo <- vcf.mds.plot(compmap.vcf) + ggtitle("Competitive Mapping")

color.listuncor <- c(
  "#27AAE1", "#006838", "#8DC63F", "#F9ED32", "#1B75BC", "#2B3990",
  "#92278F", "#EC008C", "#F15A29", "#A97C50", "#00A79D", "blue", "gray40"
)

mds.haploid <- ggplot(data = MDS.col.hap, mapping = aes(x = V1, y = -V2, col = Species)) +
  geom_point(cex = 4, alpha = 1, aes(col = Species)) +
  labs(x = "Dimension 1", y = "Dimension 2", col = "Species") +
  theme_bw(base_size = 20) +
  scale_color_manual(values = color.listuncor) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))


setEPS()
postscript("./figures/MDS_compare.eps", width = 40, height = 10)


ggarrange(mds.contam + theme(legend.position = "none") + ggtitle("Original data"),
  mds.redo + theme(legend.position = "none") + ggtitle("Competitive mapping"),
  mds.haploid + theme(legend.position = "none") + ggtitle("Competitive mapping + Haploid"),
  mds.ADR0525 + theme(legend.position = "none") + ggtitle("Competitive mapping + ADR"),
  ncol = 4, nrow = 1, common.legend = T,
  legend.grob = get_legend(mds.contam),
  legend = "right", labels = c("A", "B", "C", "D"), font.label = list(size = 20)
)
dev.off()

# FIGURE 4. MDS ADR species correct ####
allvcf.mds.plot <- function(vcf) {
  gl <- vcfR2genlight(vcf)
  d <- dist(gl)
  m <- as.matrix(d)
  rm(Lasi_)
  rm(`Lasi Chtono`)
  rm(`Lasi jaune`)
  rm(Lasi.Chtono)
  rm(Lasi.jaune)
  rm(Lasi_alie)
  rm(`Lasi_alie gr`)
  rm(Lasi_alie.gr)
  rm(Lasi_brun)
  rm(Lasi_emar)
  rm(Lasi_flav)
  rm(Lasi_fuli)
  rm(Lasi_mixt)
  rm(Lasi_myop)
  rm(Lasi_nige)
  rm(Lasi_umbr)
  rm(Lasi_nige.plat)
  rm(`Lasi_nige/plat`)
  rm(Lasi_para)
  rm(Lasi_plat)
  rm(Lasi_psam)
  ants <- read.csv("data/antmetadata.csv")
  colnames(ants)[1] <- "CATALOGUENUMBER"

  id.list <- read.table("data/all_individuals.txt")
  lasius <- subset(ants, grepl("Las", SPECIESSUMMARY))

  id.las <- lasius[lasius[, "CATALOGUENUMBER"] %in% unlist(id.list$V1), ]
  id.las <- as.data.frame(cbind(id.las$CATALOGUENUMBER, id.las$SPECIESSUMMARY))

  ir <- split(id.las, id.las$V2)

  for (name in (names(ir))) {
    list <- (ir[[name]]$V1)
    assign(name, list)
  }
  Lasi_plat <- c(Lasi_plat, "4457")
  Lasi_flav <- Lasi_flav[Lasi_flav != "9990946"]

  MDS_all <- isoMDS(d = d, k = 2, maxit = 100)

  MDS.p <- as.data.frame(MDS_all$points)
  MDS.col <- MDS.p %>%
    mutate(colors = case_when(
      try((rownames(MDS.p) %in% Lasi) ~ "grey40"),
      try((rownames(MDS.p) %in% `Lasi Chtono`) ~ "grey40"),
      try((rownames(MDS.p) %in% `Lasi jaune`) ~ "grey40"),
      try((rownames(MDS.p) %in% Lasi_alie) ~ "#27AAE1"),
      try((rownames(MDS.p) %in% `Lasi_alie gr`) ~ "grey40"),
      try((rownames(MDS.p) %in% Lasi_brun) ~ "#006838"),
      try((rownames(MDS.p) %in% Lasi_emar) ~ "#8DC63F"),
      try((rownames(MDS.p) %in% Lasi_flav) ~ "#F9ED32"),
      try((rownames(MDS.p) %in% Lasi_fuli) ~ "#1B75BC"),
      try((rownames(MDS.p) %in% Lasi_mixt) ~ "#2B3990"),
      try((rownames(MDS.p) %in% Lasi_myop) ~ "#92278F"),
      try((rownames(MDS.p) %in% Lasi_nige) ~ "#EC008C"),
      try((rownames(MDS.p) %in% `Lasi_nige/plat`) ~ "grey40"),
      try((rownames(MDS.p) %in% Lasi_para) ~ "#F15A29"),
      try((rownames(MDS.p) %in% Lasi_plat) ~ "#A97C50"),
      try((rownames(MDS.p) %in% Lasi_psam) ~ "#00A79D"),
      try((rownames(MDS.p) %in% Lasi_umbr) ~ "#C2B59B")
    ))

  MDS.col <- MDS.col %>%
    mutate(Species = case_when(
      try((rownames(MDS.p) %in% Lasi) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% `Lasi Chtono`) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% `Lasi jaune`) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% Lasi_alie) ~ "L.alienus"),
      try((rownames(MDS.p) %in% `Lasi_alie gr`) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% Lasi_brun) ~ "L.brunneus"),
      try((rownames(MDS.p) %in% Lasi_emar) ~ "L.emarginatus"),
      try((rownames(MDS.p) %in% Lasi_flav) ~ "L.flavus"),
      try((rownames(MDS.p) %in% Lasi_fuli) ~ "L.fuliginosus"),
      try((rownames(MDS.p) %in% Lasi_mixt) ~ "L.mixtus"),
      try((rownames(MDS.p) %in% Lasi_myop) ~ "L.myops"),
      try((rownames(MDS.p) %in% Lasi_nige) ~ "L.niger"),
      try((rownames(MDS.p) %in% `Lasi_nige/plat`) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% Lasi_para) ~ "L.paralienus"),
      try((rownames(MDS.p) %in% Lasi_plat) ~ "L.platythorax"),
      try((rownames(MDS.p) %in% Lasi_psam) ~ "L.neglectus"),
      try((rownames(MDS.p) %in% Lasi_umbr) ~ "L.umbratus")
    ))

  color.listuncor <- c(
    "#27AAE1", "#006838", "#8DC63F", "#F9ED32", "#1B75BC", "#2B3990",
    "#92278F", "#EC008C", "#F15A29", "#A97C50", "#00A79D", "blue", "gray40"
  )

  uncorr.mds <- ggplot(data = MDS.col, mapping = aes(x = V1, y = V2, col = Species)) +
    geom_point(cex = 6, alpha = 0.4, aes(col = Species)) +
    labs(x = "Dimension 1", y = "Dimension 2", col = "Species") +
    theme_bw(base_size = 20) +
    scale_color_manual(values = color.listuncor) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))

  # corrected ids##
  rm(Lasi_)
  rm(`Lasi Chtono`)
  rm(`Lasi jaune`)
  rm(Lasi.Chtono)
  rm(Lasi.jaune)
  rm(Lasi_alie)
  rm(`Lasi_alie gr`)
  rm(Lasi_alie.gr)
  rm(Lasi_brun)
  rm(Lasi_emar)
  rm(Lasi_flav)
  rm(Lasi_fuli)
  rm(Lasi_mixt)
  rm(Lasi_myop)
  rm(Lasi_nige)
  rm(Lasi_umbr)
  rm(Lasi_nige.plat)
  rm(`Lasi_nige/plat`)
  rm(Lasi_para)
  rm(Lasi_plat)
  rm(Lasi_psam)

  correctedID <- read.csv("data/correctID1.csv", header = T)
  for (name in (colnames(correctedID))) {
    list <- na.exclude(correctedID[[name]])
    assign(name, list)
  }

  MDS.col <- MDS.col %>%
    mutate(New.Species = case_when(
      try((rownames(MDS.p) %in% Lasi) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% Lasi_alie) ~ "L.alienus"),
      try((rownames(MDS.p) %in% Lasi_brun) ~ "L.brunneus"),
      try((rownames(MDS.p) %in% Lasi_emar) ~ "L.emarginatus"),
      try((rownames(MDS.p) %in% Lasi_flav) ~ "L.flavus"),
      try((rownames(MDS.p) %in% Lasi_fuli) ~ "L.fuliginosus"),
      try((rownames(MDS.p) %in% Lasi_mixt) ~ "L.mixtus"),
      try((rownames(MDS.p) %in% Lasi_myop) ~ "L.myops"),
      try((rownames(MDS.p) %in% Lasi_nige) ~ "L.niger"),
      try((rownames(MDS.p) %in% Lasi_para) ~ "L.paralienus"),
      try((rownames(MDS.p) %in% Lasi_plat) ~ "L.platythorax"),
      try((rownames(MDS.p) %in% Lasi_psam) ~ "L.neglectus"),
      try((rownames(MDS.p) %in% Lasi_umbr) ~ "L.umbratus")
    ))
  MDS.col <- MDS.col %>%
    mutate(colors = case_when(
      try((rownames(MDS.p) %in% Lasi_alie) ~ "#27AAE1"),
      try((rownames(MDS.p) %in% Lasi_brun) ~ "#006838"),
      try((rownames(MDS.p) %in% Lasi_emar) ~ "#8DC63F"),
      try((rownames(MDS.p) %in% Lasi_flav) ~ "#F9ED32"),
      try((rownames(MDS.p) %in% Lasi_fuli) ~ "#1B75BC"),
      try((rownames(MDS.p) %in% Lasi_mixt) ~ "#2B3990"),
      try((rownames(MDS.p) %in% Lasi_myop) ~ "#92278F"),
      try((rownames(MDS.p) %in% Lasi_nige) ~ "#EC008C"),
      try((rownames(MDS.p) %in% Lasi_para) ~ "#F15A29"),
      try((rownames(MDS.p) %in% Lasi_plat) ~ "#A97C50"),
      try((rownames(MDS.p) %in% Lasi_psam) ~ "#00A79D"),
      try((rownames(MDS.p) %in% Lasi_umbr) ~ "#C2B59B")
    ))

  MDS.col <- MDS.col %>%
    mutate(Spec.change = ifelse(MDS.col$Species == MDS.col$New.Species, "Same as Morphological Identification", "Different from Morphological Identification"))

  MDS.col$idnum <- rownames(MDS.col)

  MDS.col <- MDS.col %>%
    mutate(alph = case_when(
      try((Spec.change == "Same as Morphological Identification") ~ 0.4),
      try((Spec.change == "Different from Morphological Identification") ~ 1)
    ))


  col <- as.character(MDS.col$colors)
  names(col) <- as.character(MDS.col$New.Species)

  corr.mds <- ggplot(data = MDS.col, mapping = aes(x = V1, y = V2)) +
    geom_point(cex = 4, aes(shape = Spec.change, col = New.Species), alpha = 1) +
    labs(
      x = "Dimension 1", y = "Dimension 2", alpha = "Final Species Identification",
      shape = "Final Species Identification", col = "Species"
    ) +
    theme_bw(base_size = 20) +
    scale_color_manual(values = col) +
    scale_shape_manual(values = c(17, 16)) +
    guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))

  return(corr.mds)
}
allmds.ADR0525 <- allvcf.mds.plot(ADR.vcf0525)


ggarrange(mds.contam + theme(legend.position = "none") + ggtitle("Original data"),
  mds.redo + theme(legend.position = "none") + ggtitle("Competitive mapping"),
  mds.haploid + theme(legend.position = "none") + ggtitle("Competitive mapping + Haploid"),
  mds.ADR0525 + theme(legend.position = "none") + ggtitle("Competitive mapping + ADR"),
  allmds.ADR0525 + theme(legend.position = "none") + ggtitle("Competitive mapping + ADR - Species correction"),
  ncol = 2, nrow = 3, common.legend = T,
  legend.grob = get_legend(mds.contam),
  legend = "right", labels = c("A", "B", "C", "D", "E"), font.label = list(size = 20)
)


# FIGURE 5. COI + ADMIX comparison ####
correctedID <- read.csv("data/correctID1.csv", header = T)
for (name in (colnames(correctedID))) {
  list <- na.exclude(correctedID[[name]])
  assign(name, list)
}


# contam compmap adr haploid
admix.plot <- function(DECON) {
  if (DECON == "contam") {
    order <- read.table("data/contam/cleaned.fam", as.is = T, header = F)
    qfile <- read.table("data/contam/cleaned.9.Q", as.is = T, header = F)
    species.list <- c(
      "spL.platythorax", "spL.myops", "spL.mixtus", "spL.emarginatus",
      "spL.brunneus", "spL.fuliginosus", "spL.flavus", "spL.niger", "spL.paralienus", "id"
    )
  }
  if (DECON == "compmap") {
    order <- read.table("data/compmap/cleaned.fam", as.is = T, header = F)
    qfile <- read.table("data/compmap/cleaned.9.Q", as.is = T, header = F)
    species.list <- c(
      "spL.platythorax", "spL.fuliginosus", "spL.flavus", "spL.niger",
      "spL.mixtus", "spL.paralienus", "spL.emarginatus", "spL.myops", "spL.brunneus", "id"
    )
  }
  if (DECON == "adr") {
    order <- read.table("data/adr/cleaned.fam", as.is = T, header = F)
    qfile <- read.table("data/adr/cleaned.9.Q", as.is = T, header = F)
    species.list <- c(
      "spL.paralienus", "spL.brunneus", "spL.platythorax", "spL.niger", "spL.emarginatus",
      "spL.mixtus", "spL.flavus", "spL.myops", "spL.fuliginosus", "id"
    )
  }
  if (DECON == "haploid") {
    order <- read.table("data/haploid/cleaned.fam", as.is = T, header = F)
    order$V2 <- gsub(".sorted.filtered", "", order$V2)
    qfile <- read.table("data/haploid/cleaned.9.Q", as.is = T, header = F)
    species.list <- c(
      "spL.fuliginosus", "spL.flavus", "spL.platythorax", "spL.emarginatus", "spL.niger",
      "spL.mixtus", "spL.paralienus", "spL.brunneus", "spL.myops", "id"
    )
  }

  species.color <- c(
    "#006838", "#8DC63F", "#F9ED32", "#1B75BC", "#2B3990", "#92278F", "#EC008C",
    "#F15A29", "#A97C50", "#00A79D", "grey40", "#27AAE1"
  )

  order$V1 <- rownames(order)
  row.names(order) <- order$V2

  ind2pop <- order %>%
    mutate(spec = case_when(
      try((rownames(order) %in% Lasi) ~ "Lasi"),
      try((rownames(order) %in% Lasi_alie) ~ "Lasi_alie"),
      try((rownames(order) %in% Lasi_brun) ~ "Lasi_brun"),
      try((rownames(order) %in% Lasi_emar) ~ "Lasi_emar"),
      try((rownames(order) %in% Lasi_flav) ~ "Lasi_flav"),
      try((rownames(order) %in% Lasi_fuli) ~ "Lasi_fuli"),
      try((rownames(order) %in% Lasi_mixt) ~ "Lasi_mixt"),
      try((rownames(order) %in% Lasi_myop) ~ "Lasi_myop"),
      try((rownames(order) %in% Lasi_nige) ~ "Lasi_nige"),
      try((rownames(order) %in% Lasi_para) ~ "Lasi_para"),
      try((rownames(order) %in% Lasi_plat) ~ "Lasi_plat"),
      try((rownames(order) %in% Lasi_psam) ~ "Lasi_negl"),
      try((rownames(order) %in% Lasi_umbr) ~ "Lasi_umbr")
    ))

  rownames(qfile) <- order$V2

  qfilesp <- qfile %>%
    mutate(spec = case_when(
      try((rownames(qfile) %in% Lasi_alie) ~ "Lasi_alie"),
      try((rownames(qfile) %in% Lasi_brun) ~ "Lasi_brun"),
      try((rownames(qfile) %in% Lasi_emar) ~ "Lasi_emar"),
      try((rownames(qfile) %in% Lasi_flav) ~ "Lasi_flav"),
      try((rownames(qfile) %in% Lasi_fuli) ~ "Lasi_fuli"),
      try((rownames(qfile) %in% Lasi_mixt) ~ "Lasi_mixt"),
      try((rownames(qfile) %in% Lasi_myop) ~ "Lasi_myop"),
      try((rownames(qfile) %in% Lasi_nige) ~ "Lasi_nige"),
      try((rownames(qfile) %in% Lasi_para) ~ "Lasi_para"),
      try((rownames(qfile) %in% Lasi_plat) ~ "Lasi_plat"),
      try((rownames(qfile) %in% Lasi_psam) ~ "Lasi_negl"),
      try((rownames(qfile) %in% Lasi_umbr) ~ "Lasi_umbr")
    ))

  colnames(qfilesp) <- species.list

  qfilesp$id.num <- rownames(qfilesp)
  longqfile <- pivot_longer(qfilesp,
    cols = starts_with("sp"),
    names_to = "Species",
    names_prefix = "sp",
    values_to = "prop",
  )


  id.orders <- unique(longqfile[order(longqfile$id), ]$id.num)

  campoqfile <- data.frame(rbind(
    c("Campo", "Campo", "L.brunneus", 0),
    c("Campo", "Campo", "L.emarginatus", 0),
    c("Campo", "Campo", "L.flavus", 0),
    c("Campo", "Campo", "L.fuliginosus", 0),
    c("Campo", "Campo", "L.mixtus", 0),
    c("Campo", "Campo", "L.myops", 0),
    c("Campo", "Campo", "L.niger", 0),
    c("Campo", "Campo", "L.paralienus", 0),
    c("Campo", "Campo", "L.platythorax", 0)
  ))
  colnames(campoqfile) <- colnames(longqfile)

  tree <- read.tree("data/LasiusCOI.MAFFT.treefile")
  tips <- as.data.frame(tree$tip.label)
  subbed.longQ <- filter(longqfile, id.num %in% tree$tip.label)

  colnames(campoqfile) <- colnames(subbed.longQ)
  subbed.longQ <- rbind(subbed.longQ, campoqfile)
  subbed.longQ$prop <- as.numeric(subbed.longQ$prop)

  id.orders <- unique(subbed.longQ[order(subbed.longQ$id), ]$id.num)
  subbed.longQ$id.num <- (factor(subbed.longQ$id.num, level = id.orders))

  qfilem <- qfilesp
  qfilem$max <- pmax(
    qfilem$spL.mixtus, qfilem$spL.flavus, qfilem$spL.platythorax, qfilem$spL.myops, qfilem$spL.niger,
    qfilem$spL.brunneus, qfilem$spL.fuliginosus, qfilem$spL.paralienus, qfilem$spL.emarginatus
  )
  blackhyb <- qfilem[qfilem$max < 0.96875, ] # 1,2,3,4 backcrosses away
  grayhyb <- qfilem[between(qfilem$max, 0.96875, 0.984375), ] # 11 5 backcrosses away
  blackhyb <- blackhyb[, !names(blackhyb) %in% c("max")]
  grayhyb <- grayhyb[, !names(grayhyb) %in% c("max")]
  hybrids <- c(rownames(blackhyb), rownames(grayhyb))
  hybrids <- hybrids[!hybrids %in% c("20350", "12831", "10008")] # removes alie and negl

  subbed.longQ <- subbed.longQ %>%
    mutate(hybrid = ifelse(subbed.longQ$id.num %in% hybrids, "hybrid", NA))


  HYBPERC <- round((100 * length(hybrids)) / length(rownames(qfilesp)), digits = 2)

  # negl and alie not hybrids
  subbed.longQ$hybrid[which(subbed.longQ$id == "Lasi_negl")] <- NA
  subbed.longQ$hybrid[which(subbed.longQ$id == "Lasi_alie")] <- NA


  Hfulladmix <- ggplot(
    data = subbed.longQ,
    aes(
      fill = Species, y = prop, x = id.num,
      label = ifelse(hybrid == "hybrid", "*", "")
    )
  ) +
    geom_bar(position = "stack", stat = "identity", width = 1) +
    theme_pubr(base_size = 20) +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 15), plot.subtitle = element_text(size = 15),
      axis.title.x = element_text(size = 15)
    ) +
    xlab("") +
    ylab("Proportion of ancestry") +
    scale_fill_manual(values = unique(species.color)) +
    coord_flip() +
    geom_text(size = 10, y = 0) +
    labs(subtitle = paste(HYBPERC, "% admixed individuals", sep = ""))

  return(Hfulladmix)
}


contam.admix <- admix.plot("contam")
compmap.admix <- admix.plot("compmap")
adr.admix <- admix.plot("adr")
haploid.admix <- admix.plot("haploid")


# tree####
tree <- read.tree("data/LasiusCOI.MAFFT.treefile")

tips <- as.data.frame(tree$tip.label)
tips <- tips %>%
  mutate(spec = case_when(
    try((tips$`tree$tip.label` %in% Lasi_alie) ~ "Lasi_alie"),
    try((tips$`tree$tip.label` %in% Lasi_brun) ~ "Lasi_brun"),
    try((tips$`tree$tip.label` %in% Lasi_emar) ~ "Lasi_emar"),
    try((tips$`tree$tip.label` %in% Lasi_flav) ~ "Lasi_flav"),
    try((tips$`tree$tip.label` %in% Lasi_fuli) ~ "Lasi_fuli"),
    try((tips$`tree$tip.label` %in% Lasi_mixt) ~ "Lasi_mixt"),
    try((tips$`tree$tip.label` %in% Lasi_myop) ~ "Lasi_myop"),
    try((tips$`tree$tip.label` %in% Lasi_nige) ~ "Lasi_nige"),
    try((tips$`tree$tip.label` %in% Lasi_para) ~ "Lasi_para"),
    try((tips$`tree$tip.label` %in% Lasi_plat) ~ "Lasi_plat"),
    try((tips$`tree$tip.label` %in% Lasi_psam) ~ "Lasi_negl"),
    try((tips$`tree$tip.label` %in% Lasi_umbr) ~ "Lasi_umbr"),
    try((tips$`tree$tip.label` == "Campo") ~ "Campo")
  ))


dropped <- setdiff(tree$tip.label, MDS.col.contam$idnum)
subtree <- drop.tip(tree, dropped)

treedata <- as.data.frame(cbind(
  ID = MDS.col.contam$idnum,
  Species = MDS.col.contam$New.Species, Colors = MDS.col.contam$colors
))

p <- ggtree(subtree) %<+% treedata

color.list3 <- c(
  "#27AAE1", "#006838", "#8DC63F", "#F9ED32", "#1B75BC", "#2B3990", "#92278F", "#00A79D",
  "#EC008C", "#F15A29", "#A97C50"
)
p2 <- p + geom_tippoint(aes(color = Species), size = 1, shape = 15) +
  scale_color_manual(values = unique(color.list3)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 16), legend.title = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 5)))
p2


# redo tree ####
gt <- extract.gt(ADR.vcf0525, "GT")
individuals <- dimnames(gt)[[2]]
order <- read.table("data/adr/cleaned.fam", as.is = T, header = F)
qfile <- read.table("data/adr/cleaned.9.Q", as.is = T, header = F)
species.list <- c(
  "spL.paralienus", "spL.brunneus", "spL.platythorax", "spL.niger", "spL.emarginatus",
  "spL.mixtus", "spL.flavus", "spL.myops", "spL.fuliginosus", "id"
)
order$V1 <- rownames(order)
row.names(order) <- order$V2

ind2pop <- order %>%
  mutate(spec = case_when(
    try((rownames(order) %in% Lasi) ~ "Lasi"),
    try((rownames(order) %in% Lasi_alie) ~ "Lasi_alie"),
    try((rownames(order) %in% Lasi_brun) ~ "Lasi_brun"),
    try((rownames(order) %in% Lasi_emar) ~ "Lasi_emar"),
    try((rownames(order) %in% Lasi_flav) ~ "Lasi_flav"),
    try((rownames(order) %in% Lasi_fuli) ~ "Lasi_fuli"),
    try((rownames(order) %in% Lasi_mixt) ~ "Lasi_mixt"),
    try((rownames(order) %in% Lasi_myop) ~ "Lasi_myop"),
    try((rownames(order) %in% Lasi_nige) ~ "Lasi_nige"),
    try((rownames(order) %in% Lasi_para) ~ "Lasi_para"),
    try((rownames(order) %in% Lasi_plat) ~ "Lasi_plat"),
    try((rownames(order) %in% Lasi_psam) ~ "Lasi_negl"),
    try((rownames(order) %in% Lasi_umbr) ~ "Lasi_umbr")
  ))

rownames(qfile) <- order$V2

qfilesp <- qfile %>%
  mutate(spec = case_when(
    try((rownames(qfile) %in% Lasi_alie) ~ "Lasi_alie"),
    try((rownames(qfile) %in% Lasi_brun) ~ "Lasi_brun"),
    try((rownames(qfile) %in% Lasi_emar) ~ "Lasi_emar"),
    try((rownames(qfile) %in% Lasi_flav) ~ "Lasi_flav"),
    try((rownames(qfile) %in% Lasi_fuli) ~ "Lasi_fuli"),
    try((rownames(qfile) %in% Lasi_mixt) ~ "Lasi_mixt"),
    try((rownames(qfile) %in% Lasi_myop) ~ "Lasi_myop"),
    try((rownames(qfile) %in% Lasi_nige) ~ "Lasi_nige"),
    try((rownames(qfile) %in% Lasi_para) ~ "Lasi_para"),
    try((rownames(qfile) %in% Lasi_plat) ~ "Lasi_plat"),
    try((rownames(qfile) %in% Lasi_psam) ~ "Lasi_negl"),
    try((rownames(qfile) %in% Lasi_umbr) ~ "Lasi_umbr")
  ))

colnames(qfilesp) <- species.list

qfilesp$id.num <- rownames(qfilesp)
qfilesp.subbed <- qfilesp[(rownames(qfilesp) %in% subtree[["tip.label"]]), ]

dropped <- setdiff(tree$tip.label, MDS.col.adr$idnum)
subtree <- drop.tip(tree, dropped)


setEPS()
postscript("./figures/g_phyloadmix_contam.eps", width = 100, height = 30)


par(mfrow = c(2, 1), mai = c(0.1, 0.6, 0.1, 0.1))


barplot(t(as.matrix(qfilesp.subbed))[, match(subtree$tip.label, rownames(qfilesp.subbed))],
  las = 2, border = NA, space = 0, cex.names = 0.5,
  col = c("#F15A29", "#006838", "#A97C50", "#EC008C", "#8DC63F", "#2B3990", "#F9ED32", "#92278F", "#1B75BC")
)
plot.phylo(subtree,
  align.tip.label = T, cex = 0.4, direction = "upwards", root.edge = T,
  tip.color = MDS.col.contam$oldcolors[match(subtree$tip.label, MDS.col.contam$idnum)]
)


dev.off()


# combo####


dev.off()

setEPS()
postscript("./figures/treeadmix_compare.eps", width = 25, height = 18)


print((contam.admix + labs(tag = "B")) |>
  insert_left((p2 + labs(tag = "A"))) |>
  insert_right((compmap.admix + labs(tag = "C"))) |>
  insert_right((haploid.admix + labs(tag = "D"))) |>
  insert_right((adr.admix + labs(tag = "E")))) & theme(
  legend.position = "bottom",
  plot.tag = element_text(size = 20, face = "bold")
)

dev.off()


# SUPP Cryptic Species MDS ####
MDS.col <- vcf.mds(ADR.vcf0525)
gl <- vcfR2genlight(ADR.vcf0525)
d <- dist(gl)
m <- as.matrix(d)

crytpic.mds.plot <- function(SPECIES) {
  cluster <- m |>
    as.data.frame() |>
    filter(rownames(m) %in% SPECIES)
  cluster <- cluster |>
    select(rownames(cluster)) |>
    as.matrix()

  clusterMDS_all <- isoMDS(d = cluster, k = 2, maxit = 100)
  MDS.p <- as.data.frame(clusterMDS_all$points)


  rm(Lasi_)
  rm(`Lasi Chtono`)
  rm(`Lasi jaune`)
  rm(Lasi.Chtono)
  rm(Lasi.jaune)
  rm(Lasi_alie)
  rm(`Lasi_alie gr`)
  rm(Lasi_alie.gr)
  rm(Lasi_brun)
  rm(Lasi_emar)
  rm(Lasi_flav)
  rm(Lasi_fuli)
  rm(Lasi_mixt)
  rm(Lasi_myop)
  rm(Lasi_nige)
  rm(Lasi_umbr)
  rm(Lasi_nige.plat)
  rm(`Lasi_nige/plat`)
  rm(Lasi_para)
  rm(Lasi_plat)
  rm(Lasi_psam)

  correctedID <- read.csv("data/correctID1.csv", header = T)
  for (name in (colnames(correctedID))) {
    list <- na.exclude(correctedID[[name]])
    assign(name, list)
  }


  MDS.col <- MDS.p %>%
    mutate(New.Species = case_when(
      try((rownames(MDS.p) %in% Lasi) ~ "Unidentified Lasius"),
      try((rownames(MDS.p) %in% Lasi_alie) ~ "L.alienus"),
      try((rownames(MDS.p) %in% Lasi_brun) ~ "L.brunneus"),
      try((rownames(MDS.p) %in% Lasi_emar) ~ "L.emarginatus"),
      try((rownames(MDS.p) %in% Lasi_flav) ~ "L.flavus"),
      try((rownames(MDS.p) %in% Lasi_fuli) ~ "L.fuliginosus"),
      try((rownames(MDS.p) %in% Lasi_mixt) ~ "L.mixtus"),
      try((rownames(MDS.p) %in% Lasi_myop) ~ "L.myops"),
      try((rownames(MDS.p) %in% Lasi_nige) ~ "L.niger"),
      try((rownames(MDS.p) %in% Lasi_para) ~ "L.paralienus"),
      try((rownames(MDS.p) %in% Lasi_plat) ~ "L.platythorax"),
      try((rownames(MDS.p) %in% Lasi_psam) ~ "L.neglectus"),
      try((rownames(MDS.p) %in% Lasi_umbr) ~ "L.umbratus")
    ))
  MDS.col <- MDS.col %>%
    mutate(colors = case_when(
      try((rownames(MDS.p) %in% Lasi_alie) ~ "#27AAE1"),
      try((rownames(MDS.p) %in% Lasi_brun) ~ "#006838"),
      try((rownames(MDS.p) %in% Lasi_emar) ~ "#8DC63F"),
      try((rownames(MDS.p) %in% Lasi_flav) ~ "#F9ED32"),
      try((rownames(MDS.p) %in% Lasi_fuli) ~ "#1B75BC"),
      try((rownames(MDS.p) %in% Lasi_mixt) ~ "#2B3990"),
      try((rownames(MDS.p) %in% Lasi_myop) ~ "#92278F"),
      try((rownames(MDS.p) %in% Lasi_nige) ~ "#EC008C"),
      try((rownames(MDS.p) %in% Lasi_para) ~ "#F15A29"),
      try((rownames(MDS.p) %in% Lasi_plat) ~ "#A97C50"),
      try((rownames(MDS.p) %in% Lasi_psam) ~ "#00A79D"),
      try((rownames(MDS.p) %in% Lasi_umbr) ~ "#C2B59B")
    ))

  MDS.col$idnum <- rownames(MDS.col)

  corr.mds <- ggplot(data = MDS.col, mapping = aes(x = V1, y = V2)) +
    geom_point(cex = 4, alpha = 1, color = MDS.col$colors) +
    labs(
      x = "Dimension 1", y = "Dimension 2", col = "Species"
    ) +
    theme_bw(base_size = 20) +
    ggtitle(MDS.col$New.Species[1])


  return(corr.mds)
}

cryp.brun <- crytpic.mds.plot(Lasi_brun)
cryp.emar <- crytpic.mds.plot(Lasi_emar)
cryp.flav <- crytpic.mds.plot(Lasi_flav)
cryp.fuli <- crytpic.mds.plot(Lasi_fuli)

cryp.nige <- crytpic.mds.plot(Lasi_nige)
cryp.para <- crytpic.mds.plot(Lasi_para)
cryp.plat <- crytpic.mds.plot(Lasi_plat)


setEPS()
postscript("./figures/SUPP2_cryptic_species_MDS.eps", width = 40, height = 20)

ggarrange(cryp.brun, cryp.emar, cryp.flav, cryp.fuli,
  cryp.nige, cryp.para, cryp.plat,
  ncol = 4, nrow = 2, legend = "none"
)
dev.off()


# Supplementary Fig 1 Comp map genome mapping ####
summap <- read.table(file = "data/LASIallsummap.txt", header = T)


name <- deparse(substitute(summap))
colnames(summap) <- gsub("X", "", as.character(colnames(summap)))

summap1 <- melt(summap, id.vars = "spec")
summap1[is.na(summap1)] <- 0 # change NAs to 0
summap1$value <- as.numeric(summap1$value)

summapP <- ggplot(summap1, aes(spec, value)) +
  geom_jitter(size = 4, alpha = 1, width = 0.3) +
  labs(
    title = "Summary of Competitive Mapping",
    x = "Genus", y = "Proportion of Mapped Reads"
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1)) +
  theme_bw(base_size = 20)


setEPS()
postscript("./figures/SUPP1_summap.eps", width = 25, height = 18)

summapP

dev.off()


# SUPP full COI tree####

tree <- read.tree("data/LasiusCOI.MAFFT.treefile")
allindv <- read.table("./data/all_individuals.txt")

dropped <- setdiff(tree$tip.label, allindv$V1)
subtree <- drop.tip(tree, dropped)

sub.tips <- as.data.frame(subtree$tip.label)
sub.tips <- sub.tips %>%
  mutate(spec = case_when(
    try((sub.tips$`subtree$tip.label` %in% Lasi_alie) ~ "L.alienus"),
    try((sub.tips$`subtree$tip.label` %in% Lasi_brun) ~ "L.brunneus"),
    try((sub.tips$`subtree$tip.label` %in% Lasi_emar) ~ "L.emarginatus"),
    try((sub.tips$`subtree$tip.label` %in% Lasi_flav) ~ "L.flavus"),
    try((sub.tips$`subtree$tip.label` %in% Lasi_fuli) ~ "L.fuliginosus"),
    try((sub.tips$`subtree$tip.label` %in% Lasi_mixt) ~ "L.mixtus"),
    try((sub.tips$`subtree$tip.label` %in% Lasi_myop) ~ "L.myops"),
    try((sub.tips$`subtree$tip.label` %in% Lasi_nige) ~ "L.niger"),
    try((sub.tips$`subtree$tip.label` %in% Lasi_para) ~ "L.paralienus"),
    try((sub.tips$`subtree$tip.label` %in% Lasi_plat) ~ "L.platythorax"),
    try((sub.tips$`subtree$tip.label` %in% Lasi_psam) ~ "L.neglectus"),
    try((sub.tips$`subtree$tip.label` %in% Lasi_umbr) ~ "L.umbratus"),
    try((sub.tips$`subtree$tip.label` == "Campo") ~ "Campo")
  ))

sub.tips$spec <- replace_na(sub.tips$spec, "Unidentified Lasius")

treedata <- as.data.frame(cbind(
  ID = sub.tips$`subtree$tip.label`,
  Species = sub.tips$spec
))

p0 <- ggtree(subtree) %<+% treedata

color.list3 <- c(
  "#27AAE1", "#006838", "#8DC63F", "#F9ED32", "#1B75BC", "#2B3990", "#92278F", "#EC008C",
  "#F15A29", "#A97C50", "#00A79D", "blue", "gray40"
)
p3 <- p0 + geom_tippoint(aes(color = Species), size = 2, shape = 15) +
  scale_color_manual(values = unique(color.list3)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 16), legend.title = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 5)))
p3

setEPS()
postscript("./figures/SUPP3_fullCOI.eps", width = 25, height = 18)

p3
dev.off()


# SUPP hybrid only admix####

hyb.admix.plot <- function(DECON) {
  if (DECON == "contam") {
    order <- read.table("data/contam/cleaned.fam", as.is = T, header = F)
    qfile <- read.table("data/contam/cleaned.9.Q", as.is = T, header = F)
    species.list <- c(
      "spL.platythorax", "spL.myops", "spL.mixtus", "spL.emarginatus",
      "spL.brunneus", "spL.fuliginosus", "spL.flavus", "spL.niger", "spL.paralienus", "id"
    )
  }
  if (DECON == "compmap") {
    order <- read.table("data/compmap/cleaned.fam", as.is = T, header = F)
    qfile <- read.table("data/compmap/cleaned.9.Q", as.is = T, header = F)
    species.list <- c(
      "spL.platythorax", "spL.fuliginosus", "spL.flavus", "spL.niger",
      "spL.mixtus", "spL.paralienus", "spL.emarginatus", "spL.myops", "spL.brunneus", "id"
    )
  }
  if (DECON == "adr") {
    order <- read.table("data/adr/cleaned.fam", as.is = T, header = F)
    qfile <- read.table("data/adr/cleaned.9.Q", as.is = T, header = F)
    species.list <- c(
      "spL.paralienus", "spL.brunneus", "spL.platythorax", "spL.niger", "spL.emarginatus",
      "spL.mixtus", "spL.flavus", "spL.myops", "spL.fuliginosus", "id"
    )
  }
  if (DECON == "haploid") {
    order <- read.table("data/haploid/cleaned.fam", as.is = T, header = F)
    order$V2 <- gsub(".sorted.filtered", "", order$V2)
    qfile <- read.table("data/haploid/cleaned.9.Q", as.is = T, header = F)
    species.list <- c(
      "spL.fuliginosus", "spL.flavus", "spL.platythorax", "spL.emarginatus", "spL.niger",
      "spL.mixtus", "spL.paralienus", "spL.brunneus", "spL.myops", "id"
    )
  }

  species.color <- c(
    "#006838", "#8DC63F", "#F9ED32", "#1B75BC", "#2B3990", "#92278F", "#EC008C",
    "#F15A29", "#A97C50", "#00A79D", "grey40", "#27AAE1"
  )

  order$V1 <- rownames(order)
  row.names(order) <- order$V2

  ind2pop <- order %>%
    mutate(spec = case_when(
      try((rownames(order) %in% Lasi) ~ "Lasi"),
      try((rownames(order) %in% Lasi_alie) ~ "Lasi_alie"),
      try((rownames(order) %in% Lasi_brun) ~ "Lasi_brun"),
      try((rownames(order) %in% Lasi_emar) ~ "Lasi_emar"),
      try((rownames(order) %in% Lasi_flav) ~ "Lasi_flav"),
      try((rownames(order) %in% Lasi_fuli) ~ "Lasi_fuli"),
      try((rownames(order) %in% Lasi_mixt) ~ "Lasi_mixt"),
      try((rownames(order) %in% Lasi_myop) ~ "Lasi_myop"),
      try((rownames(order) %in% Lasi_nige) ~ "Lasi_nige"),
      try((rownames(order) %in% Lasi_para) ~ "Lasi_para"),
      try((rownames(order) %in% Lasi_plat) ~ "Lasi_plat"),
      try((rownames(order) %in% Lasi_psam) ~ "Lasi_negl"),
      try((rownames(order) %in% Lasi_umbr) ~ "Lasi_umbr")
    ))

  rownames(qfile) <- order$V2

  qfilesp <- qfile %>%
    mutate(spec = case_when(
      try((rownames(qfile) %in% Lasi_alie) ~ "Lasi_alie"),
      try((rownames(qfile) %in% Lasi_brun) ~ "Lasi_brun"),
      try((rownames(qfile) %in% Lasi_emar) ~ "Lasi_emar"),
      try((rownames(qfile) %in% Lasi_flav) ~ "Lasi_flav"),
      try((rownames(qfile) %in% Lasi_fuli) ~ "Lasi_fuli"),
      try((rownames(qfile) %in% Lasi_mixt) ~ "Lasi_mixt"),
      try((rownames(qfile) %in% Lasi_myop) ~ "Lasi_myop"),
      try((rownames(qfile) %in% Lasi_nige) ~ "Lasi_nige"),
      try((rownames(qfile) %in% Lasi_para) ~ "Lasi_para"),
      try((rownames(qfile) %in% Lasi_plat) ~ "Lasi_plat"),
      try((rownames(qfile) %in% Lasi_psam) ~ "Lasi_negl"),
      try((rownames(qfile) %in% Lasi_umbr) ~ "Lasi_umbr")
    ))

  colnames(qfilesp) <- species.list

  qfilesp$id.num <- rownames(qfilesp)
  longqfile <- pivot_longer(qfilesp,
    cols = starts_with("sp"),
    names_to = "Species",
    names_prefix = "sp",
    values_to = "prop",
  )


  id.orders <- unique(longqfile[order(longqfile$id), ]$id.num)

  campoqfile <- data.frame(rbind(
    c("Campo", "Campo", "L.brunneus", 0),
    c("Campo", "Campo", "L.emarginatus", 0),
    c("Campo", "Campo", "L.flavus", 0),
    c("Campo", "Campo", "L.fuliginosus", 0),
    c("Campo", "Campo", "L.mixtus", 0),
    c("Campo", "Campo", "L.myops", 0),
    c("Campo", "Campo", "L.niger", 0),
    c("Campo", "Campo", "L.paralienus", 0),
    c("Campo", "Campo", "L.platythorax", 0)
  ))
  colnames(campoqfile) <- colnames(longqfile)

  tree <- read.tree("data/LasiusCOI.MAFFT.treefile")
  tips <- as.data.frame(tree$tip.label)
  subbed.longQ <- filter(longqfile, id.num %in% tree$tip.label)

  colnames(campoqfile) <- colnames(subbed.longQ)
  subbed.longQ <- rbind(subbed.longQ, campoqfile)
  subbed.longQ$prop <- as.numeric(subbed.longQ$prop)

  id.orders <- unique(subbed.longQ[order(subbed.longQ$id), ]$id.num)
  subbed.longQ$id.num <- (factor(subbed.longQ$id.num, level = id.orders))

  qfilem <- qfilesp
  qfilem$max <- pmax(
    qfilem$spL.mixtus, qfilem$spL.flavus, qfilem$spL.platythorax, qfilem$spL.myops, qfilem$spL.niger,
    qfilem$spL.brunneus, qfilem$spL.fuliginosus, qfilem$spL.paralienus, qfilem$spL.emarginatus
  )
  blackhyb <- qfilem[qfilem$max < 0.96875, ] # 1,2,3,4 backcrosses away
  grayhyb <- qfilem[between(qfilem$max, 0.96875, 0.984375), ] # 11 5 backcrosses away
  blackhyb <- blackhyb[, !names(blackhyb) %in% c("max")]
  grayhyb <- grayhyb[, !names(grayhyb) %in% c("max")]
  hybrids <- c(rownames(blackhyb), rownames(grayhyb))
  hybrids <- hybrids[!hybrids %in% c("20350", "12831", "10008")] # removes alie and negl

  subbed.longQ <- subbed.longQ %>%
    mutate(hybrid = ifelse(subbed.longQ$id.num %in% hybrids, "hybrid", NA))


  HYBPERC <- round((100 * length(hybrids)) / length(rownames(qfilesp)), digits = 2)

  # negl and alie not hybrids
  subbed.longQ$hybrid[which(subbed.longQ$id == "Lasi_negl")] <- NA
  subbed.longQ$hybrid[which(subbed.longQ$id == "Lasi_alie")] <- NA

  subbed.longQ |> filter(id == "Lasi_alie")
  subbed.longQ <- subbed.longQ |> filter(hybrid == "hybrid")

  Hfulladmix <- ggplot(
    data = subbed.longQ,
    aes(
      fill = Species, y = prop, x = id.num,
    )
  ) +
    geom_bar(position = "stack", stat = "identity", width = 1) +
    theme_pubr(base_size = 20) +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 15), plot.subtitle = element_text(size = 15),
      axis.title.x = element_text(size = 15)
    ) +
    xlab("") +
    ylab("Proportion of ancestry") +
    scale_fill_manual(values = unique(species.color)) +
    coord_flip()


  return(Hfulladmix)
}


hyb.contam.admix <- hyb.admix.plot("contam")
hyb.compmap.admix <- hyb.admix.plot("compmap")
hyb.adr.admix <- hyb.admix.plot("adr")
hyb.haploid.admix <- hyb.admix.plot("haploid")


setEPS()
postscript("./figures/SUPP4_hybadmix.eps", width = 25, height = 18)


print((hyb.contam.admix + labs(tag = "A")) |> insert_right(((hyb.compmap.admix + labs(tag = "B")))) |>
  insert_right(((hyb.haploid.admix + labs(tag = "C")))) |>
  insert_right(((hyb.adr.admix + labs(tag = "D"))))) & theme(legend.position = "bottom")
dev.off()


# ingorup mapping####


allprop <- read.csv("data/all_mapprops.csv", header = F)
colnames(allprop) <- c("ID", "pre_mapped", "pre_unmapped", "post_mapped")


# Pre filtering #####
rm(Lasi_)
rm(`Lasi Chtono`)
rm(`Lasi jaune`)
rm(Lasi.Chtono)
rm(Lasi.jaune)
rm(Lasi_alie)
rm(`Lasi_alie gr`)
rm(Lasi_alie.gr)
rm(Lasi_brun)
rm(Lasi_emar)
rm(Lasi_flav)
rm(Lasi_fuli)
rm(Lasi_mixt)
rm(Lasi_myop)
rm(Lasi_nige)
rm(Lasi_umbr)
rm(Lasi_nige.plat)
rm(`Lasi_nige/plat`)
rm(Lasi_para)
rm(Lasi_plat)
rm(Lasi_psam)
ants <- read.csv("data/antmetadata.csv")
colnames(ants)[1] <- "CATALOGUENUMBER"

id.list <- read.table("data/all_individuals.txt")
lasius <- subset(ants, grepl("Las", SPECIESSUMMARY))

id.las <- lasius[lasius[, "CATALOGUENUMBER"] %in% unlist(id.list$V1), ]
id.las <- as.data.frame(cbind(id.las$CATALOGUENUMBER, id.las$SPECIESSUMMARY))

ir <- split(id.las, id.las$V2)

for (name in (names(ir))) {
  list <- (ir[[name]]$V1)
  assign(name, list)
}

Lasi_plat <- c(Lasi_plat, "4457")
Lasi_flav <- Lasi_flav[Lasi_flav != "9990946"]

allprop <- allprop %>%
  mutate(OldSpecies = case_when(
    try((allprop$ID %in% Lasi) ~ "Unidentified Lasius"),
    try((allprop$ID %in% `Lasi Chtono`) ~ "Unidentified Lasius"),
    try((allprop$ID %in% `Lasi jaune`) ~ "Unidentified Lasius"),
    try((allprop$ID %in% Lasi_alie) ~ "L.alienus"),
    try((allprop$ID %in% `Lasi_alie gr`) ~ "Unidentified Lasius"),
    try((allprop$ID %in% Lasi_brun) ~ "L.brunneus"),
    try((allprop$ID %in% Lasi_emar) ~ "L.emarginatus"),
    try((allprop$ID %in% Lasi_flav) ~ "L.flavus"),
    try((allprop$ID %in% Lasi_fuli) ~ "L.fuliginosus"),
    try((allprop$ID %in% Lasi_mixt) ~ "L.mixtus"),
    try((allprop$ID %in% Lasi_myop) ~ "L.myops"),
    try((allprop$ID %in% Lasi_nige) ~ "L.niger"),
    try((allprop$ID %in% `Lasi_nige/plat`) ~ "Unidentified Lasius"),
    try((allprop$ID %in% Lasi_para) ~ "L.paralienus"),
    try((allprop$ID %in% Lasi_plat) ~ "L.platythorax"),
    try((allprop$ID %in% Lasi_psam) ~ "L.neglectus"),
    try((allprop$ID %in% Lasi_umbr) ~ "L.umbratus"),
    try((allprop$ID %in% c("2198", "4172B")) ~ "Unidentified Lasius")
  ))

allprop$preperc <- (100 * (allprop$pre_mapped)) / (allprop$pre_mapped + allprop$pre_unmapped)
allprop$postperc <- (100 * (allprop$post_mapped)) / (allprop$post_mapped + allprop$pre_unmapped)
#####
pre <- ggplot(data = allprop, aes(OldSpecies, preperc)) +
  geom_boxplot() +
  geom_jitter(width = 0.25, aes(size = (pre_mapped + pre_unmapped))) +
  labs(
    title = "Pre-filtering mapping to Lasius niger",
    x = "Species", y = "Percentage of reads mapping"
  ) +
  theme_pubr(x.text.angle = 90)


# Post filtering #####
rm(Lasi_)
rm(`Lasi Chtono`)
rm(`Lasi jaune`)
rm(Lasi.Chtono)
rm(Lasi.jaune)
rm(Lasi_alie)
rm(`Lasi_alie gr`)
rm(Lasi_alie.gr)
rm(Lasi_brun)
rm(Lasi_emar)
rm(Lasi_flav)
rm(Lasi_fuli)
rm(Lasi_mixt)
rm(Lasi_myop)
rm(Lasi_nige)
rm(Lasi_umbr)
rm(Lasi_nige.plat)
rm(`Lasi_nige/plat`)
rm(Lasi_para)
rm(Lasi_plat)
rm(Lasi_psam)

correctedID <- read.csv("data/correctID1.csv", header = T)
for (name in (colnames(correctedID))) {
  list <- na.exclude(correctedID[[name]])
  assign(name, list)
}
allprop <- allprop %>%
  mutate(NewSpecies = case_when(
    try((allprop$ID %in% Lasi) ~ "Unidentified Lasius"),
    try((allprop$ID %in% Lasi_alie) ~ "L.alienus"),
    try((allprop$ID %in% Lasi_brun) ~ "L.brunneus"),
    try((allprop$ID %in% Lasi_emar) ~ "L.emarginatus"),
    try((allprop$ID %in% Lasi_flav) ~ "L.flavus"),
    try((allprop$ID %in% Lasi_fuli) ~ "L.fuliginosus"),
    try((allprop$ID %in% Lasi_mixt) ~ "L.mixtus"),
    try((allprop$ID %in% Lasi_myop) ~ "L.myops"),
    try((allprop$ID %in% Lasi_nige) ~ "L.niger"),
    try((allprop$ID %in% Lasi_para) ~ "L.paralienus"),
    try((allprop$ID %in% Lasi_plat) ~ "L.platythorax"),
    try((allprop$ID %in% Lasi_psam) ~ "L.neglectus"),
    try((allprop$ID %in% Lasi_umbr) ~ "L.umbratus"),
    try((allprop$ID %in% c("18536", "2198", "3771", "4172B", "518", "7739B")) ~ "Unidentified Lasius")
  ))
#####
post<- ggplot(data = allprop, aes(NewSpecies, postperc)) +
  geom_boxplot() +
  geom_jitter(width = 0.25, aes(size = (post_mapped + pre_unmapped))) +
  labs(
    title = "Post-filtering mapping to Lasius niger",
    x = "Species", y = "Percentage of reads mapping"
  ) +
  theme_pubr(x.text.angle = 90)


setEPS()
postscript("./figures/SUPP5_ingroup.eps", width = 25, height = 18)

ggarrange(pre, post, ncol = 1)

dev.off()
