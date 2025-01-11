# Plot the lognormal read variance for various
# Swedish IBA data and spike-in calibrations

library(reshape2)
library(ggplot2)
library(patchwork)

old_dir <- getwd()
setwd("~/dev/ms-repos-iba/utils/")
source("get_iba_co1_data_fxn.R")
setwd(old_dir)


# Define plot function 
vio_plot <- function(D,x_labels) {
    ggplot(data=D, aes(x=cluster, y=log10(reads))) + 
        theme_minimal() +
        geom_violin() +
        ylim(c(-0.05,6.5)) +
        scale_x_discrete(labels=x_labels)
}


# Get homogenate and lysate data
data_path <- "~/dev/figshare-repos/iba/processed_data/v3/"
metadata_path <- "~/dev/figshare-repos/iba/raw_data/v4/"
L <- get_iba_co1_data(
                      data_path=data_path,
                      metadata_path=metadata_path,
                      country="SE",
                      dataset="lysate",
                      calibrate=FALSE,
                      remove_spikes=FALSE
                    )
H <- get_iba_co1_data(
                      data_path=data_path,
                      metadata_path=metadata_path,
                      country="SE",
                      dataset="homogenate",
                      calibrate=FALSE,
                      remove_spikes=FALSE
                     )

L <- data.frame(L)
H <- data.frame(H)


# Find sample columns
colnames(L) -> x
sample_cols_L <- which(grepl("P",x) & grepl("_",x))
colnames(H) -> x
sample_cols_H <- which(grepl("P",x) & grepl("_",x))

# Find spikeins
prop_samples_L <- rowMeans(L[,sample_cols_L]>0)
prop_samples_H <- rowMeans(H[,sample_cols_H]>0)

spikeins_L <- which(prop_samples_L > 0.95)
spikeins_H <- which(prop_samples_H > 0.95)

LS <- L[spikeins_L,c(which(colnames(L)=="cluster"),sample_cols_L)]
HS <- H[spikeins_H,c(which(colnames(H)=="cluster"),sample_cols_H)]

# Only keep columns containing all spikeins
LSA <- LS[,c(TRUE,colMeans(LS[,-1]>0)>0.95)]
HSA <- HS[,c(TRUE,colMeans(HS[,-1]>0)>0.95)]

# Make homogenate data version with only bio spikeins
HSAB <- HSA[!grepl("synth",HSA$cluster),]

# Make homogenate data version with only synthetic spikeins
HSAS <- HSA[grepl("synth",HSA$cluster),]

# List names of spikeins
bio_spikes <- c(
                "Blattidae_cluster1",
                "Gryllidae_cluster1",
                "Gryllidae_cluster2",
                "Drosophilidae_cluster1",
                "Drosophilidae_cluster2",
                "Drosophilidae_cluster3"
                )
synth_spikes <- c("Callio-synth","tp53-synth")

# Define x axis labels
x_labels <- c("Sh.la","Gr.bi","Gr.su","Dr.bi","Dr.se","Dr.ja")


# Violin plots for L counts (raw/calibrated with biological spike-ins)
# -------------------------------------------------------------------

# Raw values
D <- reshape2::melt(LSA,id="cluster",value.name="reads")

# Calibrate using geometric mean (leads to smaller variances than arithmetic mean)
# Leave-One-Out (LOO) calibration
LSA_LOO <- LSA[FALSE,]
for (i in 1:nrow(LSA)) {
    x <- log10(colSums(LSA[-i,-1]))
    f <- mean(x) - x
    LSA_LOO <- rbind(LSA_LOO,as.list(c(cluster=LSA$cluster[i],(10^f)*as.numeric(LSA[i,-1]))))
}
E <- reshape2::melt(LSA_LOO,id="cluster",value.name="reads")

D$cluster <- factor(D$cluster,levels=bio_spikes)
E$cluster <- factor(E$cluster,levels=bio_spikes)
D$reads <- as.numeric(D$reads)
E$reads <- as.numeric(E$reads)

p1 <- vio_plot(D, x_labels)
p2 <- vio_plot(E, x_labels)
ggsave(file="Fig_lysate_cals.jpg", height=3.5, width=7, plot = p1 + p2)


# Violin plots for H counts (raw/calibrated using syn/bio/both spike-ins) 
# -----------------------------------------------------------------------

# Raw values
D1 <- reshape2::melt(HSAB,id="cluster",value.name="reads")

# Calibrate using synthetic spike-ins
HSAB_S <- HSAB[FALSE,]
for (i in 1:nrow(HSAB)) {
    x <- log10(colSums(HSAS[,-1]))
    f <- mean(x) - x
    HSAB_S <- rbind(HSAB_S,as.list(c(cluster=HSAB$cluster[i],(10^f)*as.numeric(HSAB[i,-1]))))
}
D2 <- reshape2::melt(HSAB_S,id="cluster",value.name="reads")


# Calibrate using biological spike-ins
# Leave-One-Out (LOO) calibration
HSAB_LOO <- HSAB[FALSE,]
for (i in 1:nrow(HSAB)) {
    x <- log10(colSums(HSAB[-i,-1]))
    f <- mean(x) - x
    HSAB_LOO <- rbind(HSAB_LOO,as.list(c(cluster=HSAB$cluster[i],(10^f)*as.numeric(HSAB[i,-1]))))
}
D3 <- reshape2::melt(HSAB_LOO,id="cluster",value.name="reads")


# Calibrate using both synthetic and biological spike-ins,
# in a linear model (of log effects) fitted to observed
# read numbers. This is a theoretical maximum in variance
# reduction; an LOO exercise, analogous to deriving a calibration
# curve for (n-1) samples  before predicting the nth sample, is
# likely to yield similar results
HSAB_S_LOO_lm <- HSAB[FALSE,]
y <- log10(colSums(HSAS[,-1]))
z_y <- (y - mean(y)) / sd(y)
for (i in 1:nrow(HSAB)) {
    d <- log10(colSums(HSAB[-i,-1])) - y # This represents the biological spike-in info after dividing by the PCR factor evidenced by synthetic spike-ins (y)
    z_d <- (d - mean(d)) / sd(d)
    x <- log10(colSums(HSAB[i,-1]))
    z_x <- (x - mean(x)) / sd(x)
    fit <- lm(as.numeric(x) ~ as.numeric(y) + as.numeric(d))
#    print (fit)
#    print(summary(fit))
    f <- mean(x) + fit$residuals
    names(f) <- colnames(HSAB)[-1]
    HSAB_S_LOO_lm <- rbind(HSAB_S_LOO_lm,as.list(c(cluster=HSAB$cluster[i],10^f)))
}
D4 <- reshape2::melt(HSAB_S_LOO_lm,id="cluster",value.name="reads")


# Calibrate using first synthetic then biological spike-ins
# Here we combine the evidence from the synthetic spike-ins
# and biological spike-ins by converting the deviation to
# standard normal terms, and then scaling the deviation up
# to the same standard deviation as the two component
# distributions. This gives a slightly better result than
# just combining the two deviations, which is exactly the
# same thing as just using the biological spike-in information.
# The standard deviation of the two component distributions is
# fairly similar, so the synthetic spike-ins do not add a lot
# more information to that of the biological spike-ins
# It is possible to lower the sd of some species a little more
# by weighting the deviations from the two components differently
# but the method below lowers the sd of all but one species.
HSAB_S_LOO <- HSAB[FALSE,]
y <- log10(colSums(HSAS[,-1]))
z_y <- (y - mean(y)) / sd(y)
for (i in 1:nrow(HSAB)) {
    x <- log10(colSums(HSAB[-i,-1]))
    d <- x - y  # This represents the biological spike-in info after dividing by the PCR factor evidenced by synthetic spike-ins (y)
    z_d <- (d - mean(d)) / sd(d)
    z <- z_y + z_d  # Now we have desired balanced deviation for each sample
    f <- sqrt( (sd(d)^2 + sd(y)^2)/2 ) * z  # z-based deviation, lowering sd for 5 of 6 species
#    f <- 1.0*((x-y)-mean(x-y)) + 1.0*(y - mean(y)) # For comparison, this is the naive combination, which is identical to just using biological spike-ins
#    f <- ((2*sd(y)^2)/(sd(d)^2 + sd(y)^2))*(d-mean(d)) + ((2*sd(d)^2)/(sd(d)^2 + sd(y)^2))*(y-mean(y)) # Empirically found good combination (lowers sd for 3 of 6 spp)
#    f <- 0.85*((x-y)-mean(x-y)) + 1.0*(y - mean(y)) # Good combination
    HSAB_S_LOO <- rbind(HSAB_S_LOO,as.list(c(cluster=HSAB$cluster[i],HSAB[i,-1]/as.numeric(10^f))))
}
#D4 <- reshape2::melt(HSAB_S_LOO,id="cluster",value.name="reads")

D1$cluster <- factor(D1$cluster,levels=bio_spikes)
D2$cluster <- factor(D2$cluster,levels=bio_spikes)
D3$cluster <- factor(D3$cluster,levels=bio_spikes)
D4$cluster <- factor(D4$cluster,levels=bio_spikes)
D1$reads <- as.numeric(D1$reads)
D2$reads <- as.numeric(D2$reads)
D3$reads <- as.numeric(D3$reads)
D4$reads <- as.numeric(D4$reads)

# Uncomment this code to print standard deviations for all calibrations
#for (spike in bio_spikes) {
#    cat("Mean and stdevs (log scale) for",spike,"\n")
#    cat("raw:",mean(log10(D1$reads[D1$cluster==spike])),"--",sd(log10(D1$reads[D1$cluster==spike])),"\n")
#    cat("syn-calibrated:",mean(log10(D2$reads[D2$cluster==spike])),"--",sd(log10(D2$reads[D2$cluster==spike])),"\n")
#    cat("bio-calibrated:",mean(log10(D3$reads[D3$cluster==spike])),"--",sd(log10(D3$reads[D3$cluster==spike])),"\n")
#    cat("syn-bio-calibrated:",mean(log10(D4$reads[D4$cluster==spike])),"--",sd(log10(D4$reads[D4$cluster==spike])),"\n")
#}

p3 <- vio_plot(D1, x_labels)
p4 <- vio_plot(D2, x_labels)
p5 <- vio_plot(D3, x_labels)
p6 <- vio_plot(D4, x_labels)
ggsave(file="Fig_homogenate_cals.jpg", height=7, width=7, plot = (p3 + p4)/(p5 + p6))

# Join all calibrations in one figure
ggsave(file="Fig_calibrations.jpg", height=10.5, width=7, plot = (p1 + p2 + p3 + p4 + p5 + p6) + plot_layout(axis_titles="collect",ncol=2) + plot_annotation(tag_levels="A"))
