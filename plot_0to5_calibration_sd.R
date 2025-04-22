# Plot the log read sd for lysate and homogenate
# Swedish IBA data and 0 to 5 spike-in calibrations

# Libraries needed
library(reshape2)
library(ggplot2)
library(patchwork)

old_dir <- getwd()
setwd("~/dev/ms-repos-iba/utils/")
source("get_iba_co1_data_fxn.R")
setwd(old_dir)


# Define a function that generates all possible bool vectors
generate_includes <- function(len) {

    includes <- rep(FALSE,times=len)
    include_list <- list(includes)
    
    includes <- increment(includes)
    while (sum(includes)!=0) {
        include_list[[length(include_list)+1]] <- includes
        includes <- increment(includes)
    }

    return (include_list)
}

# Define a function that increments a bool vector
increment <- function(includes) {

    for (i in length(includes):1) {
        if (includes[i]==FALSE)
            break
    }
    if (includes[i]==TRUE)
        j <- i
    else {
        includes[i] <- TRUE
        j <- i+1
    }
    while (j<=length(includes)) {
        includes[j] <- FALSE
        j <- j + 1
    }
    return (includes)
}


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
metadata_path <- "~/dev/figshare-repos/iba/raw_data/v6/"
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

# List names of spikeins
bio_spikes <- c(
                "Blattidae_cluster1",
                "Gryllidae_cluster1",
                "Gryllidae_cluster2",
                "Drosophilidae_cluster1",
                "Drosophilidae_cluster2",
                "Drosophilidae_cluster3"
                )

# Define x axis labels
x_labels <- c("Sh.la","Gr.bi","Gr.su","Dr.bi","Dr.se","Dr.ja")


# Compute calibration sd for L counts with 0 to 5 spike-ins
# -------------------------------------------------------------------

# Set up data frames
D <- data.frame(list(
                    spikein=character(),
                    num_includes=numeric(),
                    Sh.la=logical(),
                    Gr.bi=logical(),
                    Gr.su=logical(),
                    Dr.bi=logical(),
                    Dr.se=logical(),
                    Dr.ja=logical(),
                    stddev=numeric()
                    )
                )
E <- D

# Iterate over biological spikeins
for (i in 1:length(bio_spikes)) {

    spikein <- bio_spikes[i]
    x <- log10(LSA[i,-1])

    # Get all possible bool vectors
    include_list <- generate_includes(6)

    for (j in 1:length(include_list)) {

        includes <- include_list[[j]]
        if (includes[i]==TRUE)
            next

        num_includes <- sum(includes)
        if (num_includes==0) {
            stddev <- sd(x)
        } else {
            y <- log10(colSums(LSA[includes,-1]))
            f <- y - mean(y)
            stddev <- sd(x - f)
        }
        D <- rbind(D,
                   list(spikein=spikein,
                        num_includes=num_includes,
                        Sh.la=includes[1],
                        Gr.bi=includes[2],
                        Gr.su=includes[3],
                        Dr.bi=includes[4],
                        Dr.se=includes[5],
                        Dr.ja=includes[6],
                        stddev=stddev
                        )
                   )
    }
}


# Compute calibration sd for H counts with 0 to 5 spike-ins
# -------------------------------------------------------------------

# Iterate over biological spikeins
for (i in 1:length(bio_spikes)) {

    spikein <- bio_spikes[i]
    x <- log10(HSAB[i,-1])

    # Get all possible bool vectors
    include_list <- generate_includes(6)

    for (j in 1:length(include_list)) {

        includes <- include_list[[j]]
        if (includes[i]==TRUE)
            next
 
        num_includes <- sum(includes)
        if (num_includes==0) {
            stddev <- sd(x)
        } else {
            y <- log10(colSums(HSAB[includes,-1]))
            f <- y - mean(y)
            stddev <- sd(x - f)
        }
        E <- rbind(E,
                   list(spikein=spikein,
                        num_includes=num_includes,
                        Sh.la=includes[1],
                        Gr.bi=includes[2],
                        Gr.su=includes[3],
                        Dr.bi=includes[4],
                        Dr.se=includes[5],
                        Dr.ja=includes[6],
                        stddev=stddev
                        )
                   )
    }
}


# Make plots
pdf("Fig_0to5_lysate_cal_stddev.pdf")
plot(D$stddev~D$num_includes,xlab="No. spike-ins",ylab="Std dev (log10)",main="Lysate calibration")
dev.off()

pdf("Fig_0to5_homogenate_cal_stddev.pdf")
plot(E$stddev~E$num_includes,xlab="No. spike-ins",ylab="Std dev (log10)",main="Homogenate calibration")
dev.off()


