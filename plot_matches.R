# Plot matches of lysate and homogenate reads

library(reshape2)
library(ggplot2)
library(patchwork)

# define plotting functions
plot_lys <- function(D,lbls) {
    ggplot(data=D, aes(x=cluster, y=log10(lysate_reads+1))) +
        theme_minimal() +
        geom_jitter(color="black",size=0.4,alpha=0.5) +
#        geom_dotplot(binaxis="y", stackdir="center", dotsize=0.4, alpha=0.5, binwidth=0.1) +
        geom_violin(alpha=0.5) +
        ylim(c(-0.05,6.5)) +
        scale_x_discrete(labels=lbls)
}

plot_hom <- function(D,lbls) {
    ggplot(data=D, aes(x=cluster, y=log10(homogenate_reads+1))) +
        theme_minimal() +
        geom_jitter(color="black",size=0.4,alpha=0.5) +
#        geom_dotplot(binaxis="y", stackdir="center", dotsize=0.4, alpha=0.5, binwidth=0.1) +
        geom_violin(alpha=0.5) +
        ylim(c(-0.05,6.5)) +
        scale_x_discrete(labels=lbls)
}

# Get matched data
D <- read.delim("matched_lysate_homogenate.tsv")
D <- D[D$lysate_reads >0 | D$homogenate_reads > 0,]

# Keep only the 20 most frequent clusters
x <- sort(table(D$cluster),decreasing=TRUE)
E <- D[D$cluster %in% names(x[1:20]),]
E$cluster <- factor(E$cluster,levels=names(x[1:20]))

# Define x axis labels
x_labels <- substr(names(x[1:20]),1,3)

p1 <- plot_lys(E, x_labels)
p2 <- plot_hom(E, x_labels)
ggsave(file="Fig_matched_most_occurring.jpg", height=7, width=14, plot = p1 + p2 + plot_layout(ncol=1) + plot_annotation(tag_levels="A"))

# Keep the clusters occurring in more than 100 samples
F <- D[D$cluster %in% names(x[x>100]),]
F$cluster <- factor(F$cluster,levels=names(x[x>100]))

# Aggregate data for clusters so we can find meaningful subgroups

# Mean reads is larger for lysates than homogenates; find out by how much
lys_sample_mean <- sum(D$lysate_reads)/length(unique(D$sample))
hom_sample_mean <- sum(D$homogenate_reads)/length(unique(D$sample))

L <- aggregate(F$lysate_reads[F$lysate_reads!=0]~F$cluster[F$lysate_reads!=0],FUN=mean)
H <- aggregate(F$homogenate_reads[F$homogenate_reads!=0]~F$cluster[F$homogenate_reads!=0],FUN=mean)
#L <- aggregate(F$lysate_reads~F$cluster,FUN=mean)
#H <- aggregate(F$homogenate_reads~F$cluster,FUN=mean)
colnames(L) <- c("cluster","mean_lys_reads")
colnames(H) <- c("cluster","mean_hom_reads")
LH <- merge(L, H, by="cluster")

set1 <- LH$cluster[order(LH$mean_lys_reads/lys_sample_mean+LH$mean_hom_reads/hom_sample_mean, decreasing=TRUE)]
set2 <- LH$cluster[order(LH$mean_lys_reads/LH$mean_hom_reads,decreasing=TRUE)]
set3 <- LH$cluster[order(LH$mean_lys_reads/LH$mean_hom_reads)]
set4 <- LH$cluster[order(LH$mean_lys_reads/lys_sample_mean+LH$mean_hom_reads/hom_sample_mean)]

cutoff <- 20
set1 <- set1[1:cutoff]
set2 <- set2[1:cutoff]
set3 <- set3[1:cutoff]
set4 <- set4[1:cutoff]

lbl_len <- 5
x_labels <- substr(set1,1,lbl_len)
p3 <- plot_lys(D[D$cluster %in% set1,], x_labels)
p4 <- plot_hom(D[D$cluster %in% set1,], x_labels)
ggsave(file="Fig_matched_most_reads.jpg", height=7, width=14, plot = p3 + p4 + plot_layout(ncol=1) + plot_annotation(tag_levels="A"))

x_labels <- substr(set2,1,lbl_len)
p5 <- plot_lys(D[D$cluster %in% set2,], x_labels)
p6 <- plot_hom(D[D$cluster %in% set2,], x_labels)
ggsave(file="Fig_matched_largest_lysate_surplus.jpg", height=7, width=14, plot = p5 + p6 + plot_layout(ncol=1) + plot_annotation(tag_levels="A"))

x_labels <- substr(set3,1,lbl_len)
p7 <- plot_lys(D[D$cluster %in% set3,], x_labels)
p8 <- plot_hom(D[D$cluster %in% set3,], x_labels)
ggsave(file="Fig_matched_largest_homogenate_surplus.jpg", height=7, width=14, plot = p7 + p8 + plot_layout(ncol=1) + plot_annotation(tag_levels="A"))

x_labels <- substr(set4,1,lbl_len)
p9 <- plot_lys(D[D$cluster %in% set4,], x_labels)
p10 <- plot_hom(D[D$cluster %in% set4,], x_labels)
ggsave(file="Fig_matched_least_reads.jpg", height=7, width=14, plot = p9 + p10 + plot_layout(ncol=1) + plot_annotation(tag_levels="A"))

