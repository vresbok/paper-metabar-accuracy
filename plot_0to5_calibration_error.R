# Plot the standard deviation of the relative error in abundance
# estimates for lysate and homogenate Swedish IBA data for each target
# spike-in, using 0 to 5 of the other spikeins for calibration

# Libraries needed
library(ggplot2)
library(patchwork)

source("get_0to5_calibration_error.R")

D <- lysate_cal_error
E <- homogenate_cal_error

# Define plot function
make_plot <- function(D) {
    
    ggplot(data=D, aes(x=num_includes, y=stdev_rel)) +
            theme_minimal() +
#            theme(axis.line=element_line())
            scale_x_continuous(name="No. calibrator spike-ins", position="bottom") +
            geom_point() +
            scale_y_continuous(name="Relative error (std dev)", limits=c(1,4.5))
}

# Make plots

p1 <- make_plot(D)

p2 <- make_plot(E)

ggsave(file="Figures/Fig_0to5_calibration_error.jpg", height=3.5, width=7, plot = (p1 + p2) + plot_layout(axis_titles="collect_y",ncol=2) + plot_annotation(tag_levels="A"))

