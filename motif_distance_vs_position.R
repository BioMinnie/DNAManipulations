# Import library
library(ggplot2)

# Read in the data
# motif, start position, end position, distance, region in genome
motif_reg = read.csv(file.choose(), header=T, sep='\t')

# Attach variable (column) names
attach(motif_reg)

# Check first few observations
head(motif_reg)

# Summary statistics
dist_sd <- sd(motif_reg$distance)      # Standard deviation
dist_mean <- mean(motif_reg$distance)  # Mean
cutoff <- dist_mean + (3 * dist_sd)    # Generate a cutoff at mean + 3x stdev
cutoff
hist(motif_reg$distance)               # Visualise a histogram of the distance integers

# Normal R scatterplot - no colour:
distance = motif_reg$distance
position = motif_reg$start
plot(position, distance)

# ggplot2 scatterplot in colour:
# Using my own colours
motif_plot <- ggplot(motif_reg, aes(start, distance, color=factor(region))) + # Colour according to region
  xlab("Position in genome (bp)") +                                            # X label
  ylab("Distance between motifs (bp)") +                                 # Y label
  geom_point() + scale_color_manual(values=c("#0033CC", "#BABABA", "#00CC00", "#FF0000", "#FF9933", "#9900FF", "#FF99CC",
                                             "#0099FF", "#EE03A8", "#003300", "#993300", "#006666", "#66CCFF", "#660066"), 
                                             # Added the colours for each region (genome is grey)
                                    name = "Region", # Name of legend
                                    breaks=c("genome", "GI-thrW", "Prophage_1", "Prophage_2", "Prophage_3", 
                                               "Prophage_4", "Prophage_5", "Prophage_6", "HPI", "Cryptic_Phage", 
                                               "Prophage_7", "GI-pheV", "GI-selC", "GI-LeuX"), 
                                    labels=c("Genome", "GI-ThrW", "Prophage 1", "Prophage 2", "Prophage 3", "Progphage 4", "Prophage 5",
                                             "Prophage 6", "HPI", "Cryptic Phage", "Prophage 7", "GI-PheV", "GI-SelC", "GI-LeuX")) + 
                                             # Updated names for each region in the legend
  geom_hline(yintercept=1119.55, linetype="dashed") + # Add in a dashed line representing mean + 3 x stdev
  theme_bw(base_size = 8, base_family = "Helvetica") + theme(panel.border = element_rect(color = "dark grey", size = .3), 
                                                              panel.grid.major.y = element_line(size = .2, color = "grey"), 
                                                              panel.grid.minor.y = element_line(size = .2, color = "#E8E7E7"),
                                                              panel.grid.major.x = element_line(size = .2, color = "#E8E7E7"), 
                                                              panel.grid.minor.x = element_blank(),
                                                              axis.line = element_line(size=.7, color = "black"),
                                                              text = element_text(size=14)) 
                                                              # White theme, black border, grid lines, text size etc
motif_plot

# Save the file as a PDF
# This can be altered to any image file type
pdf(file = "filename.pdf",
    width = 12, 
    height = 8, 
    pointsize = 12, 
    res = 600) 
# Specifying larger sizes (height and width) will yield a PDF with smaller text and more space between points
# 600 dpi for publication quality resolution
motif_plot # Print the plot
dev.off()  # Turn plot off
