###############################################################################
## AUTHOR: Melinda Ashcroft
## AFFILIATIONS: University of Queensland
## STRATEGY: Plot a bar plot of counts of variable per genomic segment.
## INPUT: Takes in a tab delimited file of segment number, count of variable 
## in that segment and genome position (core-genome or MGE).
###############################################################################

# Import required libraries
library(ggplot2)
library(scales)
library(RSvgDevice)

# Read in the data
motif_reg = read.csv(file.choose(), header=T, sep='\t')

# Attach variable (column) names
attach(motif_reg)

# Check first few observations
head(motif_reg)

# Summary statistics
dist_sd <- sd(motif_reg$distance)      # standard deviation
dist_mean <- mean(motif_reg$distance)  # mean
cutoff <- dist_mean + (3 * dist_sd)    # generate cut-off = mean + (3x stdev)
cutoff
hist(motif_reg$distance)               # visualise distribution of distance

# Normal scatterplot - no colour:
distance = motif_reg$distance
position = motif_reg$start
plot(position, distance)

# bar graph test
segment = motif_reg$Segment
GATC = motif_reg$GATC

# Legend is not correct - have manually edited it in Adobe Illustrator CS6
s.plot <- ggplot(data=motif_reg, aes(x=segment, y=GATC, color=factor(Region))) 
  + geom_bar(stat="identity") + 
  xlab("Segments in EC958 genome (1000bp, 250bp overlap)") +
  ylab("GATC Motif Counts") +
  scale_colour_manual(values=c("#BABABA", "#4169E1", "#DB3D95"), 
                      name = "Region",
                      breaks=c("Genome", "Genomic_Island", "Prophage"), 
                      labels=c("Genome", "Genomic Island", "Prophage")) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.border = element_rect(color = "dark grey", size = .3), 
        axis.line = element_line(size=.7, color = "black"),
        axis.text.x = element_text(size=12, angle=90),
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(vjust=0), 
        axis.title.y = element_text(vjust=1))

# Plot image
s.plot

# This section is used to force the axis to the limits of the data only
t.plot <- s.plot + scale_x_continuous(expand = c(0,0), 
breaks = c(0,1000,2000,3000,4000,5000,6000,7000)) + 
scale_y_continuous(expand = c(0,0))

# Plot updated image
t.plot

# Save the updated image in an svg file
devSVG(file="/Users/Mel/Desktop/GATC_per_segment_GI_Phi.svg", width=7, 
height=3.85) # Open plotting device
  t.plot     # Plot image
dev.off()    # Turn off plotting device
