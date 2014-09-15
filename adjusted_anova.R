#######################################################################################################
#                                                                                                     #
## AUTHOR: Melinda M Ashcroft                                                                         #
## AFFILIATION: Beatson Lab | SCMB - University of Queensland                                         #
#                                                                                                     #
## PURPOSE: Statistical evaluation of distribution of methylation motifs between the genome and MGEs  #
## Adjusted Analysis of Variance                                                                      #
## Date: 15 September 2014                                                                            #
#                                                                                                     #
## PUBLISHED: The following script is tailored to suit the data published in Forde et al. 2014        #
#                                                                                                     #
#######################################################################################################

# Read in libraries required
library(mvtnorm)
library(survival)
library(TH.data)
library(multcomp)
library(sandwich)

# Read in the data
methyl_counts = read.csv(file.choose(), header=T, sep='\t')

# Attach variable (column) names
attach(methyl_counts)

# Check first few observations
head(methyl_counts)

#### FOR GATC motifs ####
#########################

# To visualise the heteroscedasticity
# View a boxplot of each methylation motif count per segment, for each region
par(cex.axis=0.7)
boxplot(methyl_counts$GATC ~ methyl_counts$Region, las = 2)

# First, fit a common anova model and view the output
amod_GATC <- aov(GATC~Region, data=methyl_counts)
summary(amod_GATC)
# Then pass the fitted anova model to glht function, where hypothesis test = multiple comparisons of means
# vcov = vcovHC specifies use of heteroscedastic consistent covariance estimation HC3, accounting for 
# non-homogenous variances
# Use pairwise comparisons of just the genome to each MGE accounting for non-homogenous variances:
amod_GATC_2 <- glht(amod_GATC, mcp(Region=c('Genome - GI_ThrW = 0', 'Genome - Prophage_1 = 0',
                                       'Genome - Prophage_2 = 0', 'Genome - Prophage_3 = 0',
                                       'Genome - Prophage_4 = 0', 'Genome - Prophage_5 = 0', 
                                       'Genome - Prophage_6 = 0', 'Genome - HPI = 0', 
                                       'Genome - Cryptic_Phage = 0', 'Genome - Prophage_7 = 0', 
                                       'Genome - GI_PheV = 0', 'Genome - GI_SelC = 0',
                                       'Genome - GI_LeuX = 0')), vcov=vcovHC)

# Obtain adjusted p values assuring that the familywise error rate is not larger than the alpha significance region
# for all pairwise comparisons between groups
summary(amod_GATC_2)

# Compute the confidence levels for each difference of means using the pairwise comparisons of the genome to each MGE
confint(amod_GATC_2)
# Visualise the confidence intervals
par(mar=c(5,7,4,2)) 
plot(confint(amod_GATC_2), cex.axis=0.5)

#### FOR GAGACC motifs ####
###########################

# To visualise the heteroscedasticity
# View a boxplot of each methylation motif count per segment, for each region
par(cex.axis=0.7)
boxplot(methyl_counts$GAGACC ~ methyl_counts$Region, las = 2)

# First, fit a common anova model
amod_GAGACC <- aov(GAGACC~Region, data=methyl_counts)
summary(amod_GAGACC)
# Then pass the fitted anova model to glht function, where hypothesis test = multiple comparisons of means
# vcov = vcovHC specifies use of heteroscedastic consistent covariance estimation HC3, accounting for 
# non-homogenous variances
# Use pairwise comparisons of just the genome to each MGE accounting for non-homogenous variances:
amod_GAGACC_2 <- glht(amod_GAGACC, mcp(Region=c('Genome - GI_ThrW = 0', 'Genome - Prophage_1 = 0',
                                            'Genome - Prophage_2 = 0', 'Genome - Prophage_3 = 0',
                                            'Genome - Prophage_4 = 0', 'Genome - Prophage_5 = 0', 
                                            'Genome - Prophage_6 = 0', 'Genome - HPI = 0', 
                                            'Genome - Cryptic_Phage = 0', 'Genome - Prophage_7 = 0', 
                                            'Genome - GI_PheV = 0', 'Genome - GI_SelC = 0',
                                            'Genome - GI_LeuX = 0')), vcov=vcovHC)
summary(amod_GAGACC_2)

# Compute the confidence levels for each difference of means
confint(amod_GAGACC_2)
# Visualise the confidence intervals
plot(confint(amod_GAGACC_2), cex.axis=0.5)

#### FOR RTACNNNNGTG motifs ####
################################

# To visualise the heteroscedasticity
# View a boxplot of each methylation motif count per segment, for each region
par(cex.axis=0.7)
boxplot(methyl_counts$RTACNNNNGTG ~ methyl_counts$Region, las = 2)

# First, fit a common anova model
amod_RTACNNNNGTG <- aov(RTACNNNNGTG~Region, data=methyl_counts)
summary(amod_RTACNNNNGTG)
# Then pass the fitted anova model to glht function, where hypothesis test = multiple comparisons of means
# vcov = vcovHC specifies use of heteroscedastic consistent covariance estimation HC3, accounting for 
# non-homogenous variances
# Use pairwise comparisons of just the genome to each MGE accounting for non-homogenous variances:
amod_RTACNNNNGTG_2 <- glht(amod_RTACNNNNGTG, mcp(Region=c('Genome - GI_ThrW = 0', 'Genome - Prophage_1 = 0',
                                            'Genome - Prophage_2 = 0', 'Genome - Prophage_3 = 0',
                                            'Genome - Prophage_4 = 0', 'Genome - Prophage_5 = 0', 
                                            'Genome - Prophage_6 = 0', 'Genome - HPI = 0', 
                                            'Genome - Cryptic_Phage = 0', 'Genome - Prophage_7 = 0', 
                                            'Genome - GI_PheV = 0', 'Genome - GI_SelC = 0',
                                            'Genome - GI_LeuX = 0')), vcov=vcovHC)
summary(amod_RTACNNNNGTG_2)

# Compute the confidence levels for each difference of means
confint(amod_RTACNNNNGTG_2)
# Visualise the confidence intervals
plot(confint(amod_RTACNNNNGTG_2), cex.axis=0.5)

#### FOR CANCATC motifs ####
############################

# To visualise the heteroscedasticity
# View a boxplot of each methylation motif count per segment, for each region
par(cex.axis=0.7)
boxplot(methyl_counts$CANCATC ~ methyl_counts$Region, las = 2)

# First, fit a common anova model
amod_CANCATC <- aov(CANCATC~Region, data=methyl_counts)
summary(amod_CANCATC)
# Then pass the fitted anova model to glht function, where hypothesis test = multiple comparisons of means
# vcov = vcovHC specifies use of heteroscedastic consistent covariance estimation HC3, accounting for 
# non-homogenous variances
# Use pairwise comparisons of just the genome to each MGE accounting for non-homogenous variances:
amod_CANCATC_2 <- glht(amod_CANCATC, mcp(Region=c('Genome - GI_ThrW = 0', 'Genome - Prophage_1 = 0',
                                                          'Genome - Prophage_2 = 0', 'Genome - Prophage_3 = 0',
                                                          'Genome - Prophage_4 = 0', 'Genome - Prophage_5 = 0', 
                                                          'Genome - Prophage_6 = 0', 'Genome - HPI = 0', 
                                                          'Genome - Cryptic_Phage = 0', 'Genome - Prophage_7 = 0', 
                                                          'Genome - GI_PheV = 0', 'Genome - GI_SelC = 0',
                                                          'Genome - GI_LeuX = 0')), vcov=vcovHC)
summary(amod_CANCATC_2)

# Compute the confidence levels for each difference of means
confint(amod_CANCATC_2)
# Visualise the confidence intervals
plot(confint(amod_CANCATC_2), cex.axis=0.5)

#### FOR AACNNNNCTTT motifs ####
################################

# To visualise the heteroscedasticity
# View a boxplot of each methylation motif count per segment, for each region
par(cex.axis=0.7)
boxplot(methyl_counts$AACNNNNCTTT ~ methyl_counts$Region, las = 2)

# First, fit a common anova model
amod_AACNNNNCTTT <- aov(AACNNNNCTTT~Region, data=methyl_counts)
summary(amod_AACNNNNCTTT)
# Then pass the fitted anova model to glht function, where hypothesis test = multiple comparisons of means
# vcov = vcovHC specifies use of heteroscedastic consistent covariance estimation HC3 accounting for 
# non-homogenous variances
# Use pairwise comparisons of just the genome to each MGE accounting for non-homogenous variances:
amod_AACNNNNCTTT_2 <- glht(amod_AACNNNNCTTT, mcp(Region=c('Genome - GI_ThrW = 0', 'Genome - Prophage_1 = 0',
                                                          'Genome - Prophage_2 = 0', 'Genome - Prophage_3 = 0',
                                                          'Genome - Prophage_4 = 0', 'Genome - Prophage_5 = 0', 
                                                          'Genome - Prophage_6 = 0', 'Genome - HPI = 0', 
                                                          'Genome - Cryptic_Phage = 0', 'Genome - Prophage_7 = 0', 
                                                          'Genome - GI_PheV = 0', 'Genome - GI_SelC = 0',
                                                          'Genome - GI_LeuX = 0')), vcov=vcovHC)
summary(amod_AACNNNNCTTT_2)

# Compute the confidence levels for each difference of means
confint(amod_AACNNNNCTTT_2)
# Visualise the confidence intervals
plot(confint(amod_AACNNNNCTTT_2), cex.axis=0.5)
