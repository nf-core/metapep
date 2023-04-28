#!/usr/bin/env Rscript
####################################################################################################
#
# Author: Sabrina Krakau
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
####################################################################################################

library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
scores <- args[1]                               # Input file containing: prediction_score, condition_name, weight_sum
args_alleles <- args[2]                         # Input file containing: allele_id, allele_name
conditions <- args[3]                           # Input file containing: microbiome_id, condition_id, condition_name
allele_id <- args[4]                            # allele_id
method <- args[5]                               # Epitope prediction method used
syfpeithi_score_threshold <- args[6]            # Syfpeithi score threshold
mhcflurry_mhcnuggets_score_threshold <- args[7] # MHCflurry and MHCnuggets score thresholds

data <- fread(scores)
alleles <- fread(args_alleles)
match <- (alleles$allele_id == allele_id)
allele_name <- alleles[match, ]$allele_name
allele_str <- str_replace_all(allele_name, '\\*', '_')
allele_str <- str_replace_all(allele_str, '\\:', '_')

# Keep only rows with weight > 0
data <- data[data$weight_sum>0, ]

if (method == "syfpeithi"){
    score_threshold <- syfpeithi_score_threshold
} else {
    score_threshold <- mhcflurry_mhcnuggets_score_threshold
}

data$condition_name <- as.factor(data$condition_name)

p <- ggplot(data, aes(x=condition_name, y=prediction_score, weight = weight_sum, fill=condition_name)) +
    ylab("Epitope prediction score") +
    xlab("Condition") +
    ggtitle(allele_name) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), adjust = 0.2) +
    scale_fill_brewer(palette="Dark2") +
    geom_hline(yintercept=score_threshold) +
    theme_classic() +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))
# TODO which bandwidth to use best for smoothing?
# (somehow also affects printed quantiles)

ggsave(paste0("prediction_score_distribution.", allele_str, ".pdf"), height=5, width=5)

