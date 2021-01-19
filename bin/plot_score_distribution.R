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
library(argparser, quietly=TRUE)
library(stringr)


parser <- arg_parser("Description")
parser <- add_argument(parser, "--scores", nargs=1, help="Input file containing: prediction_score, condition_name, weight_sum.")
parser <- add_argument(parser, "--alleles", nargs=1, help="Input file containing: allele_id, allele_name.")
parser <- add_argument(parser, "--conditions", nargs=1, help="Input file containing: microbiome_id, condition_id, condition_name.")
parser <- add_argument(parser, "--allele_id", nargs=1, help="allele_id.")
parser <- add_argument(parser, "--method", nargs=1, help="Epitope prediction method used.")
args <- parse_args(parser)

data <- fread(args$scores)
alleles <- fread(args$alleles)

allele_name <- alleles[alleles$allele_id == args$allele_id, ]$allele_name
allele_str <- str_replace_all(allele_name, '\\*', '_')
allele_str <- str_replace_all(allele_str, '\\:', '_')

if (args$method == "syfpeithi"){
  score_threshold <- 0.50
} else {
  score_threshold <- 500
}

data$condition_name <- as.factor(data$condition_name)

p <- ggplot(data, aes(x=condition_name, y=prediction_score, weight = weight_sum, fill=condition_name)) +
    ylab("Epitope prediction score") +
    xlab("Condition") +
    ggtitle(allele_name) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_brewer(palette="Dark2") + 
    geom_hline(yintercept=score_threshold) +
    theme_classic() +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))

ggsave(paste0("prediction_score_distribution.", allele_str, ".pdf"), height=5, width=5)

