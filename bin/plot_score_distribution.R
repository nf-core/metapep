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


parser <- arg_parser("Description")
parser <- add_argument(parser, "--scores", nargs=1, help="Input file containing: prediction_score, microbiome_id, weight_sum.")
parser <- add_argument(parser, "--alleles", nargs=1, help="Input file containing: allele_id, allele_name.")
parser <- add_argument(parser, "--conditions", nargs=1, help="Input file containing: microbiome_id, condition_id, condition_name.")
parser <- add_argument(parser, "--allele_id", nargs=1, help="allele_id.")
parser <- add_argument(parser, "--method", nargs=1, help="Epitope prediction method used.")
parser <- add_argument(parser, "--output", nargs=1, help="Output file name.")
args <- parse_args(parser)

data <- fread(args$scores)
alleles <- fread(args$alleles)
conditions <- fread(args$conditions)    # TODO use condition name for plot

allele_name <- alleles[alleles$allele_id == args$allele_id, ]$allele_name

if (args$method == "syfpeithi"){
  score_threshold <- 0.50
} else {
  score_threshold <- 500
}

data$microbiome_id <- as.factor(data$microbiome_id)

p <- ggplot(data, aes(x=microbiome_id, y=prediction_score, weight = weight_sum, fill=microbiome_id)) +
    ylab("Epitope prediction score") +
    xlab("Microbiome ID") +
    ggtitle(allele_name) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_brewer(palette="Dark2") + 
    geom_hline(yintercept=score_threshold) +
    theme_classic() +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))

ggsave(args$output, height=5, width=5)
