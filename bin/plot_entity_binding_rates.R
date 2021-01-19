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
parser <- add_argument(parser, "--binding_rates", nargs=1, help="Input file containing: condition_name, binding_rate, entity_weight.")
parser <- add_argument(parser, "--alleles", nargs=1, help="Input file containing: allele_id, allele_name.")
parser <- add_argument(parser, "--allele_id", nargs=1, help="allele_id.")
parser <- add_argument(parser, "--output", nargs=1, help="Output file name.")
args <- parse_args(parser)

data <- fread(args$scores)
alleles <- fread(args$alleles)
allele_name <- alleles[alleles$allele_id == args$allele_id, ]$allele_name

data$condition_name <- as.factor(data$condition_name)
p <- ggplot(data, aes(x=condition_name, y=binding_rate, fill=condition_name)) +
    ylab("Entity-wise binding rate") +
    xlab("Condition") +
    ggtitle(allele_name) +
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    scale_fill_brewer(palette="Dark2") + 
    theme_classic() +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))

# TODO use weights for violin plot

ggsave(args$output, height=5, width=5)
