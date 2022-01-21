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
binding_rates <- args[1]                # Input file containing: condition_name, binding_rate, entity_weight.
args_alleles <- args[2]                      # Input file containing: allele_id, allele_namee
allele_id <- args[3]                    # allele_id

data <- fread(binding_rates)
alleles <- fread(args_alleles)
match <- (alleles$allele_id == allele_id)
allele_name <- alleles[match, ]$allele_name
allele_str <- str_replace_all(allele_name, '\\*', '_')
allele_str <- str_replace_all(allele_str, '\\:', '_')

data$condition_name <- as.factor(data$condition_name)
p <- ggplot(data, aes(x=condition_name, y=binding_rate, fill=condition_name)) +
    ylab("Entity-wise binding ratio") +
    xlab("Condition") +
    ggtitle(allele_name) +
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    scale_fill_brewer(palette="Dark2") +
    theme_classic() +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))

ggsave(paste0("entity_binding_ratios.with_points.", allele_str, ".pdf"), height=5, width=5)


p2 <- ggplot(data, aes(x=condition_name, y=binding_rate, fill=condition_name)) +
    ylab("Entity-wise binding ratio") +
    xlab("Condition") +
    ggtitle(allele_name) +
    geom_boxplot() +
    scale_fill_brewer(palette="Dark2") +
    theme_classic() +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))

ggsave(paste0("entity_binding_ratios.", allele_str, ".pdf"), height=5, width=5)

# TODO use weights for violin plot