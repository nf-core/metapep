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

parser <- add_argument(parser, "--predictions", nargs=1, help="Input file containing: peptide_id, allele_id, prediction_score.")
parser <- add_argument(parser, "--proteins_peptides", nargs=1, help="Input file containing: protein_id, peptide_id, count.")
parser <- add_argument(parser, "--proteins_microbiomes", nargs=1, help="Input file containing: microbiome_id, protein_id, protein_weight.")
parser <- add_argument(parser, "--conditions", nargs=1, help="Input file containing: microbiome_id, condition_id, condition_name.")
parser <- add_argument(parser, "--conditions_alleles", nargs=1, help="Input file containing: condition_id, allele_id.")
parser <- add_argument(parser, "--alleles", nargs=1, help="Input file containing: allele_id, allele_name.")
#parser <- add_argument(parser, "--output", nargs=1, help="Output file name.")


args <- parse_args(parser)

predictions <- fread(args$predictions)
proteins_peptides <- fread(args$proteins_peptides)
proteins_microbiomes <- fread(args$proteins_microbiomes)
conditions <- fread(args$conditions)
conditions_alleles <- fread(args$conditions_alleles)
alleles <- fread(args$alleles)


data <-
  conditions %>%
  inner_join(conditions_alleles) %>%
  inner_join(alleles) %>%
  inner_join(predictions) %>%
  inner_join(proteins_peptides) %>%
  inner_join(proteins_microbiomes) %>%
  group_by(condition_name, peptide_id, prediction_score, allele_name) %>% 
  summarise(weight_sum = sum(protein_weight))


data$condition_name <- as.factor(data$condition_name)

# one distribution including all alleles
p1 <- ggplot(data, aes(x=condition_name, y=prediction_score, fill=condition_name)) +
    ylab("Epitope prediction score") +
    xlab("Condition") +
    geom_violin() +
    scale_fill_brewer(palette="Dark2") +
    geom_boxplot(width=0.05, fill="white") +
    theme_classic()+
    theme(legend.position="none")

ggsave("prediction_score_distribution.pdf", height=5, width=5)

# individual distributions for alleles
p2 <- ggplot(data, aes(x=condition_name, y=prediction_score, fill=condition_name)) +
    facet_grid(. ~ allele_name) +
    ylab("Epitope prediction score") +
    xlab("Condition") +
    geom_violin() +
    scale_fill_brewer(palette="Dark2") +  # "Spectral") +
    geom_boxplot(width=0.05, fill="white") +
    theme_classic() +
    theme(legend.position="none")

nc <- ceiling(nrow(alleles)/3)
if (nc == 1){
  nr <- nrow(alleles) %% 3
} else {
  nr <- 3
}

ggsave("prediction_score_distribution.alleles.pdf", height=5*nc, width=5*nr)
