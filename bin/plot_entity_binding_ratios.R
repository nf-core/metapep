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

library(optparse)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(ggpubr)

# Parse command line arguments
option_list <- list(
    make_option(c("-r", "--ratios"     ), type="character", default=NULL    , metavar="path"   , help="Input file containing the precomputed binding ratios: condition_name, binding_rate, entity_weight."),
    make_option(c("-a", "--alleles"    ), type="character", default=NULL    , metavar="path"   , help="Alleles input file: allele_id, allele_name."                                                       ),
    make_option(c("-i", "--allele_id"  ), type="integer"  , default=0       , metavar="integer", help="Allele ID."                                                                                        ),
    make_option(c("-d", "--hide_pvalue"), type="logical"  , default=FALSE   , metavar="boolean", help="Disable mean comparison and do not show p-values."                                                          )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$ratios)){
    print_help(opt_parser)
    stop("Please provide the input file containing the precomputed binding ratios.", call.=FALSE)
}
if (is.null(opt$alleles)){
    print_help(opt_parser)
    stop("Please provide the alleles input file.", call.=FALSE)
}
if (is.null(opt$allele_id)){
    print_help(opt_parser)
    stop("Please provide the allele id.", call.=FALSE)
}

# Process data
data <- fread(opt$ratios)
alleles <- fread(opt$alleles)
match <- (alleles$allele_id == opt$allele_id)
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
    {if(!opt$hide_pvalue) stat_compare_means()} +
    scale_fill_brewer(palette="Dark2") +
    theme_classic() +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))

ggsave(paste0("entity_binding_ratios.with_points.", allele_str, ".pdf"), height=5, width=5)


p2 <- ggplot(data, aes(x=condition_name, y=binding_rate, fill=condition_name)) +
    ylab("Entity-wise binding ratio") +
    xlab("Condition") +
    ggtitle(allele_name) +
    geom_boxplot() +
    {if(!opt$hide_pvalue) stat_compare_means()} +
    scale_fill_brewer(palette="Dark2") +
    theme_classic() +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))

ggsave(paste0("entity_binding_ratios.", allele_str, ".pdf"), height=5, width=5)

# TODO use weights for violin plot
