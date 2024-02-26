
library(tidyverse)
library(data.table)
library(ggpubr)
library(openxlsx)
library(colorspace)
library(dbscan)
library(factoextra)
library(ggrepel)
library(broom)
library(nlme)
library(ggVennDiagram)
library(svglite)
library(GGally)

## Functions-------------------------------------------

### The gene fitness is calculated as described in the following publication:
###
### Wetmore KM, Price MN, Waters RJ, Lamson JS, He J, Hoover CA, et al.
### Rapid quantification of mutant fitness in diverse bacteria by sequencing
### randomly bar-coded transposons. mBio. 2015;6(3):e00306-15.
### DOI 10.1128/mBio.00306-15


#' Calculate strain fitness
#'
#' @param n1 barcode abundance in experiment
#' @param n0 barcode abundance in control
#' @param psi pseudocount: sum of all barcode abundances in the experiment
#'            divided by sum of all barcode abundances in the control,
#'            multiplied by gene factor (see below)
#'
#' @return numeric value

strain_fitness <- function(n1, n0, psi){
  log2(n1 + sqrt(psi)) - log2(n0 + sqrt(1/psi))
}

#' Calculate strain variance based on Poisson noise
#'
#' @param n1 barcode abundance in experiment
#' @param n0 barcode abundance in control
#'
#' @return numeric value
#' 
strain_variance <- function(n1, n0){
  (1/(1+n1) + 1/(1+n0))/log(4)
}

# maximal allowed value of the variance
var_ceiling <- 1/strain_variance(20, 20)

#' Calculate the mode (most often value) in the distribution
#' 
#' @param distribution 
#'
#' @return numeric value
#' 
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Merges 2 strings together with underscore as a spacer
#'
#' @param x string
#' @param y string
#' @param ... additional arguments for paste
#'
#' @return string

paste_ <- function(x, y, ...){
  paste(x, y, sep = "_", ...)
}

#' Title
#'
#' @param df dataframe containing barcode counts data and sample information 
#' @param col name of column containing the experiment sample
#'
#' @return dataframe with calculated fitness values

analyze_condition <- function(df, col){
  ## calculate strain fitness and variance
  df[, prelimfit := strain_fitness(get(col), `control_0`, 1)]
  df[, varS := strain_variance(get(col), `control_0`)]

  ## generate weights for each strain
  df[, w_i := pmin(var_ceiling, 1/varS)]

  ## preliminary calculation of gene fitness (median of strain fitnesses) to calculate the gene factor
  df[, gene_fit_prelim := median(prelimfit), by = c("gene")][, gene_fit_prelim := gene_fit_prelim - median(gene_fit_prelim)]

  ## some genes are large and have many mutations, some genes are small and have just 1 ot 2 mutations.
  ## this discrepancy leads to outliers. to counter that, we use this value of gene factor
  df[, gene_factor := 2^(gene_fit_prelim)]
  df[, n := .N, by = c("gene")][n<=2, gene_factor:= 1]

  ## generate pseudocount values (total barcode ratio experiment/control, adjusted for strain)
  df[, psi := gene_factor*sum(get(col))/sum(`control_0`)]

  ## use the pseudocount data to get more accurate strain fitness values
  df[, fit := strain_fitness(get(col), `control_0`, psi)]

  ## calculate unnormalized gene fitness using the strain fitness values and the strain variances
  df[, gene_fitness_unnorm := weighted.mean(fit, w_i), by = c("gene")]

  ## generate a helper dataset for normalization, extracting only gene, position and gene fitness
  support <- df[, .(position, gene, gene_fitness_unnorm)]
  support <- support[, .(position = median(position)), .(gene, gene_fitness_unnorm)][order(position)]

  ## use running median to remove the local biases
  support[, bias := runmed((gene_fitness_unnorm), k = 251, endrule = "constant")]
  support[, geneFitness := gene_fitness_unnorm - bias]

  ## subtract the mode to center on 0
  M <- Mode(support$geneFitness)
  support[, geneFitness := geneFitness-M]

  ## merge the initial dataset to the helper dataset
  df[support, geneFitness := i.geneFitness,
     on = c("gene", "gene_fitness_unnorm")]

  ## remove the helper columns and rename the resulting one
  df[, c("prelimfit", "gene_factor", "n", "gene_fit_prelim", "psi", "fit", "varS", "w_i",
         "gene_fitness_unnorm") := NULL]
  setnames(df, old = c("geneFitness"), new = paste(col, c("geneFitness"), sep = "_"))
}

#' helper function to remove raw data from the dataset and bring it to the long format 
#'
#' @param x dataset
#' @param c vector with colnames to be removed
#'
#' @return dataframe

reshape_bigdata <- function(x, c){
  x[, c(c, "control_0") := NULL]
  x <- melt.data.table(x, measure.vars = patterns("F"))
  x[, c("rep", "sample", "day", "var") := tstrsplit(variable, "_", fixed=TRUE)]
  x[, variable:=NULL]
  x <- dcast(x, ... ~ var, value.var = "value")
}

#' Main function: calculate gene fitness for a biological replicate
#'
#' @param df dataframe containing barcode counts data and sample information
#' @param conds list of samples
#' @param threshold default value 3, to omit the barcodes/genes with too low abundance
#'
#' @return dataframe

just_do_it <- function(df, conds, threshold = 3){
  ## reshape dataset in the wide format
  ## each column is a sample, each row is a barcode
  ## fill NA with 1 to avoid division by 0
  a <- dcast(df, ... ~ rep+sample+day, value.var = "abundance", fill = 1)
  
  ## add gene information
  a <- merge(a, anno, by = "seq", all.x = T)
  
  ## QC: remove barcodes not corresponding to any genes
  a <- a[is.na(gene)==F]
  
  ## summarize the values in 3 control samples to improve detection
  a[, `control_0`:= sum(.SD), .SDcols = names(a) %like% "control", 1:nrow(a)]
  a[, c("1_control_0", "2_control_0", "3_control_0"):=NULL]
  
  ## QC: remove barcodes with abundance < threshold in control
  ## and genes with total abundance < 10*threshold
  a <- a[, genesum :=sum(`control_0`), "gene"]
  a <- a[genesum>threshold*10][, genesum:=NULL]
  a <- a[`control_0`>=threshold]
  
  setorder(a, position)

  ## calculate gene fitness values for each sample in the dataset
  ## adds gene fitness column for each condition
  map_df(conds, ~analyze_condition(a, .))

  ## remove raw data from the dataset and bring it to the long format
  a <- reshape_bigdata(a, conds)
  a <- unique(a[, .(gene, sample, day, rep, geneFitness)])
  a <- a[, .(geneFitness = mean(geneFitness)), .(gene, sample, day, rep)]
  a
}


#' combine results from 3 biological replicates to 1 dataframe
#'
#' @param df dataframe
#' @param ... 
#'
#' @return

byExp2full <- function(df, ...){
  df <- map(df, as.data.table)
  map2(df, 1:length(df), ~.x[, exp := rep(.y)])
  df <- rbindlist(df, ...)
}

## read annotation: barcode position
anno <- fread("anno_appended_upd.csv")
## read annotation: gene position
anno2 <- fread("gene_positions_W3110_NC.csv")

## generate 
anno <- merge(anno, anno2[, V1:=NULL], by = "gene")

anno[, pos_rel := (position - start)*100/(end-start)]
rm(anno2)

anno2 <- fread("D:/2020-10. Trying to bind it all/uniprot-proteome_UP000000625 (1).tab", na.strings = NULL)
setnames(anno2, old = c("Gene names  (primary )", "Gene names  (synonym )",
                        "Gene ontology (biological process)", "Gene ontology (cellular component)",
                        "Gene ontology (GO)", "Gene ontology (molecular function)"),
         new = c("gene", "Synonyms", "GO_bio", "GO_cell", "GO", "GO_mol"))
anno2 <- anno2[, .(gene, `Protein names`)]

### files

## read the file containing barcode counts
data <- fread("./seq/results.csv")

## sample list
conditions <- c("F_2", "NF_2", "F_4", "NF_4", "F_8", "NF_8")
conditions_wide <- outer(1:3, conditions, paste_) %>% as.vector()

data <- unique(data[, V1:=NULL])

## remove one bad quality sample
data <-  data[!(biol_repl==1&sample=="F"&day==4&rep==1)]
conditions_wide <- conditions_wide[conditions_wide!="1_F_4"]

data <- list(data[biol_repl==1],
             data[biol_repl==2],
             data[biol_repl==3])

## generate the dataset with the fitness values
output <- map(data, ~just_do_it(., conditions_wide, 3))
