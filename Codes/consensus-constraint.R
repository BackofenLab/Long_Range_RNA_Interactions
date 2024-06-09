
suppressWarnings(
  suppressPackageStartupMessages({
    library(optparse)
    # tidyverse packages
    library(readr)
    library(magrittr)
    library(tidyr)
    library(dplyr)
    library(purrr)
    library(stringr)
  })
)

################################################
# alignment file
clustalw.alignment.file <- NA
# constraint file
locarna.constraint.file <- NA
# constraint type
constraint.type <- "S"
# minimal number of constraints to be considered in consensus
min.constraints <- 2
# minimal fraction of similar constraints to be considered in consensus
min.fraction <- 0.7
################################################


################################################
# check if script is run within rstudio
################################################
if (isNamespaceLoaded("rstudioapi")) {

  ################################################
  # set local files for testing
  ################################################

    # change directory to script file location
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  # set alignment file
  clustalw.alignment.file <- 
    "result.aln"
  
  # set constraint file
  locarna.constraint.file <- 
    # "in-cm.fa"
    "input-constraints.fa"

  # set constraint type  
  constraint.type="S"
  
} else {
  
  ################################################
  # parse command line arguments
  ################################################
  
  myArgs <-
    # create an option parser with general tool information message
    # produces automatically the "--help" and "-h" option
    OptionParser( description =
                    str_c(
                    "\n",
                    "Reads a LocARNA input FASTA file with structure (#S) or fixed\n",
                    "structure (#FS) constraints and a corresponding ClustalW alignment\n",
                    "file produced by LocARNA using the given constraints.\n",
                    "Individual sequence constraints are mapped to respective alignment positions\n",
                    "to generate a consensus constraint that can be used with RNAalifold\n",
                    "to predict a constraint consensus RNA secondary structure of the alignment.\n",
                    "\n",
                    "NOTE: requires a NESTED structure constraint FOR EACH aligned sequence!\n",
                    "\n",
                    "Supported constraint encodings are '()<>.x', see LocARNA documentation.\n",
                    "All other encodings are silently ignored."
                    ),
                  add_help_option = T
                  ) |> 
    # define file option to parse
    add_option( c("-a","--alignment"),
                type="character",
                metavar = "clustalw.file",
                help = "LocARNA Clustal-w alignment output file to map the constaints to") |> 
    add_option( c("-c","--constraint"),
                type="character",
                metavar = "fasta.file",
                help = str_c(
                  "LocARNA FASTA input file with structure constraints used to generate the alignment"
                )) |> 
    add_option( c("-t","--type"),
                type="character",
                metavar = "S|FS",
                default = constraint.type,
                help = str_c(
                  "LocARNA structure constraint type to be used: (S)tructure constraint or (FS) = fixed structure constraint.",
                  "Default: '",constraint.type,"'"
                )) |> 
    add_option( c("-m","--min"),
                type="integer",
                metavar = "(>0)",
                default = as.character(min.constraints),
                help = str_c(
                  "Minimal number of similar constraints per position to be considered in consensus.",
                  "Default: ", min.constraints
                )) |> 
    add_option( c("-f","--fraction"),
                type="double",
                metavar = "[0.6..1]",
                default = as.character(min.fraction),
                help = str_c(
                  "Minimal fraction of similar constraints per position to be considered in consensus.",
                  "Default: ", min.fraction
                ))  
    
  
  # parse and store command line arguments
  # prints help message and exits if "--help" or "-h" found among arguments
  args <- parse_args(myArgs, args = commandArgs(trailingOnly = T) )
  
  # store parsing results
  clustalw.alignment.file <- args$alignment
  locarna.constraint.file <- args$constraint
  constraint.type <- args$type
  min.constraints <- args$min
  min.fraction <- args$fraction
  
  # sanity check
  if (is.null(clustalw.alignment.file) | is.null(locarna.constraint.file)) {
    stop("Alignment and constraint file must be provided.")
  }
}

# input validation
if (min.constraints < 1) {
  stop("Minimal number of constraints to be considered in consensus must be greater than 0.")
}
if ( ! (constraint.type %in% c("S","FS")) ) {
  stop("Constraint type must be either 'S' or 'FS'")
}
if ( min.fraction>1 | min.fraction<0.6) {
  stop("Minimal fraction of similar constraints to be considered in consensus must be between 0.6 and 1.")
}
if ( ! file.exists(clustalw.alignment.file) ) {
  stop(str_c("Alignment file '",clustalw.alignment.file,"' does not exist."))
}
if ( ! file.exists(locarna.constraint.file) ) {
  stop(str_c("Constraint file '",locarna.constraint.file,"' does not exist."))
}

# read alignment 
alignment <- 
  read_tsv(clustalw.alignment.file, skip=1, col_names = F, show_col_types = F) |> 
  separate(1, into = c("genome","alignment"), sep = "\\s+") |> 
  filter( str_detect(genome,"^[^#]") ) |> 
  group_by(genome) |> 
  summarize( alignment = str_c(alignment, collapse=""))

# read constraints
constraints <-
  read_tsv(
    locarna.constraint.file,
           col_names="dat", show_col_types = F ) |> 
  filter( str_detect(dat,str_c("\\s#",constraint.type,"\\s*$")) | str_detect(dat,"^>") ) |> 
  # number ids and constraints consecutively (assuming one constraint FOR EACH sequence)
  mutate( isID = str_detect(dat,"^>") ) |>
  group_by( isID ) |>
  mutate( j = row_number() ) |> 
  ungroup() |>
  group_by(j) |> 
  pivot_wider( names_from = isID, values_from = 1, names_prefix = "i" ) |> 
  ungroup() |> 
  transmute( genome = str_sub(iTRUE,2,-1),
             constraint = str_remove(iFALSE, "\\s+\\S+\\s*$" ) )

if( inner_join(alignment, constraints, by="genome") |> nrow() != nrow(alignment)) {
  stop("Number of constraints does not match number of sequences in alignment or sequence IDs differ. Check input files.")
}

##################################################################
# update constraints to alignment length
##################################################################


# function that returns provides for each position in the input string the
# for positions with open brackets the count of opening brackets up to that position
# and for positions with closing brackets the count of closing brackets up to that position

get_open_bracket_count_per_position <- function(x) {
  x <- strsplit(x, "")[[1]]
  open_bracket_count <- 0
  close_bracket_count <- 0
  result <- numeric(length(x))
  for (i in seq_along(x)) {
    if (x[i] == "(") {
      open_bracket_count <- open_bracket_count + 1
      result[i] <- open_bracket_count
    } else if (x[i] == ")") {
      result[i] <- open_bracket_count - close_bracket_count
      close_bracket_count <- close_bracket_count +1
    } else if (x[i] == "<") {
      result[i] <- -1
    } else if (x[i] == ">") {
      result[i] <- -2
    } else if (x[i] == "x") {
      result[i] <- 0
    } else {
      result[i] <- NA
    }
  }
  if (open_bracket_count != close_bracket_count) {
    stop(str_c("Unbalanced opening brackets in input string. Open: ",open_bracket_count," Close: ",close_bracket_count), call. = F)
  }
  result
}

# function that provides for each vector element the index of the element 
# that shows the same value if the value is > 0
get_matching_bracket_index <- function(x) {
  result <- numeric(length(x))
  for (i in seq_along(x)) {
    if (!is.na(x[i]) & x[i] > 0) {
      result[i] <- which(x == x[i] & seq_along(x) != i)[1]
    } else {
      result[i] <- x[i]
    }
  }
  result
}

# con <- ".((<.x.))>"
# get_open_bracket_count_per_position(con)
# get_matching_bracket_index(get_open_bracket_count_per_position(con))


get_alignment_position <- function(x) {
  x <- strsplit(x, "")[[1]]
  result <- if_else(x == "-", NA, 1:length(x))
  result[!is.na(result)]
}

# constraint mapping including catching of error messages
map_constraint_to_alignment <- safely(function(aln, con) {
  alnPos <- get_alignment_position(aln)
  alnCon <- rep(NA,str_length(aln[1]))
  alnCon[alnPos] <- get_open_bracket_count_per_position(con)
  alnCon <- get_matching_bracket_index(alnCon)
  alnCon
})
# aln <- "A--AAAA-AAA-AA--"
# alnPos <- get_alignment_position(aln)
# 
# alnCon <- rep(NA,str_length(aln[1]))
# conPos <- get_matching_bracket_index(get_open_bracket_count_per_position(con))
# alnCon[alnPos] <- conPos


# map constraints to alignment
mapped_constraints <- 
  left_join( alignment, constraints, by="genome" ) |> 
  transmute( aln = alignment, con = constraint) |> 
  rowwise() |>
  pmap(.f = map_constraint_to_alignment) |> 
  set_names( alignment$genome ) |> 
  transpose() 

# handle mapping errors
if (mapped_constraints$error |> compact() |> length() > 0) {
  mapping_errors <-
    mapped_constraints$error |> 
    compact() |> 
    unlist() |> 
    bind_rows() |> 
    pivot_longer(everything(), names_to = "sequence", values_to = "error.message") |> 
    mutate( sequence = str_remove(sequence,".message$"))
  #print error message to stderr
  write(str_c("Errors of mapping constraints to alignment:\n",
                   transmute(mapping_errors, out = str_c("- ",sequence,": ",error.message)) |> 
                     pull(out) |> 
                     str_c(collapse = "\n")
              ),
        stderr())
  # stop execution
  stop("Please correct the errors and try again.")
}

# provides the maximal frequency of any value among the provided ones
getMaxFrequency <- function(n=NA, ...) {
  print(c(...))
  return( max( table(c(...)), rm.na=T ) / n )
}
#getMaxFrequency(1,3,23,4,234,234234,2,2,2,n=1)

# generate consensus constraint encoding
mapped_constraints$result |> 
  bind_cols( ) |> 
  mutate( pos = row_number() ) |>
  rowwise() |> 
  transmute(pos, all = list( c_across(-pos))) |> 
  mutate( 
    # get count of distinct constraint encodings at each position
    consCount = unlist(all) |> unique() |> length(),
    `c(` = sum(na.omit(unlist(all)) > pos ),
    # max number of equal downstream base pairing partner
    `c(max` = ifelse(`c(`==0,0, na.omit(unlist(all))[na.omit(unlist(all)) > pos] |> table() |> max()), 
    `c)` = sum(na.omit(unlist(all)) < pos & na.omit(unlist(all)) > 0 ),
    # max number of equal downstream base pairing partner
    `c)max` = ifelse(`c)`==0,0, na.omit(unlist(all))[na.omit(unlist(all)) < pos & na.omit(unlist(all)) > 0 ] |> table() |> max()), 
    `c<` = sum(na.omit(unlist(all)) == -1 ),
    `c>` = sum(na.omit(unlist(all)) == -2 ),
    `cx` = sum(na.omit(unlist(all)) == 0 ),
    fMax = max(c(`c(max`, `c)max`, `c>`, `c<` ))/length(all), # maximum frequency of any pair constraint
    cons = ifelse(consCount==1,
                  # agreement in constraint
                  case_when(
                    `c(`>0 ~ "(", # exactly the same pairing partner within alignment
                    `c)`>0 ~ ")", # exactly the same pairing partner within alignment
                    `c<`>0 ~ "<",
                    `c>`>0 ~ ">",
                    `cx`>0 ~ "x"
                  ),
                  ifelse( `cx`>0,
                          # at least one block constraint 'x' at this position
                          "x",
                          ifelse( fMax >= min.fraction,
                                  # maximally frequent constraint
                                  case_when(
                                    near(`c(max`/length(all),fMax) ~ "(", # mostly the same pairing partner within alignment
                                    near(`c)max`/length(all),fMax) ~ ")", # mostly the same pairing partner within alignment
                                    near(`c<`/length(all),fMax) ~ "<",
                                    near(`c>`/length(all),fMax) ~ ">"
                                  ),
                                  ifelse( sum(`c<`,`c(`)>=min.constraints & sum(`c>`,`c)`)==0,
                                          # only upstream pair constraints at this position
                                          "<",
                                          ifelse( sum(`c>`,`c)`)>=min.constraints & sum(`c<`,`c(`)==0,
                                                  # only downstream pair constraints at this position
                                                  ">",
                                                  # mixed pair constraints at this position
                                                  NA
                                          )))))
  ) |> 
  mutate( cons = ifelse( is.na(cons), ".", cons ) ) |>
  # view()
  pull(cons) |>
  str_c(collapse="") |> 
  writeLines(sep="") # write consensus constraint to console


