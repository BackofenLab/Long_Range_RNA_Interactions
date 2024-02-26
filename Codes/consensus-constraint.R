#####install.packages("optparse", repos = "http://cran.us.r-project.org")
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
                    "Reads a LocARNA input fasta file with structure constraints (#S)\n",
                    "and a corresponding ClustalW alignment LocARNA output file\n",
                    "and maps the constraints to the alignment positions to\n",
                    "generate a consensus constraint to be used with RNAalifold.\n",
                    "\n",
                    "NOTE: requires a structure constraint FOR EACH aligned sequence!"
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
                  "LocARNA FASTA input file with '#S' structure constraints used to generate the alignment"
                  )) 
  
  # parse and store command line arguments
  # prints help message and exits if "--help" or "-h" found among arguments
  args <- parse_args(myArgs, args = commandArgs(trailingOnly = T) )
  
  # store parsing results
  clustalw.alignment.file <- args$alignment
  locarna.constraint.file <- args$constraint
  
  if (is.na(clustalw.alignment.file) | is.na(locarna.constraint.file)) {
    stop("Alignment and constraint file must be provided.")
  }
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
  filter( str_detect(dat,"\\s#S\\s*$") | str_detect(dat,"^>") ) |> 
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
      open_bracket_count <- open_bracket_count + 1
      result[i] <- -open_bracket_count
    } else if (x[i] == ">") {
      # result[i] <- -2
      result[i] <- -(open_bracket_count - close_bracket_count)
      close_bracket_count <- close_bracket_count +1
    } else if (x[i] == "x") {
      result[i] <- 0
    } else {
      result[i] <- NA
    }
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

map_constraint_to_alignment <- function(aln, con) {
  alnPos <- get_alignment_position(aln)
  alnCon <- rep(NA,str_length(aln[1]))
  alnCon[alnPos] <- get_open_bracket_count_per_position(con)
  alnCon <- get_matching_bracket_index(alnCon)
  alnCon
}
# aln <- "A--AAAA-AAA-AA--"
# alnPos <- get_alignment_position(aln)
# 
# alnCon <- rep(NA,str_length(aln[1]))
# conPos <- get_matching_bracket_index(get_open_bracket_count_per_position(con))
# alnCon[alnPos] <- conPos


# checks whether a given symbol is within the constraint list
contains_constraint <- function( x, con=".") {
  sum(na.omit(unlist(x)) == con ) > 0
}
# contains_constraint( x = (c("A","B","C",NA)),"x")

left_join( alignment, constraints, by="genome" ) |> 
  transmute( aln = alignment, con = constraint) |> 
  rowwise() |>
  pmap(.f = map_constraint_to_alignment) |> 
  set_names( alignment$genome ) |>
  bind_cols( ) |> 
  mutate( pos = row_number() ) |>
  rowwise() |> 
  transmute(pos, all = list( c_across(-pos))) |> 
  mutate( 
    # get count of distinct constraint encodings at each position
    consCount = unlist(all) |> unique() |> length(),
    `c(` = sum(na.omit(unlist(all)) > pos ) > 0,
    `c)` = sum(na.omit(unlist(all)) < pos & na.omit(unlist(all)) > 0 ) > 0,
    `c<` = sum(na.omit(unlist(all)) == -1 ) > 0,
    `c>` = sum(na.omit(unlist(all)) == -2 ) > 0,
    `cx` = sum(na.omit(unlist(all)) == 0 ) > 0,
    cons = ifelse(consCount==1,
                  # agreement in constraint
                  case_when(
                    `c(` ~ "(", # exactly the same pairing partner within alignment
                    `c)` ~ ")", # exactly the same pairing partner within alignment
                    `c<` ~ "<",
                    `c>` ~ ">",
                    `cx` ~ "x"
                  ),
                  ifelse( `cx`,
                          # at least one block constraint 'x' at this position
                          "x",
                          ifelse( (`c<` | `c(`) & !(`c>` | `c)`),
                                  # only upstream pair constraints at this position
                                  "<",
                                  ifelse( (`c>` | `c)`) & !(`c<` | `c(`),
                                          # only downstream pair constraints at this position
                                          ">",
                                          # mixed pair constraints at this position
                                          NA)))) ) |>
  mutate( cons = ifelse( is.na(cons), ".", cons ) ) |>
  # view()
  pull(cons) |>
  str_c(collapse="") |> 
  writeLines() # write consensus constraint to console


