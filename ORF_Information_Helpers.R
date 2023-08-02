# ORF_Information_Helpers.R

#----- Load Libraries -----#
library("Biostrings")
library("tidyverse")

#----- Defininitions -----#
setClass("Codon",
         representation = representation(sequence_vector = "character",
                                         protein = "character",
                                         gene_name = "character",
                                         protein_position = "character"),
         prototype = prototype(sequence_vector = NA_character_,
                               protein = NA_character_,
                               gene_name = NA_character_,
                               protein_position = NA_character_)
)

setClass("ORF",
         representation = representation(name = "character",
                                         range = "character",
                                         LUT = "character",
                                         codons = "list"),
         prototype = prototype(name = NA_character_,
                               range = NA_character_,
                               LUT = NA_character_))
