#' @title get_ppi_for_molecules
#'
#' @description This function gets the PPI data from STRINGdb for the protein
#' molecules provided.
#'
#' @rdname get_ppi_for_molecules
#' @name get_ppi_for_molecules
#'
#' @details This function gets the PPI data from STRINGdb for the protein
#' molecules provided.
#'
#' @return This function returns a data frame of the PPI for the molecules.
#'
#' @param RP.protein A vector containg the receptor (RP) proteins.
#' @param KN.protein A vector containg the kinase (KN) proteins.
#' @param TF.protein A vector containg the transcription factor (TF) proteins.
#' @param species The species name, either "hsapiens" or "mmusculus".
#' @param score An interger value for the STRINGdb PPI score threshold cutoff.
#'   Default is 700.
#'
#' @importFrom STRINGdb STRINGdb
#'
#' @export
#'
#' @examples
#' ## Need two folder at working directory for downloading stringdb files for each species - stringdb_mouse, stringdb_human.
#' ## It takes some time to download the data, and then can reuse the downloaded data.
#' ## Here we will use RP.protein, KN.protein, TF.protein protein parameters. These data are automatically loaded with the package. You can modify these parameters.
#' ## And we will use the species as "mmusculus".
#'
#' ## Get PPI data for the protein molecules of species "mmusculus".
#' mm.ppi <- get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species = "mmusculus")
#' head(mm.ppi)
get_ppi_for_molecules <- function(RP.protein, KN.protein, TF.protein, species, score = 700) {
  ## get ppi interactions for molecules
  if (species == "mmusculus") {
    fs::dir_create("stringdb_mouse")
    # initiate the connection, id  10090 for mouse
    string_db_mouse <- STRINGdb$new(version = "10", species = 10090, score_threshold = 0, input_directory = "stringdb_mouse")
    # now combine all the protein
    all.protein <- unique(c(RP.protein, KN.protein, TF.protein))
    # make a data frame from all the protein
    all.protein.df <- data.frame("gene" = all.protein)
    # mapping gene names to string ids
    all.protein.mapped <- string_db_mouse$map(all.protein.df, "gene", takeFirst = T, removeUnmappedRows = TRUE)
    # get interactions information
    all.protein.mapped.interactions <- string_db_mouse$get_interactions(all.protein.mapped$STRING_id)
    # get only interactions and score
    all.protein.mapped.interactions.score <- all.protein.mapped.interactions[, c(1, 2, 16)]
  }
  else if (species == "hsapiens") {
    fs::dir_create("stringdb_human")
    # initiate the connection, id  9606 for human
    string_db_human <- STRINGdb$new(version = "10", species = 9606, score_threshold = 0, input_directory = "stringdb_human")
    # now combine all the protein and make uppercase
    all.protein <- toupper(unique(c(RP.protein, KN.protein, TF.protein)))
    # make a data frame from all the protein
    all.protein.df <- data.frame("gene" = all.protein)
    # mapping gene names to string ids
    all.protein.mapped <- string_db_human$map(all.protein.df, "gene", takeFirst = T, removeUnmappedRows = TRUE)
    # get interactions information
    all.protein.mapped.interactions <- string_db_human$get_interactions(all.protein.mapped$STRING_id)
    # get only interactions and score
    all.protein.mapped.interactions.score <- all.protein.mapped.interactions[, c(1, 2, 16)]
  }
  else {
    stop("Do not support other species at this moment.")
  }
  ##
  
  
  ## from STRING_id to gene name conversion
  all.factor.M <- all.protein.mapped.interactions.score
  all.factor.N <- all.protein.mapped
  all.factor.M[, 1] <- all.factor.N[match(all.factor.M$from, all.factor.N$STRING_id), 1]
  all.factor.M[, 2] <- all.factor.N[match(all.factor.M$to, all.factor.N$STRING_id), 1]
  all.factor.PPI <- all.factor.M
  ##
  
  
  ## get only the significant interactions, here by default combined score >= 700
  all.factor.PPI.significant <- all.factor.PPI[all.factor.PPI$combined_score >= score, ]
  ##
  
  
  
  ######### To get all interactions without considering the directions
  ######### Here, we will take the highest score value for duplicates
  ## 1st get the original interactions
  all.ppi.sig.1 <- all.factor.PPI.significant
  rownames(all.ppi.sig.1) <- NULL
  ##
  
  
  #####
  ## combine the neighboring factors to treat as a single vector - original order
  comb.ppi.1 <- list()
  for (i in 1:nrow(all.ppi.sig.1)) {
    comb.ppi.1[[i]] <- paste(all.ppi.sig.1[i, 1], all.ppi.sig.1[i, 2], sep = "*")
  }
  ##
  
  ## make the first df (original order) with the combined_score
  comb.ppi.1.df <- data.frame("interaction" = unlist(comb.ppi.1), "score" = all.ppi.sig.1$combined_score)
  ##
  #####
  
  
  #####
  ## combine the neighboring factors to treat as a single vector - reverse order
  comb.ppi.2 <- list()
  for (j in 1:nrow(all.ppi.sig.1)) {
    comb.ppi.2[[j]] <- paste(all.ppi.sig.1[j, 2], all.ppi.sig.1[j, 1], sep = "*")
  }
  ##
  
  ## make the second df (reverse order) with the combined_score
  comb.ppi.2.df <- data.frame("interaction" = unlist(comb.ppi.2), "score" = all.ppi.sig.1$combined_score)
  ##
  #####
  
  
  ## Now add both the interactions' data frame - original order and reverse order
  comb.ppi.df <- rbind(comb.ppi.1.df, comb.ppi.2.df)
  ##
  
  ## order according to the score value - highest to lowest
  comb.ppi.df.ordered <- comb.ppi.df[order(comb.ppi.df$score, decreasing = T), ]
  ##
  
  ## take PPIs with the highest score valued unique one from the duplicates
  comb.ppi.df.ordered.unique <- comb.ppi.df.ordered[!duplicated(comb.ppi.df.ordered$interaction), ]
  rownames(comb.ppi.df.ordered.unique) <- NULL
  ##
  
  ## Finally return the combined PPI data frame with score value
  return(comb.ppi.df.ordered.unique)
  ##
  ##########
}

#' @title combine_mm_hs_ppi
#'
#' @description This function combines the PPI data for both "mmusculus" and
#'   "hsapiens" species created by the get_ppi_for_molecules function.
#'
#' @rdname combine_mm_hs_ppi
#' @name combine_mm_hs_ppi
#'
#' @details This function combines the PPI data for both "mmusculus" and
#'   "hsapiens" species created by the get_ppi_for_molecules function.
#'
#' @return This function returns a list consisting of the combined filtered PPI
#'   data, the RP proteins and the TF proteins of the combined filtered PPI data
#'   (RP-RP-KN-...-KN-TF and their list of RPs and TFs), to generate the pathway
#'   path.
#'
#' @param mm.ppi The PPI data for "mmusculus" species generated by the function
#'   get_ppi_for_molecules.
#' @param hs.ppi The PPI data for "hsapiens" species generated by the function
#'   get_ppi_for_molecules.
#' @param RP.protein A vector containg the same receptor (RP) proteins that are
#'   used in the function get_ppi_for_molecules.
#' @param KN.protein A vector containg the same kinase (KN) proteins that are
#'   used in the function get_ppi_for_molecules.
#' @param TF.protein A vector containg the same transcription factor (TF)
#'   proteins that are used in the function get_ppi_for_molecules.
#'
#' @export
#'
#' @examples
#' ## Need two folder at working directory for downloading stringdb files for each species - stringdb_mouse, stringdb_human.
#' ## It takes some time to download the data, and then can reuse the downloaded data.
#' ## Here we will use RP.protein, KN.protein, TF.protein protein parameters. These data are automatically loaded with the package. You can modify these parameters.
#' ## We will generate PPI data for two species - "mmusculus" and "hsapiens" by calling the function get_ppi_for_molecules two times.
#' ## Then we will combine these two PPI data sets by using the combine_mm_hs_ppi function that will be used later on to generate the pathway path data.
#'
#' ## Get PPI data for the protein molecules of species "mmusculus".
#' mm.ppi <- get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species = "mmusculus")
#' ## Get PPI data for the protein molecules of species "hsapiens".
#' hs.ppi <- get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species = "hsapiens")
#' ## Now combine and get the filtered PPI and the RP and TF proteins of the combined filtered PPI
#' comb.ppi.result <- combine_mm_hs_ppi(mm.ppi, hs.ppi, RP.protein, KN.protein, TF.protein)
#' head(summary(comb.ppi.result))
combine_mm_hs_ppi <- function(mm.ppi, hs.ppi, RP.protein, KN.protein, TF.protein) {
  ##### combine, order and take the unique PPIs with highest score
  # Combine the both PPI data
  comb.ppi <- rbind(mm.ppi, hs.ppi)
  # order according to the score value - highest to lowest
  comb.ppi.ordered <- comb.ppi[order(comb.ppi$score, decreasing = T), ]
  # take PPIs with the highest score valued unique one from the duplicates
  comb.ppi.ordered.unique <- comb.ppi.ordered[!duplicated(comb.ppi.ordered$interaction), ]
  #####
  
  #####
  # now separating the links using a list that contains all the links as vectors
  comb.ppi.interaction.split <- lapply(as.vector(comb.ppi.ordered.unique$interaction), function(x) {
    return(unlist(strsplit(x, split = "[*]")))
  })
  # making data frame from the unique split lists
  comb.ppi.interaction.split.df <- as.data.frame(do.call(rbind, lapply(comb.ppi.interaction.split, rbind)))
  # set the column names of the data frame
  colnames(comb.ppi.interaction.split.df) <- c("from", "to")
  # now add the score value as a 3rd column
  all.factor.PPI.significant <- data.frame(comb.ppi.interaction.split.df, "score" = as.vector(comb.ppi.ordered.unique$score))
  #####
  
  ##### First make all protein symbols as uppercase
  RP.protein <- toupper(RP.protein)
  KN.protein <- toupper(KN.protein)
  TF.protein <- toupper(TF.protein)
  #####
  
  
  ##### To get only the significant links exist from RP - KN - TF
  ##### FOr RP - RP, we have allowed maximum of 2 layers according to our design,
  ##### If you need different design you should change in this section according to your design.
  ## get interactions from RP to KN
  RP.to.KN.significant.ppi <- all.factor.PPI.significant[((all.factor.PPI.significant$from %in% RP.protein) &
                                                            (all.factor.PPI.significant$to %in% KN.protein)), ]
  
  ## get interactions from KN to KN - for all KNs
  KN.to.KN.significant.ppi <- all.factor.PPI.significant[((all.factor.PPI.significant$from %in% KN.protein) &
                                                            (all.factor.PPI.significant$to %in% KN.protein)), ]
  
  ## get interactions from KN to TF - for all KNs
  KN.to.TF.significant.ppi <- all.factor.PPI.significant[((all.factor.PPI.significant$from %in% KN.protein) &
                                                            (all.factor.PPI.significant$to %in% TF.protein)), ]
  
  ## get the RPs that have no direct interaction with the KNs
  RP.not.connected.with.KN <- setdiff(RP.protein, unique(RP.to.KN.significant.ppi$from))
  
  ## get the ppi from 'RP.not.connected.with.KN' to 'unique(RP.to.KN.significant.ppi$from)'
  # this will give us interaction for 2 RP layers
  # get interactions from from RP not connected with KN to RP connected with KN
  # these combined RPs will act as source to finding the paths
  RP.to.RP.significant.ppi <- all.factor.PPI.significant[((all.factor.PPI.significant$from %in% RP.not.connected.with.KN) &
                                                            (all.factor.PPI.significant$to %in% unique(RP.to.KN.significant.ppi$from))), ]
  
  ## And finally combine all the interactions from RP-RP-KN-...-KN-TF
  all.significant.filtered.ppi <- rbind(
    RP.to.RP.significant.ppi, RP.to.KN.significant.ppi,
    KN.to.KN.significant.ppi, KN.to.TF.significant.ppi
  )
  rownames(all.significant.filtered.ppi) <- NULL

  ## Now get the RP and TF of the interactions
  RPs <- unique(c(
    unique(as.vector(RP.to.RP.significant.ppi$from)),
    unique(as.vector(RP.to.KN.significant.ppi$from))
  ))
  
  TFs <- unique(as.vector(KN.to.TF.significant.ppi$to))

  ## Finally make a list of all.significant.filtered.ppi, RPs and TFs and then return
  comb.ppi.result <- list()
  comb.ppi.result[["PPI"]] <- all.significant.filtered.ppi
  comb.ppi.result[["RPs"]] <- RPs
  comb.ppi.result[["TFs"]] <- TFs
  return(comb.ppi.result)
}
