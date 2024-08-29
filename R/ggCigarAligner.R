#' Align a Read to a Reference Genome
#'
#' This function aligns a specified read to its reference genome based on its 
#' CIGAR string. 
#' The user can select a read either by its query name (`qname`) or by its index
#' in the DataFrame. If multiple reads have the same `qname`, it is recommended 
#' to use the index parameter to specify the exact read to be processed.
#'
#' @param df_gg A S4 class DataFrame containing BAM/SAM alignment data, as obtained from 
#'   the `ggbamloader` function.
#' @param qname A character string specifying the query name (optional). If 
#'   provided, the function will attempt to select the read(s) with this 
#'   `qname`. If multiple reads have the same `qname`, the `index` parameter 
#'   should be used instead.
#' @param index An integer specifying the row index in the DataFrame (optional).
#'   If provided, this will be used to select the specific read regardless of 
#'   the `qname`. This parameter takes precedence over `qname` if both are 
#'   provided.
#' @param my_reference A character string specifying the name of the reference 
#'   genome to be used (e.g., "BSgenome.Hsapiens.UCSC.hg38"). The genome must be
#'   available in the `BSgenome` package.
#' @return A list containing:
#'   \item{reference_sequence}{A `DNAString` object representing the aligned reference sequence.}
#'   \item{aligned_sequence}{A `DNAString` object representing the aligned read sequence.}
#'   \item{original_sequence}{The original sequence from the read as a character string.}
#' @details 
#' If the `qname` parameter is specified, the function will search for reads with this `qname`. If more than 
#' one read is found, the function will notify the user and return `NULL` if the `index` parameter is not used.
#' If the `index` parameter is used, it will directly select the read at the given index. The `qname` 
#' parameter will be ignored in this case. The `index` parameter must be a valid row index in the DataFrame.
#' 
#' @examples
#' # Create a small mock df_gg DataFrame
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
#'   df_gg <- ggBamLoader(system.file("extdata", "subset.bam", package = "GGCigarAligner"))
#'   
#'   # Align read using qname
#'   result1 <- ggCigarAligner(df_gg, qname = "ERR188273.4711308", my_reference = "BSgenome.Hsapiens.UCSC.hg38")
#'   
#'   # Align read using index
#'   result2 <- ggCigarAligner(df_gg, index = 3, my_reference = "BSgenome.Hsapiens.UCSC.hg38")
#'   
#' } else {
#'   message("Skipping example as BSgenome.Hsapiens.UCSC.hg38 is not installed.")
#' }
#' @export

#' @importFrom BSgenome getBSgenome getSeq
#' @importFrom Biostrings subseq DNAStrings
#' @importFrom stringr str_extract_all

ggCigarAligner <- function(df_gg, qname = NULL, index = NULL, my_reference) {
  
  # ----------------------------------------------------------
  # Extract the reference genome
  if (!requireNamespace("BSgenome", quietly = TRUE)) {
    stop("The 'BSgenome' package is required but not installed. Please install it before proceeding.")
  }
  
  if (!requireNamespace(my_reference, quietly = TRUE)) {
    stop(paste("The", my_reference, "package is required but not installed. Please install it before proceeding."))
  }
  
  # Load the reference genome
  genome <- BSgenome::getBSgenome(my_reference)
  
  
  # ----------------------------------------------------------
  # Select the read(s) based on the input
  if (!is.null(index)) {
    # Use the index to select the specific read
    if (index > nrow(df_gg) || index <= 0) {
      stop("Index is out of bounds.")
    }
    selected_read <- df_gg[index, ]
  } else if (!is.null(qname)) {
    # Filter reads by qname
    selected_reads <- df_gg[df_gg$qname == qname, ]
    
    if (nrow(selected_reads) == 0) {
      stop("No reads found with the specified qname.")
    } else if (nrow(selected_reads) > 1) {
      cat("Multiple reads found with the same qname. Please use the index parameter to specify the read.\n")
      return(NULL)
    } else {
      selected_read <- selected_reads
    }
  } else {
    stop("Either qname or index must be provided.")
  }
  
  
  # ----------------------------------------------------------
  # Extract necessary information from the selected read
  cigar <- selected_read$cigar
  strand <- selected_read$strand
  pos <- selected_read$pos
  rname <- selected_read$rname
  
  # Print diagnostic information
  cat("Strand: ", strand, "\n")
  cat("CIGAR: ", cigar, "\n")
  
  
  # ----------------------------------------------------------
  # Error checking for strand and position
  if (!(strand %in% c("+", "-"))) {  
    stop("The provided read entry presents an invalid strand notation. Only \"+\" and \"-\" allowed.")
  }
  
  if (!is.numeric(pos) || pos <= 0) {
    stop("The provided read entry start position is not a valid integer.")
  }
  
  
  # ----------------------------------------------------------
  # Parsing the CIGAR string to extract the lengths of the operations
  cigar_matches <- str_extract_all(cigar, "[0-9]+[MIDNSHP=X]")[[1]]
  
  # Calculate the length of the reference sequence based on the CIGAR string
  running <- 0
  for (op in cigar_matches) {
    length <- as.numeric(gsub("[^0-9]", "", op))
    running <- running + length
  }
  
  cat("The reference sequence will be", running, "bases long. \n")
  
  
  # ----------------------------------------------------------
  # Extract the reference sequence from the reference genome
  ref <- BSgenome::getSeq(genome, rname, pos, pos + running - 1)
  
  # Extract the read sequence
  seq <- selected_read$seq[[1]]
  
  # Initialize aligned sequences
  reference_sequence <- ""
  aligned_sequence <- ""
  modificabile <- selected_read$seq[[1]]
  
  # Temporary variables for positions
  seq_position <- 1
  ref_position <- 1
  
  
  
  # ----------------------------------------------------------
  # Real loop translator: aligns the read based on the CIGAR string
  for (op in cigar_matches) {
    
    length <- as.numeric(gsub("[^0-9]", "", op))
    operation <- substr(op, nchar(op), nchar(op))
    
    if (operation == "M" || operation == "=" || operation == "X") {
      # Matching segment: copy from both the reference and the read
      reference_sequence <- paste(reference_sequence, subseq(ref, ref_position, ref_position + length - 1), sep = "")
      
      # Update positions by current match
      seq_position <- seq_position + length
      ref_position <- ref_position + length
    } 
    
    else if (operation == "D" || operation == "N") {
      # Deletion (and skipped/splitted regions) in the read: add gaps to the read sequence
      temp <- list( subseq(modificabile, 1, seq_position - 1), 
                    subseq(modificabile, seq_position, length(modificabile)))
      
      modificabile <- DNAString(paste(temp[[1]], 
                                      paste(rep("-", length), collapse = ""), 
                                      temp[[2]], sep = ""))
      
      
      reference_sequence <- paste(reference_sequence, subseq(ref, ref_position, ref_position + length - 1), sep = "")
      
      # Update reference position only
      ref_position <- ref_position + length
      seq_position <- seq_position + length
    } 
    
    else if (operation == "I" || operation == "S" || operation == "H" || operation == "P") {
      # Insertion in the read: add gaps to the reference sequence
      reference_sequence <- paste(reference_sequence, paste(rep("-", length), collapse = ""), sep = "")
      
      # Update read position only
      # We don't want to increase the current position by the length of the
      # operation, otherwise we'll lose part of the next bases.
      seq_position <- seq_position + length
    } 
    
    else {
      stop("Unknown CIGAR operation: ", operation)
    }
  }
  
  # ----------------------------------------------------------
  # Return or print the aligned sequences
  
  return(list(reference_sequence = DNAString(reference_sequence),
              aligned_sequence = DNAString(modificabile)))
  
}

