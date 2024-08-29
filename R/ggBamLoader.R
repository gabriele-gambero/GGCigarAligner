#' Load a BAM/SAM file and convert it to a DataFrame.
#'
#' This function loads a BAM/SAM file, checks if the format is correct, and 
#' converts it to an S4 class DataFrame for easier manipulation.
#'
#'
#' @param bampath A string representing the path to the input BAM or SAM file.
#' @return An S4 DataFrame containing the mapped reads from the file.
#' @export
#' @importFrom Rsamtools scanBam asBam
#' @importFrom S4Vectors DataFrame
#' @examples
#' # Load a BAM file from the package's extdata folder
#' df_gg <- ggBamLoader(system.file("extdata", "subset.bam", package = "GGCigarAligner"))


ggBamLoader <- function(bampath) {
  
  # ----------------------------------------------------------
  # File Validation and Importing
  
  # Check if the file exists
  if (!file.exists(bampath)) {
    stop("The file does not exist. Please provide a valid file path.")
  }
  
  # Check file format
  if (endsWith(bampath, ".sam")) {
    # Convert SAM to BAM and then read
    ggbam <- Rsamtools::scanBam(Rsamtools::asBam(bampath))
  } else if (endsWith(bampath, ".bam")) {
    # Directly read BAM file
    ggbam <- Rsamtools::scanBam(bampath)
  } else {
    stop("Your file is not in .sam or .bam format. Please provide a valid file."
         )
  }
  
  # Check for empty BAM data
  if (length(ggbam) == 0) {
    stop("The BAM/SAM file appears to be empty or not properly formatted.")
  }
  
  # Convert to DataFrame
  df_gg <- S4Vectors::DataFrame(ggbam[[1]])
  
  # Output the DataFrame
  return(df_gg)
}
