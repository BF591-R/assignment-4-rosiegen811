library('tidyverse')
library('RColorBrewer')

#' Read the expression data "csv" file as a dataframe, not tibble
#'
#' @param filename (str): the path of the file to read (intensity_data)
#' @param delimiter (str): generalize the function so it can read in data with
#'   your choice of delimiter
#'
#' @return A dataframe containing the example intensity data with rows as probes
#'   and columns as samples
#' @export
#'
#' @examples
read_data <- function(intensity_data, delimiter) {
    df <- read.table(
      file = intensity_data,
      sep = delimiter,
      header = TRUE,
      row.names = 1,  # Make the probe ID column the row names.
      check.names = FALSE,  # Keep probe ID names as they appear.
      stringsAsFactors = FALSE
    )
    return(df)
}

#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples
calculate_variance_explained <- function(pca_results) {
  
    # Formula: Variance explained by each PC = (PC's SD)^2 / Total variance
    var_explained <- (pca_results$sdev^2) / sum(pca_results$sdev^2)
  
    return(var_explained)
}

#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PC names, variance explained by each PC, and the
#' cumulative sum of variance explained. These columns should be named 
#' "principal_components", "variance_explained", and "cumulative", respectively.
#' 
#'
#'
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained, and the cumulative variance explained with names described above
#' @export
#' @examples 
make_variance_tibble <- function(pca_ve, pca_results) {
  
    # Create a tibble by columns.
    var_tibble <- tibble(
      # Generate PC names based on the length of pca_ve.
      principal_components = paste0("PC", seq_along(pca_ve)),
      variance_explained = pca_ve,
      cumulative = cumsum(pca_ve)
    )
    
    return(var_tibble)
}

#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples
make_biplot <- function(metadata, pca_results) {
    # Get the PC coordinates from pca_results.
    sample_coords <-  pca_results$x
  
    # Turn PC coordinates into a data frame.
    pca_df <- as_tibble(sample_coords, rownames = "SampleID")
    
    # Bring in the metadata.
    metadata_df <- read.csv(metadata, stringsAsFactors = FALSE)
    plot_df <- left_join(pca_df, metadata_df, by = c("SampleID" = "geo_accession"))
    
    # Plot PC1 vs PC2.
    # Map SixSubTypesClassification to color.
    pc_plot <- ggplot(plot_df, aes(x = PC1, y = PC2, color = SixSubtypesClassification)) +
      geom_point() + 
      labs(
        x = "PC1",
        y = "PC2",
        color = "SixSubTypesClassification"
      )
    
    return(pc_plot)
}

#' Define a function to return a list of probeids filtered by signifiance
#'
#' @param diff_exp_tibble (tibble): A tibble containing the differential expression results
#' @param fdr_threshold (float): an appropriate FDR threshold, we will use a
#'   value of .01. This is the column "padj" in the tibble.
#'
#' @return A list with the names of the probeids passing the fdr_threshold
#' @export
#'
#' @examples
list_significant_probes <- function(diff_exp_tibble, fdr_threshold) {
    
    sig_probes <- diff_exp_tibble |> 
      # Filter rows where padj < fdr_threshold.
      filter(padj < fdr_threshold) |> 
      
      # Extract the probe IDs.
      pull(probeid)

    return(sig_probes)
}

#' Define a function that uses the list of significant probeids to return a
#' matrix with the intensity values for only those probeids.
#' @param intensity (dataframe): The dataframe of intensity data generated in
#'   part 1
#' @param sig_ids_list (list/vector): The list of differentially expressed
#'   probes generated in part 6
#'
#' @return A `matrix()` of the probe intensities for probes in the list of
#'   significant probes by FDR determined in the previous function.
#'
#' @export
#'
#' @examples
return_de_intensity <- function(intensity, sig_ids_list) {
    # Select the rows whose names match the significant probe IDs.
      # Keep the subset 2D even if there's one row.
    de_intensity <- intensity[sig_ids_list, , drop = FALSE]
    
    # Return a matrix of the probe intensities.
    return(as.matrix(de_intensity))
}

#' Define a function that takes the intensity values for significant probes and
#' creates a color-blind friendly heatmap
#'
#' @param de_intensity (matrix): The matrix of intensity values for significant
#'   differentially expressed probes returned in part 7
#' @param num_colors (int): The number of colors in a specificed RColorBrewer
#'   palette
#' @param palette (str): The name of the chosen RColorBrewer palette
#'
#' @return A heatmap displaying the intensity values for the differentially
#'   expressed probes
#' @export
#'
#' @examples
plot_heatmap <- function(de_intensity, num_colors, palette) {
  
    # Make a vector of colors.
    heatmap_colors <-  brewer.pal(num_colors, palette)
  
    # Pass the matrix into heatmap().
    hm <- heatmap(
      de_intensity,
      col = heatmap_colors
    )
}

