#'UMIQC
#'@description
#'We load the molecule_info.h5 file from CellRanger and plot the molecules by rank.
#'A molecule is defined as a combination of UMI and CB.
#'The molecules are ranked by their number of replicates.
#'@importFrom DropletUtils read10xMolInfo
#'@param samples_file_path Path to the csv file with the samples to be loaded.
#'@param molecule_info_path Path to a single molecule_info file.
#'@param umi_length The length of the UMI. Default = 12
#'@param min_reads The minimum number of reads per molecule.
#'@param sep The separator used in your samples_file. Default: ","
#'@param samples_column The column for the samples in your samples_file. Default: "sample"
#'@param molecule_info_column The column for the paths to the molecule_info files in your samples_file. Default: "molecule_info"
#'@export
UMIQC <- function(samples_file_path = NULL, molecule_info_path = NULL, umi_length = 12, min_reads = 100, sep = ",", samples_column = "sample", molecule_info_column = "molecule_info"){
  if(!is.null(molecule_info_path)){
    if(!file.exists(molecule_info_path)){
      stop("Your molecule_info_path does not point to a file.")
    } else{
      molecule_info <- DropletUtils::read10xMolInfo(molecule_info_path)
      molecule_info <- molecule_info$data
      molecule_info$umi <- UMI_seq(number = molecule_info$umi, umi_length = umi_length)
      molecule_info <- data.frame(cell = molecule_info$cell, umi = molecule_info$umi, molecule = paste0(molecule_info$cell, "_", molecule_info$umi), reads = molecule_info$reads)
      molecule_info <- molecule_info[base::order(molecule_info$reads, decreasing = TRUE),]
      molecule_info$ranks <- 1:nrow(molecule_info)
      molecule_info <- subset(molecule_info, reads >= min_reads)
      cells_retained <- length(unique(molecule_info[,1]))
      umis_per_cell <- lapply(molecule_info[,1], function(x){
        res <- subset(molecule_info, molecule_info$cell == x)
        umis <- nrow(res)
        names(umis) <- x
        return(umis)
      })
      umis_per_cell <- unlist(umis_per_cell)
      sample_summary <- data.frame(Cells = cells_retained, Average_UMIs_Per_Cells = mean(umis_per_cell), Average_Reads_Per_UMI = mean(molecule_info[,"reads"]))
      molecule_info <- list(sample_summary = sample_summary, molecule_info = molecule_info)
    }
  }
  if(!is.null(samples_file_path)){
    if(file.exists(samples_file_path)){
      samples_file <- read.table(samples_file_path, header = TRUE, sep = sep)
      if(!molecule_info_column %in% colnames(samples_file)){
        stop("Your molecule_info_column is not in the columns of your samples file.")
      } else{
        molecule_info <- lapply(1:nrow(samples_file), function(x){
          mol_info <- DropletUtils::read10xMolInfo(samples_file[x, molecule_info_column])$data
          mol_info$umi <- UMI_seq(number = mol_info$umi, umi_length = umi_length)
          mol_info <- data.frame(Sample = samples_file[x, samples_column], cell = mol_info$cell, umi = mol_info$umi, molecule = paste0(mol_info$cell, "_", mol_info$umi), reads = mol_info$reads)
          mol_info <- mol_info[base::order(mol_info$reads, decreasing = TRUE),]
          mol_info$ranks <- 1:nrow(mol_info)
          mol_info <- subset(mol_info, mol_info$reads >= min_reads)
          cells_retained <- length(unique(mol_info[,"cell"]))
          umis_per_cell <- lapply(mol_info[,"cell"], function(x){
            res <- subset(mol_info, mol_info$cell == x)
            umis <- nrow(res)
            names(umis) <- x
            return(umis)
          })
          umis_per_cell <- unlist(umis_per_cell)
          sample_summary <- data.frame(Sample = samples_file[x, samples_column], Cells = cells_retained, Average_UMIs_Per_Cells = mean(umis_per_cell), Average_Reads_Per_UMI = mean(mol_info[,"reads"]))
          return(list(sample_summary = sample_summary, molecule_info = mol_info))
        })
        molecule_info <- list(sample_summary = do.call("rbind", lapply(molecule_info, "[[", 1)), 
                              molecule_info = do.call("rbind", lapply(molecule_info, "[[", 2)))
      }
    }
  }
  return(molecule_info)
}

