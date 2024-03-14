# The Assemble_fastq.R file from MAESTER.
# Originally written by Peter van Galen, 200110
# Goal: append cell barcode and unique molecular identifier, _CB_UMI, from Read 1 to each Read 2 identifier


# Load libraries.
suppressMessages(library(ShortRead))
suppressMessages(library(optparse))

# We clean the session to make sure that no previous data interferes with the script.
rm(list=ls())

# Functions.
cutf <- function(x, f=1, d="/", ...) sapply(strsplit(x, d), function(i) i[f], ...)

# We get the various input variables.
option_list = list(
  optparse::make_option("--Input_Folder_Path",  type = "character", default = "", help = "The input folder.", metavar = "character"),
  optparse::make_option("--Sample",             type = "character", default = "", help = "The sample to be analyzed.", metavar = "character"),
  optparse::make_option("--Cell_Barcodes_Path", type = "character", default = "", help = "The path to the cell barcodes file.", metavar = "character"),
  optparse::make_option("--CB_Length",          type = "integer",   default = 16, help = "The length of the cell barcode. 10X uses a length of 16.", metavar = "integer"),
  optparse::make_option("--UMI_Length",         type = "integer",   default = 12, help = "The length of the UMI. 10X uses a length of 12.", metavar = "integer"),
  optparse::make_option("--Output_Folder_Path", type = "character", default = "", help = "The output folder.", metavar = "character")
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


# Values.
Input_Folder_Path  <- opt$Input_Folder_Path
Sample_Name        <- opt$Sample
Cell_Barcodes_Path <- opt$Cell_Barcodes_Path
CB_Length          <- opt$CB_Length
UMI_Length         <- opt$UMI_Length
Output_Folder_Path <- opt$Output_Folder_Path

# We create the output folder, if it does not exist yet.
dir.create(Output_Folder_Path, showWarnings = FALSE, recursive = TRUE)
setwd(Output_Folder_Path)


print("Values used.")
print(paste0("Input_Folder_Path:  ", Input_Folder_Path))
print(paste0("Sample_Name:        ", Sample_Name))
print(paste0("Cell_Barcodes_Path: ", Cell_Barcodes_Path))
print(paste0("CB_Length:          ", CB_Length))
print(paste0("UMI_Length:         ", UMI_Length))
print(paste0("Output_Folder_Path: ", Output_Folder_Path))


print("We read in the cell barcodes.")
Cell_Barcodes <- read.table(Cell_Barcodes_Path, sep = "\t", header = FALSE)[,1]


print("Find R1 fastq files (R2 is found by substitution later).")
R1.ch <- list.files(Input_Folder_Path, pattern = paste0(Sample_Name, ".*_R1_.*fastq.gz$"), full.names = TRUE)
message(Sys.time(), "\nLoading ", length(R1.ch)*2, " fastq files:")
message(cat(c(R1.ch, sub("_R1_", "_R2_", R1.ch)), sep = "\n"))
if(length(R1.ch) == 0) stop("Did not find fastq files.")


print("Filter for cell barcodes. Remove -1 from the end (if added by CellRanger count).")
cells.ch <- cutf(CellBarcodes, d = "-", f = 1)
if(length(cells.ch) == 0) stop("No cells found.")
message("Found ", length(cells.ch), " cells.\n")


# Process the fastq files.
message("Read R1 and R2 sequences, filter by ", length(cells.ch), " cell barcodes, write assembled fastq...")

# For each R1 fastq
for(f1 in R1.ch) {
  # Identify R2 fastq
  f2 <- sub("_R1_", "_R2_", f1)

  # Load file in 1E7 read increments
  message("file ", match(f1, R1.ch), "/", length(R1.ch), ": ", basename(f1), " ", appendLF = FALSE)
  strm1 <- ShortRead::FastqStreamer(f1, n=1E7) # 1M reads by default
  strm2 <- ShortRead::FastqStreamer(f2, n=1E7)

  # For every 10 million reads do the following.
  repeat{
    message("*", appendLF = FALSE)
    fq1 <- ShortRead::yield(strm1)
    fq2 <- ShortRead::yield(strm2)
    if(length(fq1) == 0 | length(fq2) == 0) break

    # Match to expected cell barcodes
    fq1.m <- ifelse(is.element(as.vector(XVector::subseq(ShortRead::sread(fq1), 1, CB_Length)), cells.ch), yes = TRUE, no = FALSE)

    # Filter unmatched reads from the ShortRead objects
    fq1.f <- fq1[fq1.m]
    fq2.f <- fq2[fq1.m]

    # Extract cell barcode and umi from Read1
    fq1.f.cell <- as.vector(XVector::subseq(ShortRead::sread(fq1.f), 1, CB_Length))
    fq1.f.umi  <- as.vector(XVector::subseq(ShortRead::sread(fq1.f), CB_Length + 1, CB_Length + UMI_Length))

    # Add cell barcode and umi to id of Read2
    fq2.f@id <- Biostrings::BStringSet(paste0(sub(" .:N:0:", "_", as.vector(fq2.f@id)), "_", fq1.f.cell, "_", fq1.f.umi))
    
    # Check if all the ids of Read1 and Read2 files match up
    if(!all(cutf(as.vector(fq1.f@id), d = " ", f = 1) == cutf(as.vector(fq2.f@id), d = "_", f = 1))) stop("Read ID mismatch")

    # Save fastq file
    ShortRead::writeFastq(fq2.f, file = paste0(Output_Folder_Path, Sample_Name, ".fastq.gz"), mode = "a")
  
    # Clean up memory
    invisible(gc())
  }

  close(strm1)
  close(strm2)
  message(" done")
  invisible(gc())
}

message("\nFinished!")
