\name{LoadingVarTrix}
\alias{LoadingVarTrix}
\title{
Loading VarTrix results for the down stream analysis.
}
\description{
This function loads the VarTrix results and outputs a list of two SummarizedExperiments objects, one for the somatic variants and on accompanying for the mitochondrial variants.
}
\usage{
LoadingVarTrix(samples_file, vcf_path, vcf_path_MT, patient)
}
\arguments{
\item{\bold{samples_file}}{Path to the input file. One sample per row. It has to contain the following columns:\cr\cr
\bold{patient}: The patient of the sample. A patient can have multiple samples (treated/untreated). All samples are merged per patient.\cr\cr
\bold{sample}: The ID of the sample.\cr\cr
\bold{resource}: Is this a VarTrix or MAEGATK input? Must be either VarTrix or MAEGATK.\cr\cr
\bold{type}: What type of data? Must be one of the following: scRNAseq_Somatic, scRNAseq_MT, Amplicon_Somatic, Amplicon_MT.\cr\cr
\bold{input_folder}: Path to the folder, where the VarTrix output is stored.\cr\cr
\bold{cells}: Path to the barcodes file from the CellRanger output.
}
\item{\bold{vcf_path}}{Path to the VCF file that contains the somatic variants you are interested in.}
\item{\bold{vcf_path_MT}}{Path to the VCF file that contains the MT variants you are interested in.}
\item{\bold{patient}}{The patient you want to load.}
}
\details{
The function loads the VarTrix for one of the patients in the input file.\cr
If a patient has multiple samples, then all samples are merged.\cr
There are two types of input data, single cell RNA sequencing and amplicon data.\cr
scRNAseq data is the result from something like a 10X assay, while amplicon data focuses on a specific part of the genome.\cr
Both are the output of the CellRanger pipeline. Amplicon data has a much higher coverage in the specific area and therefore is much better for these areas. The scRNAseq data is overwritten if amplicon data is available.\cr\cr
The output is a list of two SummarizedExperiment objects. Each object contains two assays:\cr\cr
\bold{consensus_somatic}: The consensus information per somatic variant and cell.\cr
\bold{consensus_MT}: The consensus information per MT variant and cell.\cr
\bold{fraction_somatic}: The variant allele frequency per somatic variant and cell.\cr
\bold{fraction_MT}: The variant allele frequency per MT variant and cell.\cr\cr
The consensus information is encoded in the following way:\cr\cr
\bold{0}: No Call, no reads are covering this position.\cr
\bold{1}: Reference, only reference reads cover this position.\cr
\bold{2}: Alternative, only mutated reads cover this position.\cr
\bold{3}: Both, both alt and ref reads cover this position.\cr
}
\value{
The return list has two objects.\cr\cr
\item{\bold{SE_Somatic}}{The somatic variants results.}
\item{\bold{SE_MT}}{The MT variants results.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
