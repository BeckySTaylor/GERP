# Modified GERP program
A modified Genomic Evolutionary Rate Profiling (GERP) pipeline for the analysis of genetic load in wildlife species

The original version of GERP++ is available here: https://github.com/tvkent/GERPplusplus

To facilitate the analysis of genetic load, we provide a script to automate input file production called 'Pipeline_input_alignment_GERP.sh', as well as a modified version of GERP++ which creates more user-friendly outputs and streamlines the extraction of sites of interest from a VCF file as known derived alleles using three outgroup species, a step which can be completed using an R script we provide called 'Derived_alleles_extract.R'.

Specifically, our modified version can handle missing data in the alignments and outputs the position for each score, as well as the allele for three specified sister species to enable the detection of derived alleles. The outputs also have a header row included to indicate what each row pertains to.

The protocol outlining the extact steps involved and how to do them is also included in the document named 'Taylor_etal_Protocol_modified_GERP', and an example control file is included named 'Control_example.tsv'.
