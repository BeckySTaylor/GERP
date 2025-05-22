# Modified GERP program
A modified Genomic Evolutionary Rate Profiling (GERP) pipeline for the analysis of genetic load in wildlife species

The original version of GERP++ is available here: https://github.com/tvkent/GERPplusplus

To facilitate the analysis of genetic load, we provide a script to automate input file production called 'Pipeline_input_alignment_GERP.sh', as well as a modified version of GERP++ which creates more user-friendly outputs and streamlines the extraction of sites of interest from a VCF file as known derived alleles using three outgroup species, a step which can be completed using an R script we provide called 'Derived_alleles_extract.R'.

Specifically, our modified version can handle missing data in the alignments and outputs the position for each score, as well as the allele for three specified sister species to enable the detection of derived alleles. The outputs also have a header row included to indicate what each row pertains to.

The protocol outlining the extact steps involved and how to do them is published in STAR Protocols (availabe here: https://star-protocols.cell.com/protocols/4202), a PDF version is also included here. Please cite this if using our pipeline. Some example files, including control files for both the bash and R scripts, are included in the 'Example_files' folder.

We have included the 'gerpcol.cc' file which is the part of the program we have modified so users can see the code and track any future changes. However, it is important to download 'gerp_modified.tar.gz' and compile the prgram fully for usage.
