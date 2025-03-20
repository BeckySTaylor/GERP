
#########################
#   1. Analysis Setup   #
#########################

# 1.1 Analysis Options #
# Before running this script, set the variables in this code block:
{

        path=/path/to/directory/ # The directory in which to run the analysis - this directory should include:
                # the control file; and
                # all genome fasta files to be included in the analysis.

        control=Control.tsv # The name of the control file.

        # !! If you do not wish to change the reference genome used (i.e. you are using the caribou reference
        # !! genome discussed in the paper), un-comment the following block of code to download and prepare the file:
                # wget https://caribougenome.s3.us-east-2.amazonaws.com/Dovetail_ref_May2021_fasta.tar.gz
                # mv Dovetail_ref_May2021_fasta.tar.gz $path/Dovetail_ref_May2021_fasta.tar.gz
                # cd $path
                # tar xzf Dovetail_ref_May2021_fasta.tar.gz

        # If you do wish to change the reference genome used, put the file in the directory specified by your 'path'
        # variable and set the 'ref' variable to that filename:
        ref=Reference_genome.fasta

        # !! the number of scaffolds you would like to analyze (i.e., the number of autosomes in the subject genome) minus 1 as the counting starts at 0.
	chr=34

 	#How many cores will be used, specifically for the alignment step in bwa.
        cores=10

}


# If the number of scaffolds in Section 3 is correct, and -t for the bwa mem has the correct number of cores to be used, no more changes to the script are necessary unless changes to parameters are desired.

# 1.2 Analysis Prep #
cd $path # Change working directory to the one provided above.
bwa index $ref # Index the reference genome.

#######################################################################
#   2. Align To Reference Genome And Generate Scaffold Fasta Files    #
#######################################################################

# Step 2 of the analysis is contained within the following loop:
# The loop iterates over each genome file and species name listed in the control file.
cat $control | while read genome species; do

        # 2.1 Directory Setup #
        mkdir $species # Create a directory for the current species.
        cd $species # Change working directory to the species-specific one we just created.

        # 2.2 Generate fastq file #
        reformat.sh \
                in=../${genome} \
                out1=${species}.fastq \
                qfake=40 \
                fastareadlen=5000 \
                qout=64 \
                addcolon=t \
                trimreaddescription=t \
                int=t

        # 2.3 Generate BAM file #
        bwa mem \
                -t ${cores} \
                -B 3 \
                -O 4,4 \
                ../${ref} \
                ${species}.fastq | \
        samtools view \
                -F 2048 \
                -bq 2 \
                -h | \
        samtools sort \
                -o ${species}_filtered_sorted.bam

        # 2.4 Generate one fasta file per scaffold #
samtools consensus \
                -f fasta \
                --min-MQ 30 \
                --min-BQ 30 \
                -o ${species}_filtered_sorted.fasta ${species}_filtered_sorted.bam
        csplit \
                -s -n 1 --prefix=${species}_Scaffold \
                -z ${species}_filtered_sorted.fasta '/>/' '{*}'

        # 2.5 Rename header line within the fasta files to have the species name #
        ls ${species}_Scaffold* | \
        xargs \
                -n1 \
                sed \
                         -i "1s/.*/>${species}/"
        cd $path

done


################################################
#   3. Create Scaffold Multialignment Files    #
################################################

# For each scaffold, create a file that contains the sequence belonging to all species in the analysis:
for i in $(seq 0 ${chr}); do

        cat */*_Scaffold$i | \
        awk '/^>/ {printf("\n%s\n",$0);next } { printf("%s",$0)} END {printf("\n")}'| \
        tail -n+2 > Scaffold${i}_GERP_formatted.mfa

done
