# epicypher-fiberseq-pipeline

Download the files and place them in a directory. It will contain the nextflow script (Fiberseq_docker_pipeline.nf) as well as the dockerfiles needed to create some of the docker images. After unzipping the file, move into the folder and run (you may need to use sudo as well):

Setting up the input directories
To set up the input directory for Fiber-seq in a pipeline compatible format, use the above script. It softlink and name the files according to the Run Design. You must also download the Run Design from SMRT Link. To download the Run Design, follow the steps under "Downloading the Run Design". 
If running the setup script for the first time, it will prompt you to establish a default raw data directory (i.e. where the data transfer scheme points to or where the raw data currently is) as well as an project analysis directory (i.e. where all the analysis should occur)
python3 prepFiberDir.py -r $RUN_NAME -d $RUN_DESIGN
Creating an input samplesheet
One of the required inputs is a samplesheet. The format must match that of the example below:

The first three columns can be directly copied and pasted from the run design (details on where to find that are below in the section "Downloading the Run Design"). The "Merge Sample Name" column is used to specify both which samples to merge and what to name the new merged sample. Any samples that should not be merged should be left blank. The samplesheet must be saved as a tab-separated file (.tsv).
Kicking off the Pipeline

Required parameters
--input_bam_path         : Path to the folder of the bams you would like to run the pipeline
--sample_sheet           : Path to samplesheet created from run design (must be tsv)
 
Optional Parameters
--outdir                 : Path that you would like the outputs to be written to (default:
                           "./fiberseq_output"
--reference_genome       : Reference genome to use for all samples (default: hg38, 
                           options: hg38, T2T, mm10)
--confidence_ml_val      : The confidence threshold of the methylation caller (default: 250)
--minimum_msp_dist       : The minimum distance of MSP calls to keep (default: 10)

# Optional Process Parameters
--pb_qc                  : Specify to run the PacBio QC script (default: False)
--phase_bams             : Run the variant calling and read phasing steps (default: False)
--create_bigwigs         : Run the pileup and bigwig creation steps (default: False)

# Nextflow specific parameters
-with-report             : Generate workflow report for pipeline. Shows resource usage.
-with-timeline           : Generate pipeline timeline report
-resume                  : Extremely useful. If pipeline fails, starts it up again where it         a                          left off

 
Example commands:
# The simplest way to run the pipeline. Uses hg38 reference. Will put outputs in current directory in folder "fiberseq_output". This will also not run the optional PacBio QC step.

nextflow run /path/to/Fiberseq_docker_pipeline.nf --input_bam_path /path/to/bam_directory/

# The recommended way to run the pipeline. Generates all the necessary reports and uses only the most useful parameters

nextflow run /path/to/Fiberseq_docker_pipeline.nf --input_bam_path /path/to/bam_directory/ --outdir /path/to/output_dir/ --reference_genome T2T -with-report -with-timeline --pb_qc

Downloading the Run Design
1.) Navigate to the "Runs" module on SMRT Link


2.) Select the sequencing run you would like to analyze by clicking the blue name under "Run Name":


3.) Select "View Run Design" in the top right corner of the window:


4.) Click the "Export Run Design" button to download the run design to your default browser location:



Summarizing the QC
After running the pipeline, run the command below to create a summary QC table that includes the commonly reported statistics. If you would like more statistics added to the table, reach out to Connor. This script will use the same configuration that was created during the "Setting up the input directories" step.
python3 collectPBMetrics.py -r $RUN_NAME
