
# Running GATK in Docker for Advanced Clinical Bioinformatics 2019

**These notes assume you are running a Linux environment.  The commands can be copied direct from this file but may need amending depending on your file names and locations.** 

### Setting Up & Running Docker

First make a folder to hold your files on your local machine:
 ```bash
 mkdir AdvClinBio2019
 cd AdvClinBio2019/
```

Pull down the gatk4 docker image (for instructions on installing docker please see [here](https://docs.docker.com/engine/install/ubuntu/):
```bash
docker pull broadinstitute/gatk
```
Use ```wget``` to download the GATK bundle from the Broad Institue website (this will be several gigabytes of additional data required for the analysis and may take some time to download):
```bash
wget -pr ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/
```
Download the aligned bam file from Galaxy and save it to your working directory (this should have also been processed with Picard's ```MarkDuplicates``` tool before entering the GATK pipeline).  My file is called 'MarkDuplicates_BAM_output.bam' 

We will map our current working directory to our docker image so that we have access to these files from within the docker container:
 ```bash
sudo docker run -v `pwd`:/home  -i -t broadinstitute/gatk
```

To map a local drive to the container we use the format `-v /local_directory:/container_directory`.  The command \`pwd\` basically returns our current working directory, rather than us having to specify it directly, and maps it to the `/home` directory in our docker container.  This saves us having to copy all the required files manually to the container.  The three flags we use in this command are: 
- `--volume , -v`	Bind mount a volume
- `--interactive , -i`		Keep STDIN open even if not attached
- `--tty , -t`		Allocate a pseudo-TTY

Your prompt should now change to indicate that you you have entered an interactive session with the docker container:
`(gatk) root@d6c83ac31f5a:/gatk#` 
We can get our current directory by using the `pwd` command:
```bash
pwd
/gatk
```
We are in the `/gatk` folder (this information could also be provided by the prompt depending on how your system is configured).  Lets check whether our mounted files are available to use within the container,  remember that they are in our `/home` directory.
```bash
ls -lh /home
total 525M
-rw-rw-r-- 1 1000 1000 525M May 23 12:22 MarkDuplicates_BAM_output.bam
drwxr-xr-x 2 1000 1000  16K Apr 16 09:21 hg19
```
If all your files appear in your `/home` directory then you are ready to proceed with the analysis. If at any point you need to exit the docker container you can simply type `exit` in the terminal.

### Preparing the BAM file for analysis

We require these six files from the bundle:
/home/hg19/ucsc.hg19.fasta.gz
/home/hg19/ucsc.hg19.fasta.fai.gz
/home/hg19/ucsc.hg19.dict.gz
/home/hg19/hapmap_3.3.hg19.sites.vcf.gz
/home/hg19/hapmap_3.3.hg19.sites.vcf.idx.gz

Before proceeding will first assign the name of our BAM file to a variable in bash so that the code below will work for any file through the magic of parameter expansion:
```bash
bam_file="MarkDuplicates_BAM_output.bam"
```
Substitute `"MarkDuplicates_BAM_output.bam"` for the name of your BAM file. Every time I reference `$bam_file` on the bash commandline it will now be interpreted as `MarkDuplicates_BAM_output.bam`, or whatever filename you allocated to this variable. 

Lets check our BAM file by running the `CountReads` command below:
```bash
gatk CountReads -I /home/$bam_file
```
You should see output similar to below:
```bash
Running:
    java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /gatk/gatk-package-4.1.2.0-local.jar CountReads -I /home/MarkDuplicates_BAM_output.bam
22:40:15.286 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/gatk/gatk-package-4.1.2.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
May 23, 2019 10:40:17 PM shaded.cloud_nio.com.google.auth.oauth2.ComputeEngineCredentials runningOnComputeEngine
INFO: Failed to detect whether we are running on Google Compute Engine.
22:40:17.089 INFO  CountReads - ------------------------------------------------------------
22:40:17.089 INFO  CountReads - The Genome Analysis Toolkit (GATK) v4.1.2.0
22:40:17.089 INFO  CountReads - For support and documentation go to https://software.broadinstitute.org/gatk/
22:40:17.090 INFO  CountReads - Executing as root@d6c83ac31f5a on Linux v4.8.0-53-generic amd64
22:40:17.090 INFO  CountReads - Java runtime: OpenJDK 64-Bit Server VM v1.8.0_191-8u191-b12-0ubuntu0.16.04.1-b12
22:40:17.090 INFO  CountReads - Start Date/Time: May 23, 2019 10:40:15 PM UTC
22:40:17.090 INFO  CountReads - ------------------------------------------------------------
22:40:17.090 INFO  CountReads - ------------------------------------------------------------
22:40:17.091 INFO  CountReads - HTSJDK Version: 2.19.0
22:40:17.091 INFO  CountReads - Picard Version: 2.19.0
22:40:17.091 INFO  CountReads - HTSJDK Defaults.COMPRESSION_LEVEL : 2
22:40:17.091 INFO  CountReads - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
22:40:17.091 INFO  CountReads - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
22:40:17.091 INFO  CountReads - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
22:40:17.091 INFO  CountReads - Deflater: IntelDeflater
22:40:17.091 INFO  CountReads - Inflater: IntelInflater
22:40:17.091 INFO  CountReads - GCS max retries/reopens: 20
22:40:17.092 INFO  CountReads - Requester pays: disabled
22:40:17.092 INFO  CountReads - Initializing engine
22:40:17.568 INFO  CountReads - Done initializing engine
22:40:17.568 INFO  ProgressMeter - Starting traversal
22:40:17.568 INFO  ProgressMeter -        Current Locus  Elapsed Minutes       Reads Processed     Reads/Minute
22:40:32.307 INFO  CountReads - 10584740 read(s) filtered by: WellformedReadFilter 

22:40:32.312 INFO  ProgressMeter -             unmapped              0.2                     0              0.0
22:40:32.312 INFO  ProgressMeter - Traversal complete. Processed 0 total reads in 0.2 minutes.
22:40:32.312 INFO  CountReads - Shutting down engine
[May 23, 2019 10:40:32 PM UTC] org.broadinstitute.hellbender.tools.CountReads done. Elapsed time: 0.28 minutes.
Runtime.totalMemory()=683147264
Tool returned:
0
```
From this verbose output we can see that my BAM file hase 10584740 reads.

Before we can use this BAM with HaplotypeCaller we need to ensure that all the metadata which GATK expects using Picard tool's `AddOrReplaceReadGroups`.
Required Arguments:
- `--INPUT,-I:String`          Input file (BAM or SAM or a GA4GH url).  Required. 
- `--OUTPUT,-O:File`              Output file (BAM or SAM).  Required. 
- `--RGLB,-LB:String`             Read-Group library  Required. 
- `--RGPL,-PL:String`             Read-Group platform (e.g. illumina, solid)  Required. 
- `--RGPU,-PU:String`             Read-Group platform unit (eg. run barcode)  Required. 
- `--RGSM,-SM:String`             Read-Group sample name  Required. 
The exact values for some of this metadata depends on your dataset, below I use some reasonable names for the cardiac panel I used:
```bash
gatk AddOrReplaceReadGroups -I /home/$bam_file -O /home/rg_$bam_file -LB cardio -PL illumina -PU H8LHPADXX -SM cardio2 --VALIDATION_STRINGENCY=LENIENT
```
To check that this process has created a valid BAM file we use `ValidateSamFile`:
```bash
gatk ValidateSamFile -I /home/rg_$bam_file -MODE=SUMMARY
```
You should see the tool output someting like:
```bash
Using GATK jar /gatk/gatk-package-4.1.2.0-local.jar
Running:
    java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /gatk/gatk-package-4.1.2.0-local.jar ValidateSamFile -I /home/rg_MarkDuplicates_BAM_output.bam --MODE=SUMMARY
23:25:00.896 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/gatk/gatk-package-4.1.2.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
[Thu May 23 23:25:00 UTC 2019] ValidateSamFile  --INPUT /home/rg_MarkDuplicates_BAM_output.bam --MODE SUMMARY  --MAX_OUTPUT 100 --IGNORE_WARNINGS false --VALIDATE_INDEX true --INDEX_VALIDATION_STRINGENCY EXHAUSTIVE --IS_BISULFITE_SEQUENCED false --MAX_OPEN_TEMP_FILES 8000 --SKIP_MATE_VALIDATION false --VERBOSITY INFO --QUIET false --VALIDATION_STRINGENCY STRICT --COMPRESSION_LEVEL 2 --MAX_RECORDS_IN_RAM 500000 --CREATE_INDEX false --CREATE_MD5_FILE false --GA4GH_CLIENT_SECRETS client_secrets.json --help false --version false --showHidden false --USE_JDK_DEFLATER false --USE_JDK_INFLATER false
May 23, 2019 11:25:02 PM shaded.cloud_nio.com.google.auth.oauth2.ComputeEngineCredentials runningOnComputeEngine
INFO: Failed to detect whether we are running on Google Compute Engine.
[Thu May 23 23:25:02 UTC 2019] Executing as root@d6c83ac31f5a on Linux 4.8.0-53-generic amd64; OpenJDK 64-Bit Server VM 1.8.0_191-8u191-b12-0ubuntu0.16.04.1-b12; Deflater: Intel; Inflater: Intel; Provider GCS is available; Picard version: Version:4.1.2.0
WARNING	2019-05-23 23:25:02	ValidateSamFile	NM validation cannot be performed without the reference. All other validations will still occur.
INFO	2019-05-23 23:25:43	SamFileValidator	Validated Read    10,000,000 records.  Elapsed time: 00:00:41s.  Time for last 10,000,000:   41s.  Last read position: chr20:34,170,844
No errors found
[Thu May 23 23:25:46 UTC 2019] picard.sam.ValidateSamFile done. Elapsed time: 0.76 minutes.
Runtime.totalMemory()=595066880
Tool returned:
0
```
The key output is `No errors found`. Because we now create an index for the BAM file (I believe the `MarkDuplicates` tool used below requires an index to function):
```bash
gatk BuildBamIndex -I /home/rg_$bam_file
```
Lets check that worked:
```bash
ls -lh /home/
total 1.3G
-rw-rw-r-- 1 1000 1000 525M May 23 12:22 MarkDuplicates_BAM_output.bam
drwxrwxr-x 3 1000 1000 4.0K May 23 12:28 ftp.broadinstitute.org
drwxr-xr-x 2 1000 1000  16K Apr 16 09:21 hg19
-rw-r--r-- 1 root root 5.6M May 23 23:29 rg_MarkDuplicates_BAM_output.bai
-rw-r--r-- 1 root root 711M May 23 23:21 rg_MarkDuplicates_BAM_output.bam
```
We  can then use Picard's `MarkDuplicates` tool to remove duplicate reads which represent the same genomic fragment so that we do not over estimate our coverage.
```bash
gatk MarkDuplicates -I /home/rg_$bam_file -O /home/md_rg_$bam_file -M /home/markup_metrics.txt
```
As we have changed the BAM again we need to create a new index for the BAM file:
```bash
gatk BuildBamIndex -I /home/md_rg_$bam_file
```
Lets check that we have all the files we expect:
```bash
ls -lh /home/
total 2.0G
-rw-rw-r-- 1 1000 1000 525M May 23 12:22 MarkDuplicates_BAM_output.bam
drwxrwxr-x 3 1000 1000 4.0K May 23 12:28 ftp.broadinstitute.org
drwxr-xr-x 2 1000 1000  16K Apr 16 09:21 hg19
-rw-r--r-- 1 root root 5.6M May 23 23:38 md_rg_MarkDuplicates_BAM_output.bai
-rw-r--r-- 1 root root 712M May 23 23:36 md_rg_MarkDuplicates_BAM_output.bam
-rw-r--r-- 1 root root 5.6M May 23 23:29 rg_MarkDuplicates_BAM_output.bai
-rw-r--r-- 1 root root 711M May 23 23:21 rg_MarkDuplicates_BAM_output.bam
```
### Call Variants Using GATK4 HaplotypeCaller
Now that all the files have been pre-processed we are ready to run GATK4 HaplotypeCaller (Note: We use the `basename` command to strip the `.bam` suffix from the file name to rename it as a vcf):
We need to uncompress the `ucsc.hg19.fasta.gz` and associated files before passing them to GATK:
```bash
gunzip /home/hg19/*
```
We can now run GATK4 to call variants on our BAM files:
```bash
gatk HaplotypeCaller -R /home/hg19/ucsc.hg19.fasta -I /home/md_rg_$bam_file --genotyping-mode DISCOVERY --standard-min-confidence-threshold-for-calling 30 -O /home/md_rg_$(basename $bam_file .bam).vcf
```

