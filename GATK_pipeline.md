---


---

<h1 id="running-gatk-in-docker-for-advanced-clinical-bioinformatics-2019">Running GATK in Docker for Advanced Clinical Bioinformatics 2019</h1>
<p><strong>These notes assume you are running a Linux environment.  The commands can be copied direct from this file but may need amending depending on your file names and locations.</strong></p>
<h3 id="setting-up--running-docker">Setting Up &amp; Running Docker</h3>
<p>First make a folder to hold your files on your local machine:</p>
<pre class=" language-bash"><code class="prism  language-bash"><span class="token function">mkdir</span> AdvClinBio2019
<span class="token function">cd</span> AdvClinBio2019/
</code></pre>
<p>Pull down the gatk4 docker image (for instructions on installing docker please see here???):</p>
<pre class=" language-bash"><code class="prism  language-bash">docker pull broadinstitute/gatk
</code></pre>
<p>Use <code>wget</code> to download the GATK bundle from the Broad Institue website (this will be several gigabytes of additional data required for the analysis and may take some time to download):</p>
<pre class=" language-bash"><code class="prism  language-bash"><span class="token function">wget</span> -pr ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/
</code></pre>
<p>Download the aligned bam file from Galaxy and save it to your working directory (this should have also been processed with Picard’s <code>MarkDuplicates</code> tool before entering the GATK pipeline).  My file is called ‘MarkDuplicates_on_data_16__MarkDuplicates_BAM_output.bam’</p>
<p>We will map our current working directory to our docker image so that we have access to these files from within the docker container:</p>
<pre class=" language-bash"><code class="prism  language-bash"><span class="token function">sudo</span> docker run -v <span class="token variable"><span class="token variable">`</span><span class="token function">pwd</span><span class="token variable">`</span></span>:/home  -i -t broadinstitute/gatk
</code></pre>
<p>To map a local drive to the container we use the format <code>-v /local_directory:/container_directory</code>.  The command `pwd` basically returns our current working directory, rather than us having to specify it directly, and maps it to the <code>/home</code> directory in our docker container.  This saves us having to copy all the required files manually to the container.  The three flags we use in this command are:</p>
<ul>
<li><code>--volume , -v</code>	Bind mount a volume</li>
<li><code>--interactive , -i</code>		Keep STDIN open even if not attached</li>
<li><code>--tty , -t</code>		Allocate a pseudo-TTY</li>
</ul>
<p>Your prompt should now change to indicate that you you have entered an interactive session with the docker container:<br>
<code>(gatk) root@d6c83ac31f5a:/gatk#</code><br>
We can get our current directory by using the <code>pwd</code> command:</p>
<pre class=" language-bash"><code class="prism  language-bash"><span class="token function">pwd</span>
/gatk
</code></pre>
<p>We are in the <code>/gatk</code> folder (this information could also be provided by the prompt depending on how your system is configured).  Lets check whether our mounted files are available to use within the container,  remember that they are in our <code>/home</code> directory.</p>
<pre class=" language-bash"><code class="prism  language-bash"><span class="token function">ls</span> -lh /home
total 525M
-rw-rw-r-- 1 1000 1000 525M May 23 12:22 MarkDuplicates_BAM_output.bam
drwxr-xr-x 2 1000 1000  16K Apr 16 09:21 hg19
</code></pre>
<p>If all your files appear in your <code>/home</code> directory then you are ready to proceed with the analysis. If at any point you need to exit the docker container you can simply type <code>exit</code> in the terminal.</p>
<h3 id="preparing-the-bam-file-for-analysis">Preparing the BAM file for analysis</h3>
<p>We require these six files from the bundle:<br>
/home/hg19/ucsc.hg19.fasta.gz<br>
/home/hg19/ucsc.hg19.fasta.fai.gz<br>
/home/hg19/ucsc.hg19.dict.gz<br>
/home/hg19/hapmap_3.3.hg19.sites.vcf.gz<br>
/home/hg19/hapmap_3.3.hg19.sites.vcf.idx.gz</p>
<p>Before proceeding will first assign the name of our BAM file to a variable in bash so that the code below will work for any file through the magic of parameter expansion:</p>
<pre class=" language-bash"><code class="prism  language-bash">bam_file<span class="token operator">=</span><span class="token string">"MarkDuplicates_BAM_output.bam"</span>
</code></pre>
<p>Substitute <code>"MarkDuplicates_BAM_output.bam"</code> for the name of your BAM file. Every time I reference <code>$bam_file</code> on the bash commandline it will now be interpreted as <code>MarkDuplicates_BAM_output.bam</code>, or whatever filename you allocated to this variable.</p>
<p>Lets check our BAM file by running the <code>CountReads</code> command below:</p>
<pre class=" language-bash"><code class="prism  language-bash">gatk CountReads -I /home/<span class="token variable">$bam_file</span>
</code></pre>
<p>You should see output similar to below:</p>
<pre class=" language-bash"><code class="prism  language-bash">Running:
    java -Dsamjdk.use_async_io_read_samtools<span class="token operator">=</span>false -Dsamjdk.use_async_io_write_samtools<span class="token operator">=</span>true -Dsamjdk.use_async_io_write_tribble<span class="token operator">=</span>false -Dsamjdk.compression_level<span class="token operator">=</span>2 -jar /gatk/gatk-package-4.1.2.0-local.jar CountReads -I /home/MarkDuplicates_BAM_output.bam
22:40:15.286 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/gatk/gatk-package-4.1.2.0-local.jar<span class="token operator">!</span>/com/intel/gkl/native/libgkl_compression.so
May 23, 2019 10:40:17 PM shaded.cloud_nio.com.google.auth.oauth2.ComputeEngineCredentials runningOnComputeEngine
INFO: Failed to detect whether we are running on Google Compute Engine.
22:40:17.089 INFO  CountReads - ------------------------------------------------------------
22:40:17.089 INFO  CountReads - The Genome Analysis Toolkit <span class="token punctuation">(</span>GATK<span class="token punctuation">)</span> v4.1.2.0
22:40:17.089 INFO  CountReads - For support and documentation go to https://software.broadinstitute.org/gatk/
22:40:17.090 INFO  CountReads - Executing as root@d6c83ac31f5a on Linux v4.8.0-53-generic amd64
22:40:17.090 INFO  CountReads - Java runtime: OpenJDK 64-Bit Server VM v1.8.0_191-8u191-b12-0ubuntu0.16.04.1-b12
22:40:17.090 INFO  CountReads - Start Date/Time: May 23, 2019 10:40:15 PM UTC
22:40:17.090 INFO  CountReads - ------------------------------------------------------------
22:40:17.090 INFO  CountReads - ------------------------------------------------------------
22:40:17.091 INFO  CountReads - HTSJDK Version: 2.19.0
22:40:17.091 INFO  CountReads - Picard Version: 2.19.0
22:40:17.091 INFO  CountReads - HTSJDK Defaults.COMPRESSION_LEVEL <span class="token keyword">:</span> 2
22:40:17.091 INFO  CountReads - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS <span class="token keyword">:</span> <span class="token boolean">false</span>
22:40:17.091 INFO  CountReads - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS <span class="token keyword">:</span> <span class="token boolean">true</span>
22:40:17.091 INFO  CountReads - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE <span class="token keyword">:</span> <span class="token boolean">false</span>
22:40:17.091 INFO  CountReads - Deflater: IntelDeflater
22:40:17.091 INFO  CountReads - Inflater: IntelInflater
22:40:17.091 INFO  CountReads - GCS max retries/reopens: 20
22:40:17.092 INFO  CountReads - Requester pays: disabled
22:40:17.092 INFO  CountReads - Initializing engine
22:40:17.568 INFO  CountReads - Done initializing engine
22:40:17.568 INFO  ProgressMeter - Starting traversal
22:40:17.568 INFO  ProgressMeter -        Current Locus  Elapsed Minutes       Reads Processed     Reads/Minute
22:40:32.307 INFO  CountReads - 10584740 read<span class="token punctuation">(</span>s<span class="token punctuation">)</span> filtered by: WellformedReadFilter 

22:40:32.312 INFO  ProgressMeter -             unmapped              0.2                     0              0.0
22:40:32.312 INFO  ProgressMeter - Traversal complete. Processed 0 total reads <span class="token keyword">in</span> 0.2 minutes.
22:40:32.312 INFO  CountReads - Shutting down engine
<span class="token punctuation">[</span>May 23, 2019 10:40:32 PM UTC<span class="token punctuation">]</span> org.broadinstitute.hellbender.tools.CountReads done. Elapsed time: 0.28 minutes.
Runtime.totalMemory<span class="token punctuation">(</span><span class="token punctuation">)</span><span class="token operator">=</span>683147264
Tool returned:
0
</code></pre>
<p>From this verbose output we can see that my BAM file hase 10584740 reads.</p>
<p>Before we can use this BAM with HaplotypeCaller we need to ensure that all the metadata which GATK expects using Picard tool’s <code>AddOrReplaceReadGroups</code>.<br>
Required Arguments:</p>
<ul>
<li><code>--INPUT,-I:String</code>          Input file (BAM or SAM or a GA4GH url).  Required.</li>
<li><code>--OUTPUT,-O:File</code>              Output file (BAM or SAM).  Required.</li>
<li><code>--RGLB,-LB:String</code>             Read-Group library  Required.</li>
<li><code>--RGPL,-PL:String</code>             Read-Group platform (e.g. illumina, solid)  Required.</li>
<li><code>--RGPU,-PU:String</code>             Read-Group platform unit (eg. run barcode)  Required.</li>
<li><code>--RGSM,-SM:String</code>             Read-Group sample name  Required.<br>
The exact values for some of this metadata depends on your dataset, below I use some reasonable names for the cardiac panel I used:</li>
</ul>
<pre class=" language-bash"><code class="prism  language-bash">gatk AddOrReplaceReadGroups -I /home/<span class="token variable">$bam_file</span> -O /home/rg_<span class="token variable">$bam_file</span> -LB cardio -PL illumina -PU H8LHPADXX -SM cardio2 --VALIDATION_STRINGENCY<span class="token operator">=</span>LENIENT
</code></pre>
<p>To check that this process has created a valid BAM file we use <code>ValidateSamFile</code>:</p>
<pre class=" language-bash"><code class="prism  language-bash">gatk ValidateSamFile -I /home/rg_<span class="token variable">$bam_file</span> -MODE<span class="token operator">=</span>SUMMARY
</code></pre>
<p>You should see the tool output someting like:</p>
<pre class=" language-bash"><code class="prism  language-bash">Using GATK jar /gatk/gatk-package-4.1.2.0-local.jar
Running:
    java -Dsamjdk.use_async_io_read_samtools<span class="token operator">=</span>false -Dsamjdk.use_async_io_write_samtools<span class="token operator">=</span>true -Dsamjdk.use_async_io_write_tribble<span class="token operator">=</span>false -Dsamjdk.compression_level<span class="token operator">=</span>2 -jar /gatk/gatk-package-4.1.2.0-local.jar ValidateSamFile -I /home/rg_MarkDuplicates_BAM_output.bam --MODE<span class="token operator">=</span>SUMMARY
23:25:00.896 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/gatk/gatk-package-4.1.2.0-local.jar<span class="token operator">!</span>/com/intel/gkl/native/libgkl_compression.so
<span class="token punctuation">[</span>Thu May 23 23:25:00 UTC 2019<span class="token punctuation">]</span> ValidateSamFile  --INPUT /home/rg_MarkDuplicates_BAM_output.bam --MODE SUMMARY  --MAX_OUTPUT 100 --IGNORE_WARNINGS <span class="token boolean">false</span> --VALIDATE_INDEX <span class="token boolean">true</span> --INDEX_VALIDATION_STRINGENCY EXHAUSTIVE --IS_BISULFITE_SEQUENCED <span class="token boolean">false</span> --MAX_OPEN_TEMP_FILES 8000 --SKIP_MATE_VALIDATION <span class="token boolean">false</span> --VERBOSITY INFO --QUIET <span class="token boolean">false</span> --VALIDATION_STRINGENCY STRICT --COMPRESSION_LEVEL 2 --MAX_RECORDS_IN_RAM 500000 --CREATE_INDEX <span class="token boolean">false</span> --CREATE_MD5_FILE <span class="token boolean">false</span> --GA4GH_CLIENT_SECRETS client_secrets.json --help <span class="token boolean">false</span> --version <span class="token boolean">false</span> --showHidden <span class="token boolean">false</span> --USE_JDK_DEFLATER <span class="token boolean">false</span> --USE_JDK_INFLATER <span class="token boolean">false</span>
May 23, 2019 11:25:02 PM shaded.cloud_nio.com.google.auth.oauth2.ComputeEngineCredentials runningOnComputeEngine
INFO: Failed to detect whether we are running on Google Compute Engine.
<span class="token punctuation">[</span>Thu May 23 23:25:02 UTC 2019<span class="token punctuation">]</span> Executing as root@d6c83ac31f5a on Linux 4.8.0-53-generic amd64<span class="token punctuation">;</span> OpenJDK 64-Bit Server VM 1.8.0_191-8u191-b12-0ubuntu0.16.04.1-b12<span class="token punctuation">;</span> Deflater: Intel<span class="token punctuation">;</span> Inflater: Intel<span class="token punctuation">;</span> Provider GCS is available<span class="token punctuation">;</span> Picard version: Version:4.1.2.0
WARNING	2019-05-23 23:25:02	ValidateSamFile	NM validation cannot be performed without the reference. All other validations will still occur.
INFO	2019-05-23 23:25:43	SamFileValidator	Validated Read    10,000,000 records.  Elapsed time: 00:00:41s.  Time <span class="token keyword">for</span> last 10,000,000:   41s.  Last <span class="token function">read</span> position: chr20:34,170,844
No errors found
<span class="token punctuation">[</span>Thu May 23 23:25:46 UTC 2019<span class="token punctuation">]</span> picard.sam.ValidateSamFile done. Elapsed time: 0.76 minutes.
Runtime.totalMemory<span class="token punctuation">(</span><span class="token punctuation">)</span><span class="token operator">=</span>595066880
Tool returned:
0
</code></pre>
<p>The key output is <code>No errors found</code>. Because we now create an index for the BAM file (I believe the <code>MarkDuplicates</code> tool used below requires an index to function):</p>
<pre class=" language-bash"><code class="prism  language-bash">gatk BuildBamIndex -I /home/rg_<span class="token variable">$bam_file</span>
</code></pre>
<p>Lets check that worked:</p>
<pre class=" language-bash"><code class="prism  language-bash"><span class="token function">ls</span> -lh /home/
total 1.3G
-rw-rw-r-- 1 1000 1000 525M May 23 12:22 MarkDuplicates_BAM_output.bam
drwxrwxr-x 3 1000 1000 4.0K May 23 12:28 ftp.broadinstitute.org
drwxr-xr-x 2 1000 1000  16K Apr 16 09:21 hg19
-rw-r--r-- 1 root root 5.6M May 23 23:29 rg_MarkDuplicates_BAM_output.bai
-rw-r--r-- 1 root root 711M May 23 23:21 rg_MarkDuplicates_BAM_output.bam
</code></pre>
<p>We  can then use Picard’s <code>MarkDuplicates</code> tool to remove duplicate reads which represent the same genomic fragment so that we do not over estimate our coverage.</p>
<pre class=" language-bash"><code class="prism  language-bash">gatk MarkDuplicates -I /home/rg_<span class="token variable">$bam_file</span> -O /home/md_rg_<span class="token variable">$bam_file</span> -M /home/markup_metrics.txt
</code></pre>
<p>As we have changed the BAM again we need to create a new index for the BAM file:</p>
<pre class=" language-bash"><code class="prism  language-bash">gatk BuildBamIndex -I /home/md_rg_<span class="token variable">$bam_file</span>
</code></pre>
<p>Lets check that we have all the files we expect:</p>
<pre class=" language-bash"><code class="prism  language-bash"><span class="token function">ls</span> -lh /home/
total 2.0G
-rw-rw-r-- 1 1000 1000 525M May 23 12:22 MarkDuplicates_BAM_output.bam
drwxrwxr-x 3 1000 1000 4.0K May 23 12:28 ftp.broadinstitute.org
drwxr-xr-x 2 1000 1000  16K Apr 16 09:21 hg19
-rw-r--r-- 1 root root 5.6M May 23 23:38 md_rg_MarkDuplicates_BAM_output.bai
-rw-r--r-- 1 root root 712M May 23 23:36 md_rg_MarkDuplicates_BAM_output.bam
-rw-r--r-- 1 root root 5.6M May 23 23:29 rg_MarkDuplicates_BAM_output.bai
-rw-r--r-- 1 root root 711M May 23 23:21 rg_MarkDuplicates_BAM_output.bam
</code></pre>
<h3 id="call-variants-using-gatk4-haplotypecaller">Call Variants Using GATK4 HaplotypeCaller</h3>
<p>Now that all the files have been pre-processed we are ready to run GATK4 HaplotypeCaller (Note: We use the <code>basename</code> command to strip the <code>.bam</code> suffix from the file name to rename it as a vcf):<br>
We need to uncompress the <code>ucsc.hg19.fasta.gz</code> and associated files before passing them to GATK:</p>
<pre class=" language-bash"><code class="prism  language-bash">gunzip /home/hg19/*
</code></pre>
<p>We can now run GATK4 to call variants on our BAM files:</p>
<pre class=" language-bash"><code class="prism  language-bash">gatk HaplotypeCaller -R /home/hg19/ucsc.hg19.fasta -I /home/md_rg_<span class="token variable">$bam_file</span> --genotyping-mode DISCOVERY --standard-min-confidence-threshold-for-calling 30 -O /home/md_rg_<span class="token variable"><span class="token variable">$(</span><span class="token function">basename</span> $bam_file .bam<span class="token variable">)</span></span>.vcf
</code></pre>

