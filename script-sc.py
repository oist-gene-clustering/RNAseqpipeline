######### HEADER

#usage: 
#       input:  (1) FASTQ file (raw sequence of the sequencer)
#               (2) FASTA file (reference genome data)
#               (3) GTF file (gene annotation)
#       output: count matrix of gene expression
#IMPORTANT NOTE: The GTF, FASTA and FASTQ files should have the same name
#SECOND IMPORTANT NOTE: About PE reads, please name them [FILE]-1.fq [FILE]-2.fq
#THIRD IMPORTANT NOTE: From https://github.com/griffithlab/rnaseq_tutorial/wiki/Annotation ('cause this is valid for SC-RNA-seq as well)
#"In order for your RNA-seq analysis to work, the chromosome names in your .gtf file must 
#match those in your reference genome (i.e. your reference genome fasta file). If you get a StringTie result 
#where all transcripts have an expression value of 0, you may have overlooked this. Unfortunately, Ensembl, NCBI,
#and UCSC can not agree on how to name the chromosomes in many species, so this problem may come up often. You can 
#avoid this by getting a complete reference genome and gene annotation package from the same source (e.g., Ensembl) 
#to maintain consistency. Your annotations must correspond to the same reference genome build as your reference genome 
#fasta file. e.g., both correspond to UCSC human build 'hg38', NCBI human build 'GRCh38', etc. Even if both your reference 
#genome and annotations are from UCSC or Ensembl they could still correspond to different versions of that genome."
#FOURTH IMPORTANT NOTE: There is no normalization for the data here. If the data is not normalized, results will be probably f***ed up.

#This script is for ONLY ONE-file analysis
#(can be easily extended to several-file processing though: see script-sc_multiple.sh)

#THIS SCRIPT IS FOR ANALYSIS OF SINGLE-CELL RNA-SEQ RESULTS
#see script.py for bulk RNA sequencing

#This thing is quite awesome:
#https://github.com/seandavi/awesome-single-cell
#Comparison with the bulk RNA pipeline:
#Stegle, O., Teichmann, S. A., & Marioni, J. C. (2015). Computational and analytical challenges in single-cell transcriptomics. 
#Nature Reviews Genetics, 16(3), 133-145.

#######################################################################################################################################

######### INITIALIZING SOFTWARES

#On Sango
#$ ssh [user]@login.oist.jp (if outside of the OIST)
#$ ssh [user]@sango.oist.jp (if inside of the OIST)
#$ module load fastqc/0.11.4
#$ module load Trimmomatic/0.33
#$ module load bowtie2/2.2.6
#$ module load samtools/1.2
#$ module load bamtools/2.3.0
#$ module load cufflinks/2.2.1
#$ module load kallisto/v0.43.0
#$ module load bowtie/1.1.0

#######################################################################################################################################

######### GENERAL PARAMETERS

### Do not forget to update this according to where files are installed
### The following values are the ones needed to run the program on Sango

### INPUT/OUTPUT FILES
chariot = "\n"
gtfFolder = "annotations/"
fastaFolder = "ref-genomes/"
fastqFolder = "seq/"
qcFolder = "qc-files/"
bamFolder = "bam-files/"
toolsFolder = ""
indexFolder = "index/"
samFolder = "sam-files/"
cuffFolder = "cufflinks-files/"
kallFolder = "kallisto-files/"
stringTieOutputFolder = "stringTieOutput-files/"
qualityCheckFolder = "qualityCheck-files/"
rsemRefFolder = "rsemRef/"
geneResultsFolder = "geneResults/"

### SOFTWARE ACCESS
fastqcFolder = toolsFolder + ""
trimmomaticFolder = toolsFolder + "/apps/free/Trimmomatic/0.33/lib/"
bowtie2Folder = toolsFolder + ""
samtoolsFolder = toolsFolder + ""
bamtoolsFolder = toolsFolder + ""
cufflinksFolder = toolsFolder + ""
kallistoFolder = toolsFolder + ""
stringTieFolder = toolsFolder + "stringtie/"
htseqFolder = toolsFolder + "HTSeq-0.7.2/scripts/"
sinQCFolder = toolsFolder + "SINQC_v1.5/"
rsemFolder = toolsFolder + "RSEM-1.3.0/"

#######################################################################################################################################

######### QUALITY CONTROL PARAMETERS (FastQC)
#See software manual:
#https://biof-edu.colorado.edu/videos/dowell-short-read-class/day-4/fastqc-manual

#######################################################################################################################################

######### POST-ALIGNMENT QUALITY CONTROL PARAMETERS (SinQC)

# SinQC is designed not only for detecting technical artifacts, but also for generating general 
# quality related information. This parameter is to define how many genes can be detected with minimal TPM cutoff.
tpmCutoff = 1
# "To define gene expression outliers (GEOs), SinQC calculates a list of Spearman rank correlations
# of a given cell to the rest of the cells (one-to-others), as well as pairwise correlations after 
# removing that cell. A one-sided Wilcoxon signed-rank test is calculated to assess whether the one-
# to-others is significantly lower than overall pairwise correlations. This parameter is the p-value 
# cutoff to define GEOs."
pValueCutoffSpearman = 0.001
# Pearson correlation
pValueCutoffPearson = 0.001
# "The options are AND or OR
# If AND, SinQC will define GEOs as cells with both p-values being significant
# If OR, SinQC will define GEOs as cells with either p-value being significant"
corTag = 'AND'
# "The maximal false positive allowed. SinQC estimates the quality "bottom lines" 
# by requiring that at least 1- FDR fraction of MPCs should pass both of the MQS and WCQS cutoffs. 
# Then SinQC applies them to GEOs to determine technical artifacts. The default is 0.05."
maxFPR = 0.05

#######################################################################################################################################

######### TRIMMING PARAMETERS (Trimmomatic)
## Note for sc Analysis: you have actually (as far as I understand it) to trim the data 
#for the same artifacts as bulk RNA data + some new things explained in the following parts
#You can actually estimate the parameters to choose by looking at FASTQC output
#(if the quality is good enough >20, there is no need for trimming)
#and knowing which sequencing pipeline you used:
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4728800/
#http://rnaseq.uoregon.edu/#analysis
#More options can be added: see software manual:
#http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

### The phred is a quality score for the identification of the nucleobases
### The bigger it is, the best is your file
### https://en.wikipedia.org/wiki/Phred_quality_score
### There are two options for this parameter: "phred64" and "phred33"
### in fact, it depends on the sequencing pipeline you used to get the files
phred = "phred64"
#THREADST: Number of threads used for this step
threadst = 8
#ADAPTER: The adapter used for the sequencing pipeline (should be a string)
adapter = None
#LEADING: Cut bases off the start of a read, if below a threshold quality
leadingQ = None
#TRAILING: Cut bases off the end of a read, if below a threshold quality
trailingQ = None
#SLIDINGWINDOW: Performs a sliding window trimming approach
#It starts scanning at the 5' end and clips the read once the average 
#quality within the window falls below a threshold 
#(which is the number provided by the parameter)
sliding1Q = 4
sliding2Q = 15
#CROP: Cut the read to a specified length by removing bases from the end
crop = None
#HEADCROP: Cut the specified number of bases from the start of the read
headcrop = None
#MINLEN: Drop the read if it is below a specified length
minLen = None
#AVGQUAL: Drop the read if the average quality is below the specified level
avgQ = None
#TOPHRED33: Convert quality scores to Phred-33
tophred33 = False
#TOPHRED64: Convert quality scores to Phred-64
tophred64 = False

#######################################################################################################################################

######### ALIGNMENT PARAMETERS (Bowtie2)
#More options can be added: see software manual:
#http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

#SKIP: Skipping the first <skip> reads/pairs in the input
skip = None
#THREADS: Number of threads used to run the alignment step
threads = 8
#UPTO: Stop after first <upto> reads/pairs in the input
upto = None
#TRIM5: Trim <trim5> bases from 5' left end of reads
trim5 = 0
#TRIM3: Trim <trim3> bases from 3' right end of reads
trim3 = 0
#PHRED33: If you used PHRED 33 (else in PHRED 64)
phred33 = True
#N: Max number of mismatches in seed alignment: can be 0 or 1
n = 0
#L: Length of seed substrings: 32 > l >  3
l = 22
#LOCAL: When set to True, ends might be soft clipped; otherwise, entire read must align: no clipping
local = False
#MP: Penalty for mismatches: the lower it is, the worse the quality will be
mp = 6
#NP: Penalty for non-A/C/G/Ts in read/ref
np = 1

#######################################################################################################################################

######### PARAMETERS (Samtools & Bamtools)
#See software manual
#http://samtools.sourceforge.net (Samtools)
#https://github.com/pezmaster31/bamtools/ (Bamtools)

#######################################################################################################################################

######### POST-ALIGNMENT ANALYSIS PARAMETERS (Cufflinks)
#More options can be added: see software manual:
#http://cole-trapnell-lab.github.io/cufflinks/manual/

#BIASCORRECTION: Well, self-explanatory
biasCorrection = False
#LIBPREPARATION: The library preparation used for input reads
libPreparation = None
#NORMMETHOD: Method to normalize the library sizes
normMethod = None
#FRAGLENMEAN: Average fragment length (unpaired/SE reads only)
fragLenMean = 200
#FRAGLENSTDDEV: Fragment length std deviation (unpaired/SE reads only)
fragLenStdDev = 80
#COMPATIBLEHITSNORM: Count hits compatible with reference RNAs only
compatibleHitsNorm = False
#TOTALHITSNORM: Count all hits for normalization
totalHitsNorm = True
#MINISOFORMFRACTION: Suppress transcripts below this abundance level
minIsoformFraction = 0.10
#PREMRNAFRACTIOM: Suppress intra-intronic transcripts below this level 
premRNAFraction = 0.15
#MAXINTRONLENGTH: Ignore alignments with gaps longer than this
maxIntronLength = 300000

#######################################################################################################################################

######### PARAMETERS (Kallisto)
#More options can be added: see software manual:
#https://pachterlab.github.io/kallisto/manual

#KALLISTO PARAMETERS
#BIAS: Do you want bias correction?
bias = False
#SINGLE: Are you dealing only with SE reads?
single = True
#KMER (not a CambHodian one...): Precise the k-mer (odd) length
kmer = 31
#FRAGLENMEANK: Average fragment length (unpaired/SE reads only) should be the same as for Cufflinks...
fragLenMeanK = 200
#FRAGLENSTDDEVK: Fragment length std deviation (unpaired/SE reads only) should be the same as for Cufflinks...
fragLenStdDevK = 80

#######################################################################################################################################

######### PARAMETERS (HTSeq-Count & StringTie)
#See software manual
#http://www-huber.embl.de/HTSeq/doc/count.html (HTSC)
#http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual (StringTie)
    
#######################################################################################################################################

######### PARAMETERS (RSEM)
#See software manual
#http://deweylab.biostat.wisc.edu/rsem/README.html

#PHRED: Quality score format (either 33 or 64, or 0 if you don't want quality scores)
phred = 33
#THREADSR: Number of threads
threadsr = 8
#OUTPUTGEN: Do you want to output a genome BAM file?
outputGen = True
#FRAGLENMEANR: Average fragment length (unpaired/SE reads only) should be the same as for Cufflinks/Kallisto...
fragLenMeanR = 200
#FRAGLENSTDDEVR: Fragment length std deviation (unpaired/SE reads only) should be the same as for Cufflinks/Kallisto...
fragLenStdDevR = 80
#Is input file SAM?
samInput = False
#Fasta for upstream files?
isfasta = False
    
#######################################################################################################################################

######### USEFUL PROCEDURES

### GET PATH TO FILE
def getFolder(filetype):
    if filetype == "gtf":
        return gtfFolder
    if filetype == "fasta":
        return fastaFolder
    if filetype == "fastq":
        return fastqFolder
    if filetype == "bam":
        return bamFolder
    if filetype == "sam":
        return samFolder
    if filetype == "qcf":
        return qualityCheckFolder
    if filetype == "geneResults":
        return geneResultsFolder
    else:
        raise ValueError

def getExtension(filetype):
    if filetype == "gtf":
        return ".gtf"
    if filetype == "fasta":
        return ".fa"
    if filetype == "fastq":
        return ".fq"
    if filetype == "bam":
        return ".bam"
    if filetype == "sam":
        return ".sam"
    if filetype == "qcf":
        return ".qcf"
    if filetype == "geneResults":
        return ".genes.results"
    else:
        raise ValueError

def getPath(filename, filetype):
    return (getFolder(filetype) + filename + getExtension(filetype))

### SANITIZING INPUTS IS LOVE, SANITIZING IS LIFE
def sanitize(i):
    i = [i for i in i.split(chariot) if i][0]
    i = [i for i in i.split('.0') if i][0]
    return i

### MISCELLANENOUS
def multList(l1, l2):
    res = []
    assert len(l1) == len(l2)
    for i in range(len(l1)):
        res += [l1[i]*l2[i]]
    return res

def compVar(l1, l2, meanScore):
    res = 0
    assert len(l1) == len(l2)
    for i in range(len(l1)):
        for j in range(l2[i]):
            res += (l1[i] - meanScore)*(l1[i] - meanScore)
    return res

def getQC(line):
    return [int(sanitize(i)) for i in line.split('\t') if i]

def getInteger(line):
    from re import search
    match = search(r'\d+', line)
    assert match
    return int(match.group(0))

def findLine(pattern, text):
    from re import search
    match = search(pattern, text)
    return (not(match==None))

### CALLING A BASH PROCESS

# Without looking at the output
def runCmd(cmd):
    from subprocess import call
    return call(cmd, shell=True)

# Getting the output
def getResFromCmd(cmd):
    from subprocess import check_output
    return check_output(cmd, shell=True)

#######################################################################################################################################

######### WORKFLOW STEPS

##############
### QUALITY CONTROL TEST

def qualityControl(filename, filetype):
    cmd = fastqcFolder + "fastqc --extract -q -o " + qcFolder + " " + getPath(filename, filetype)
    runCmd(cmd)
    fastqcName = qcFolder + filename + "_fastqc"
    runCmd("mv " + fastqcName + "/fastqc_data.txt " + fastqcName + "_data.txt")
    fastqcName += "_data.txt"
    qualityScores = []
    counts = []
    with open(fastqcName, "r") as f:
        #Line indice
        #Just printing some info about the quality control
        i = 1
        for line in f:
            if i == 7:
                nbSeq = getInteger(line)
            if i == 8:
                badSeq = getInteger(line)
            if findLine(r'Per sequence quality scores', line):
                break
            i += 1
        for line in f:
            if findLine(r'END_MODULE', line):
                break
            if findLine(r'#Quality', line):
                continue
            t = getQC(line)
            qualityScores += [t[0]]
            counts += [t[1]]
    assert sum(counts)
    tmp = multList(qualityScores, counts)
    meanScore = sum([i for i in tmp])/sum(counts)
    from math import sqrt
    stdScore = sqrt(compVar(qualityScores, counts, meanScore)/sum(counts))          
    print ("----- QUALITY CONTROL for file " + filename)
    print ("# Poor-quality sequences: " + str(badSeq) + " out of " + str(nbSeq))
    print ("Average per sequence quality score: " + str(meanScore))
    if meanScore > 28:
        print "(Good quality)"
    elif meanScore > 20:
        print "(Correct quality)"
    else:
        print "(Bad quality: should be trimmed)"
    print ("Standard deviation of per sequence quality score: " + str(int(stdScore)))
    print ("-----")
    
def qualityControlSplit(filename, filetype, mode):
    if mode == "SE":
        qualityControl(filename, filetype)
    else:
        qualityControl(filename + "-1", filetype)
        qualityControl(filename + "-2", filetype)
        
##############
### TRIMMING STEP

def adapterTrim(filename, mode):
    #Multithread
    cmd = "java -jar " + trimmomaticFolder + "trimmomatic-0.33.jar"
    #Single-end mode
    if mode == "SE":
        cmd += " SE " + "-threads " + threadst + " -" + phred + " " + getPath(filename, "fastq")
        cmd += " " + getPath(filename + "_processed", "fastq")
    elif mode == "PE":
        #First read of the paired couple
        cmd += " PE " + "-threads " + threadst + " -" + phred + " " + getPath(filename + "-1", "fastq")
        #Second read of the paired couple
        cmd += " " + getPath(filename + "-2", "fastq")
        #Paired output 1
        cmd += " " + getPath(filename + "_processed-p1", "fastq")
        #Unpaired output 1
        cmd += " " + getPath(filename + "_processed-up1", "fastq")
        #Paired output 2
        cmd += " " + getPath(filename + "_processed-p2", "fastq")
        #Unpaired output 2
        cmd += " " + getPath(filename + "_processed-up2", "fastq")
    else:
        raise ValueError
    if adapter:
        cmd += " ILLUMINACLIP:" + adapter
    if leadingQ:
        cmd += " LEADING:" + str(leadingQ)
    if trailingQ:
        cmd += " TRAILING:" + str(trailingQ)
    if sliding1Q and sliding2Q:
        cmd += " SLIDINGWINDOW:" + str(sliding1Q) + ":" + str(sliding2Q)
    if crop:
        cmd += " CROP:" + str(crop)
    if headcrop:
        cmd += " HEADCROP:" + str(headcrop)
    if minLen:
        cmd += " MINLEN:" + str(minLen)
    if avgQ:
        cmd += " AVGQUAL:" + str(avgQ)
    if tophred33 and not tophred64:
        cmd += " TOPHRED33"
    if tophred64 and not tophred33:
        cmd += " TOPHRED64"
    print ("----- ADAPTER TRIMMING")
    runCmd(cmd) 
    print ("-----")
    
def adapterTrimSplit(filename, mode, trim):
    if trim:
        adapterTrim(filename, mode)
    else:
        if mode == "SE":
            runCmd("cp " + getPath(filename, "fastq") + " " + getPath(filename + "_processed", "fastq"))
        else:
            runCmd("cp " + getPath(filename + "-1", "fastq") + " " + getPath(filename + "_processed-p1", "fastq"))
            runCmd("cp " + getPath(filename + "-2", "fastq") + " " + getPath(filename + "_processed-p2", "fastq"))
    
##############
### ALIGNEMENT STEP
    
def align(filename, mode):
    print ("----- ALIGNMENT")
    folder = indexFolder + filename + "-index"
    #Indexing files
    print "----- INDEXING FILES"
    runCmd("mkdir " + folder)
    runCmd(bowtie2Folder + "bowtie2-build " + getPath(filename, "fasta") + " " + folder)
    #see later for other options
    #Mapping
    print "-----"
    print "----- MAPPING"
    #bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
    cmd = bowtie2Folder + "bowtie2 -q " 
    if threads:
        cmd += "-p " + str(threads) + " "
    if skip:
        cmd += "-s " + str(skip) + " "
    if upto:
        cmd += "-u " + str(upto) + " "
    if trim5:
        cmd += "-5 " + str(trim5) + " "
    if trim3:
        cmd += "-3 " + str(trim3) + " "
    if not phred33:
        cmd += "--phred64 "
    if n == 1 or n == 0:
        cmd += "-N " + str(n) + " "
    if 3 < l and l < 32:
        cmd += "-L " + str(l) + " "
    if local:
        cmd += "--local --ma 2 "
    if mp:
        cmd += "--mp " + str(mp) + " "
    if np:
        cmd += "--np " + str(np) +  " "
    cmd += "-x " + folder + " " 
    if mode == "SE":
        cmd += "-U " + getPath(filename + "_processed", "fastq") + " -S " + getPath(filename, "sam")
    else:
        cmd += "-1 " + getPath(filename + "_processed-p1", "fastq") + " "
        cmd += "-2 " + getPath(filename + "_processed-p2", "fastq") + " -S " + getPath(filename, "sam")
    print cmd
    runCmd(cmd)
    print "-----"

##############
### ASSEMBLY STEP

def assemblyAndExpr(filename, mode):
    print ("----- ASSEMBLY STEP")
    cmd = cufflinksFolder + "cufflinks -q "
    if biasCorrection:
        cmd += "-b " + getPath(filename, "fasta") + " "
    if libPreparation:
        cmd += "--library-type " + libPreparation + " "
    if normMethod:
        cmd += "--library-norm-method " + normMethod + " "
    if compatibleHitsNorm:
        cmd += "--compatible-hits-norm "
    if totalHitsNorm:
        cmd += "--total-hits-norm "
    if minIsoformFraction:
        cmd += "--min-isoform-fraction " + str(minIsoformFraction) + " "
    if premRNAFraction:
        cmd += "--pre-mrna-fraction " + str(premRNAFraction) + " "
    if maxIntronLength:
        cmd += "--max-intron-length " + str(maxIntronLength) + " "
    if (mode == "SE" and fragLenMean):
        cmd += "--frag-len-mean " + str(fragLenMean) + " "
    if (mode == "SE" and fragLenStdDev):
        cmd += "--frag-len-std-dev " + str(fragLenStdDev) + " "
    cmd += getPath(filename, "bam") + " > " + cuffFolder + "computation.txt"
    runCmd(cmd)
    print "-----"

##############
### KALLISTO STEP

def referenceFromData(filename, mode, clean = False):
    print ("----- KALLISTO STEP")
    cleanedFilename = getPath(filename + "_cleaned", "fasta")
    if clean:
        # Cleans header from FASTA files
        cmd = "cat " + getPath(filename, "fasta") + " | perl -ne 'if ($_ =~/^\>\d+\s+\w+\s+"
        cmd += "(ERCC\S+)[\+\-]/){print \">$1\n\"}elsif($_ =~ /\d+\s+(ENST\d+)/){print \">$1\n\"}"
        cmd += "else{print $_}' > " + cleanedFilename
        runCmd(cmd)
    else:
        cmd = "cp " + getPath(filename, "fasta") + " " + cleanedFilename
        runCmd(cmd)
    print "----- INDEXING FILES"
    cmd = kallistoFolder + "kallisto index -i " + indexFolder + filename + "_kindex -k " + str(kmer) + " " + getPath(filename, "fasta")
    runCmd(cmd)
    print ("----- QUANTIFICATION ALGORITHM")
    cmd = kallistoFolder + "kallisto quant -i " + indexFolder + filename + "_kindex -o " + kallFolder + " --plaintext "
    if bias:
        cmd += "--bias "
    if single:
        cmd += "--single "
    if mode == "SE":
        cmd += getPath(filename, "fastq")
        # Generate abundance estimate
        runCmd(cmd)
    else:
        cmd1 = cmd + getPath(filename + "-1", "fastq")
        cmd2 = cmd + getPath(filename + "-2", "fastq")
        runCmd(cmd1)
        runCmd(cmd2)
    print "-----"
    
##############
### COUNT READS STEP (HTSEQ-COUNT)

#cf. code

##############
### TRANSCRIPTOME QUANTIFICATION (RSEM)

def transcriptQuant(filename, mode, do = True): 
    #Build references
    cmd = rsemFolder + "rsem-prepare-reference "
    #cmd += "--gtf " + getPath(filename, "gtf") + " \ "
    #cmd += "--bowtie2 \ "
    #if bowtie2Folder:
    #     cmd += "--bowtie2-path " + bowtie2Folder + " \ "
    cmd += getPath(filename, "fasta") + " " + rsemRefFolder + filename
    runCmd(cmd)
    #see option --alignments in calculate-expression
    #Skipping part for de novo assembled transcriptomes
    #For Gene Expression
    #If a SAM file has been computed before (to dodge the alignement part in RSEM)
    #You get a BAM file as output
    ## TODO
    #if samInput:
        #cmd = rsemFolder + "convert-sam-for-rsem " + getPath(filename, "sam") + " -o " + getPath(filename + "_for_rsem", "bam")
        #runCmd(cmd)
        #Then use --sam option
    #If no SAM file has been computed before
    cmd = rsemFolder + "rsem-calculate-expression -q "
    #cmd += "--bowtie2 \ "
    #if bowtie2Folder:
    #    cmd += "--bowtie2-path ." + bowtie2Folder + " \ "
    #if threadsr:
    #    cmd += "-p " + str(threadsr) + " \ "
    #if phred:
    #    cmd += "--phred" + str(phred) + "-quals \ "
    #if outputGen:
    #    cmd += "--output-genome-bam \ "
    #else:
    #    cmd += "--no-bam-output \ "
    #if fragLenMeanR:
    #    cmd += "--fragment-length-mean " + str(fragLenMeanR) + " \ "
    #if fragLenStdDevR:
    #    cmd += "--fragment-length-sd " + str(fragLenStdDevR) + " \ "
    #if isfasta:
    #    cmd += "--no-qualities \ "
    if mode == "PE":
        cmd += "--paired-end " + getPath(filename + "_1", "fastq") + " " + getPath(filename + "_2", "fastq") + " "
    else:
        cmd += getPath(filename, "fastq") + " "
    cmd += rsemRefFolder + filename + " "
    cmd += geneResultsFolder + filename
    print cmd
    runCmd(cmd)

##############
### POST-ALIGNMENT (SINQC)

def qualityCheck(filename):
    cmd = "python " + sinQCFolder + "SINQC.py -RSEM " + geneResultsFolder + " -SEQ " 
    cmd += fastqFolder + " -o " + qualityCheckFolder
    if tpmCutoff:
        cmd += " -TPMCutoff " + str(tpmCutoff)
    if pValueCutoffSpearman:
        cmd += " -PValueCutoff--Distinct--Spearman " + str(pValueCutoffSpearman)
    elif pValueCutoffPearson:
        cmd += " -PValueCutoff--Distinct--Pearson " + str(pValueCutoffPearson)
    if corTag == "OR":
        cmd += " -CorTag 'OR'"
    if maxFPR:
        cmd += " -Max_FPR " + str(maxFPR)
    runCmd(cmd)

##############
### NORMALIZATION (if you know the data)

# Based on:
#Comparison with the bulk RNA pipeline:
#Stegle, O., Teichmann, S. A., & Marioni, J. C. (2015). Computational and analytical challenges in single-cell transcriptomics. 
#Nature Reviews Genetics, 16(3), 133-145.

#You should normalize the counts according to the mRNA content in each cell (or using UMI, cf. article)
#Plus, you should discard cells with degraded mRNA, that is having a too low percentage of mapped read to the reference genome
#Or having a too low percentage of mapped reads to the spike-in molecules (cf. article)
        
##############
### HOUSEKEEPING STEP after the execution of the program
### IN CASE OF SINGLE FILE ANALYSIS do = True
### I think it would be too costly to create/delete the folders at every step for
### a multiple-file analysis

def initializeFolders(do = False):
    if do:
        runCmd("mkdir " + qcFolder)
        runCmd("mkdir " + bamFolder)
        runCmd("mkdir " + samFolder)
        runCmd("mkdir " + indexFolder)
        runCmd("mkdir " + cuffFolder)
        runCmd("mkdir " + kallFolder)
        runCmd("mkdir " + qualityCheckFolder)
        runCmd("mkdir " + rsemRefFolder)
        runCmd("mkdir " + geneResultsFolder)
        
def cleanFiles(filename, fastqcName, do = False):
    if do:
        runCmd("rm " + fastqcName + ".zip")
        runCmd("rm " + fastqcName + ".html")
        runCmd("rm -r " + fastqcName)
        runCmd("rm " + getPath(filename, "sam"))
        runCmd("rm " + getPath(filename, "bam"))
        runCmd("rm -r " + getPath(filename, "index"))
        runCmd("rm " + getPath(filename + "_assembled", "gtf"))
        runCmd("rm " + getPath(filename + "_processed", "fastq"))
        runCmd("rm " + getPath(filename + "_cleaned", "fasta"))
        runCmd("rm -r " + qcFolder)
        runCmd("rm -r " + bamFolder)
        runCmd("rm -r " + samFolder)
        runCmd("rm -r " + indexFolder)
        runCmd("rm -r " + cuffFolder)
        runCmd("rm -r " + kallFolder)
        runCmd("rm -r " + rsemRefFolder)
        runCmd("rm -r " + qualityCheckFolder)
        runCmd("rm -r " + geneResultsFolder)
    runCmd("rm script.pyc")
        
#######################################################################################################################################

######### sc RNA-seq ANALYSIS PIPELINE

def main(filename, mode, trim):
    initializeFolders()
    if not filename or not mode:
        raise ValueError
        
    ### SANITIZING THE INPUT FILES
    
    ### QUALITY CONTROL (with FASTQC)
    #qualityControlSplit(filename, "fastq", mode)
    #fastqcName = qcFolder + filename + "_fastqc"

    ### ADAPTER TRIM (with TRIMMOMATIC)
    #adapterTrimSplit(filename, mode, trim)
    
    ### ALIGNMENT (with BOWTIE2)
    #align(filename, mode)
    
    ### READ COUNTING (with BAMTOOLS and SAMTOOLS)
    # Convert SAM file to BAM file
    #runCmd(samtoolsFolder + "samtools view -b -o " + getPath(filename, "bam") + " " + getPath(filename, "sam"))
    # Sort BAM file
    #runCmd(bamtoolsFolder + "bamtools sort -in " + getPath(filename, "bam") + " -out " + getPath(filename, "bam"))
    # Visualize BAM file
    #runCmd(samtoolsFolder + "samtools tview " + getPath(filename, "bam") + " " + getPath(filename, "fasta"))
    #print ("----- COUNT ALIGNMENTS IN BAM FILE")
    #runCmd(bamtoolsFolder + "bamtools count -in " + getPath(filename, "bam") + " > " + bamFolder + "counts.txt")
    print ("-----")
    print ("----- COVERAGE STATS")
    #runCmd(bamtoolsFolder + "bamtools coverage -in " + getPath(filename, "bam") + " -out " + getPath(filename + "coverage", "bam"))
    print ("-----")
    print ("----- OTHER STATS")
    #runCmd(bamtoolsFolder + "bamtools stats -in " + getPath(filename, "bam") + " > " + bamFolder + "stats.txt")
    print ("-----")
    
    ### Basically, it seems that the general analysis pipeline for sc-seq data is mostly the same
    ### as the one for RNA-seq data, minus some more trimming and normalization
    
    ### TRANSCRIPTOME QUANTIFICATION (with RSEM)
    transcriptQuant(filename, mode, False)
    
    ### POST-ALIGNMENT QUALITY CONTROL (with SINQC)
    qualityCheck(filename)
    
    ### TRANSCRIPTOME ASSEMBLY and DIFFERENTIAL EXPRESSION ESTIMATION (with CUFFLINKS)
    #assemblyAndExpr(filename, mode)
    
    ### REFERENCING FREE TRANSCRIPT EXPRESSION ABUNDANCE ESTIMATION FROM DATA (with KALLISTO)
    #referenceFromData(filename, mode)
    
    #Old stuff with StringTie and HTSeq-count (actually quite useful to get a count matrix)
    #Help for StringTie: http://ccb.jhu.edu/software/stringtie/index.shtml
    #Help for HTSeq-Count: http://www-huber.embl.de/HTSeq/doc/count.html
    cmd = stringTieFolder + "stringtie -G " + getPath(filename, "gtf") + " -e -B -o " + stringTieOutputFolder + filename + "_assembled.gtf"
    cmd += " " + getPath(filename, "bam")
    runCmd(cmd)
    #Should give you the count matrix
    cmd = htseqFolder + "htseq-count --format bam --order pos --mode intersection-strict " + getPath(filename, "sam") + " " 
    cmd += stringTieOutputFolder + filename + "_assembled.gtf"
    runCmd(cmd)
    
    #cleanFiles(filename, fastqcName)
    
#######################################################################################################################################

######### DEBUGGING LINE, SHOULD BE DISMISSED
main("homosap", "PE", "N")
######### to use this pipeline, open python shell + call function "main" with the good arguments