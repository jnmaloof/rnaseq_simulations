# Helper functions for RNAseq simulation project.



#Function for running the read mapping program STAR

RunStar <- function(fasta, prefix, param.file = "../STAR.params.whitney.1", args = "", dir =".", twopassMode = "Basic", n=12) {
  #This function runs the external program STAR to map reads against the genome reference
  old.dir <- getwd()
  if (!dir.exists(dir)) dir.create(dir,recursive = TRUE) 
  setwd(dir)
  system(paste("STAR", 
               "--parametersFiles", param.file,
               "--readFilesIn", fasta,
               ifelse(grepl("\\.gz$",fasta),"--readFilesCommand 'gunzip -c'",""),
               "--outFileNamePrefix", prefix,
               "--twopassMode", twopassMode,
               "--runThreadN",n,
               args))
  system(paste("samtools index ",prefix,"Aligned.sortedByCoord.out.bam",sep=""))
  setwd(old.dir)
}

RunKallisto <- function(fastq, index, dir, prefix, threads=2) {
  # old.dir <- getwd()
  # if (!dir.exists(dir)) dir.create(dir,recursive = TRUE) 
  # setwd(dir)
  system(paste("kallisto quant --single --plaintext -l 250 -s 50",
               "-t", threads,
               "-i", index, 
               "-o ", file.path(dir,prefix), 
               file.path(dir,fastq))) #, 
  # setwd(old.dir)
}

RunBowtie2 <- function(fastq, index, dir, prefix, threads=2) {
  bampath <- file.path(dir,paste(prefix,"bam",sep="."))
  system(paste("bowtie2 --local",
               "--threads", threads,
               "-x",index,
               "-U",file.path(dir,fastq),
               "| samtools view -Sb - >", bampath))
}

RunBWAaln <- function(fastq, index, dir, prefix, threads=2) {
  bampath <- file.path(dir,paste(prefix,"bam",sep="."))
  system(paste("bwa aln -t",threads,
               index,
               file.path(dir,fastq),
               ">",file.path(dir,paste(prefix,"sai",sep="."))))
  system(paste("bwa samse",index,
               file.path(dir,paste(prefix,"sai",sep=".")),
               file.path(dir,fastq),
               "| samtools view -Sb - | samtools sort - > ",bampath))
  system(paste("samtools index",bampath))
}

RunStampy <- function(fasta, index, dir, prefix, threads=2) {
  bampath <- file.path(dir,paste(prefix,sep="."))
  system(paste("stampy.py -g",index,
               "-h",index,
               "-t",threads,
               "--inputformat=fasta",
               "-M", file.path(dir,fasta),
               "| samtools view -Sb - >", bampath))
}

#Function for retrieving mapped read counts per gene
GetMappedCounts <- function(prefix=NA,
                            dir=".",
                            type="STAR",
                            gff="~/Sequences/ref_genomes/tomato/ITAG2.4_Chromo2.5/ITAG2.4_gene_models.gff3",
                            bam=NA)
{
  #Get mapped counts from alignment output files
  if(type=="STAR") {
    file <- file.path(dir,paste(prefix,"ReadsPerGene.out.tab",sep=""))
    mapped_counts <- read.delim(file,header = FALSE, stringsAsFactors = FALSE)
    mapped_counts <- mapped_counts[,1:2] # 3 and 4 are directional counts which we don't care about
    colnames(mapped_counts) <- c("ID","count")
    mapped_counts$ID <- sub("gene:","",mapped_counts$ID,fixed=TRUE)
    return(mapped_counts)
  }
  if(grepl("featureCounts",type)) { 
    bampath <- ifelse(type=="STAR-featureCounts",
                      file.path(dir,paste(prefix,"Aligned.sortedByCoord.out.bam",sep="")),
                      file.path(dir,bam)
    )
    mapped_counts_list <- featureCounts(
      files=bampath,
      annot.ext=getAbsolutePath(gff,expandTilde = TRUE),
      isGTFAnnotationFile = TRUE,
      GTF.attrType = "Parent",
      useMetaFeatures = TRUE,
      allowMultiOverlap = TRUE,
      nthreads=2)
    mapped_counts <- data.frame(
      ID=rownames(mapped_counts_list$counts),
      count=mapped_counts_list$counts[,1])
    mapped_counts$ID <- substr(sub("mRNA:","",mapped_counts$ID,fixed=TRUE),1,16)
    mapped_counts <- with(mapped_counts_list,{
      colnames(stat) <- c("ID","count")
      rbind(stat,mapped_counts)
    })
    return(mapped_counts)
  }
  if(type=="transcripts") { #ie reads are mapped to a cDNA-type reference rather than genomic.
    if(is.na(bam)) stop(paste("bamfile needed"))
    bampath <- file.path(dir,bam)
    mapped_counts <- read.table(pipe(paste("samtools idxstats",bampath)))[,-2]
    colnames(mapped_counts) <- c("ID","count","unmapped")
    mapped_counts <- mapped_counts[grepl("Sopen.+Solyc",mapped_counts$ID),]
    mapped_counts$ID <- regmatches(
      mapped_counts$ID,regexpr(
        "Solyc[[:digit:]]{2}g[[:digit:]]{6}\\.[[:digit:]]",
        mapped_counts$ID))
    return(mapped_counts)
  }
  if(type=="kallisto-transcripts") {
    countspath <- file.path(dir,"abundance.tsv")
    mapped_counts <- read.delim(countspath,stringsAsFactors = FALSE)[,c(1,4)]
    head(mapped_counts)
    colnames(mapped_counts) <- c("ID","count")
    mapped_counts$ID <- regmatches(
      mapped_counts$ID,regexpr(
        "Solyc[[:digit:]]{2}g[[:digit:]]{6}\\.[[:digit:]]",
        mapped_counts$ID))
    return(mapped_counts)
  }
  stop(paste("Unknown result type:",type))
}


#Function to compare known and mapped counts

CompareCounts <- function(known.counts,mapped.counts=NA, title= NULL, plot = TRUE, chrom.separate = TRUE, correlation=TRUE, return.merged.table=!is.na(mapped.counts[[1]][1])) {
  if(!is.na(mapped.counts[[1]][1])) {
    count.comparison <- merge(known.counts,mapped.counts,by.x="lyc.id",by.y="ID",suffixes=c(".known",".mapped"))
    count.comparison$chrom <- substr(count.comparison$lyc.id,6,7)
  } else {
    count.comparison <- known.counts
  }
  read.correlation.pearson <- round(cor(count.comparison$count.known,count.comparison$count.mapped),3)
  read.correlation.spearman <- round(cor(count.comparison$count.known,count.comparison$count.mapped,method = "spearman"),3)
  if (correlation) print(paste("Pearson correlation:", read.correlation.pearson, "  Spearman correlation:",read.correlation.spearman))
  if(plot & !chrom.separate) {
    plot(count.comparison$count.known,count.comparison$count.mapped,
         main=paste(title,"\nPearson:", 
                    read.correlation.pearson,
                    "Spearman",
                    read.correlation.spearman))
  }
  if(plot & chrom.separate) {
    print(qplot(x=count.known,y=count.mapped,data=count.comparison) + 
            facet_wrap(~ chrom) + 
            ggtitle(paste(title,"\nCorrelation:", 
                          read.correlation.pearson,
                          "Spearman",
                          read.correlation.spearman)))
  }
  if(return.merged.table) {
    count.comparison$count.diff.known.mapped <- count.comparison$count.known - count.comparison$count.mapped
    count.comparison$count.ratio.known.mapped <- count.comparison$count.known / count.comparison$count.mapped
    count.comparison$pen.gene.length <- pen.gene.length[count.comparison$id]
    count.comparison$lyc.gene.length <- lyc.gene.length[count.comparison$id]
    count.comparison$gene.length.diff.pen.lyc <- count.comparison$pen.gene.length - count.comparison$lyc.gene.length
    count.comparison$gene.length.ratio.pen.lyc <- count.comparison$pen.gene.length / count.comparison$lyc.gene.length
    
    return(count.comparison)
  }
}

PSL2Granges <- function(file) {
  psl.data <- read.delim(file,header=FALSE,skip=5)
  psl.header <- make.names(apply(read.delim(file,header=FALSE,skip=2,nrow=2),2,paste,collapse=" "))
  psl.header <- gsub("[ \\.]+",".",psl.header)
  psl.header <- make.names(gsub("(\\.$)|(^\\.)","",psl.header))
  colnames(psl.data) <- psl.header
  psl.data$pen.ID <- substr(psl.data$Q.name,1,16)
  psl.data$lyc.ID <- substr(psl.data$Q.name,18,33)
  with(psl.data,
       GRanges(seqnames = T.name,
               ranges = IRanges(start=T.start,end=T.end,names=Q.name),
               strand=strand,
               matches=match,
               mismatches=mis.match,
               Q.size=Q.size,
               pen.ID=pen.ID,
               lyc.ID=lyc.ID)
  )
}

sem <- function(x, na.rm=TRUE) {
  if(na.rm) x <- na.omit(x)
  sd(x) / sqrt(length(x))
}

  