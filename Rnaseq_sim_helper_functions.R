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


#Function for retrieving mapped read counts per gene
GetMappedCounts <- function(prefix,
                            dir=".",
                            type="STAR",
                            gff="~/Sequences/ref_genomes/tomato/ITAG2.4_Chromo2.5/ITAG2.4_gene_models.gff3")
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
  if(type=="STAR-featureCounts") {
    mapped_counts_list <- featureCounts(
      files=file.path(dir,paste(prefix,"Aligned.sortedByCoord.out.bam",sep="")),
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
  stop(paste("Unknown result type:",type))}


#Function to compare known and mapped counts

CompareCounts <- function(known.counts,mapped.counts, title= NULL, plot = TRUE, chrom.separate = TRUE, correlation=TRUE, return.merged.table=TRUE) {
  count.comparison <- merge(known.counts,mapped.counts,by.x="lyc.id",by.y="ID",suffixes=c(".known",".mapped"))
  count.comparison$chrom <- substr(count.comparison$lyc.id,6,7)
  read.correlation <- cor(count.comparison$count.known,count.comparison$count.mapped)
  if (correlation) print(read.correlation)
  if(plot & !chrom.separate) {
    plot(count.comparison$count.known,count.comparison$count.mapped,
         main=paste(title,"\nCorrelation:", 
                    round(read.correlation,3)))
  }
  if(plot & chrom.separate) {
    print(qplot(x=count.known,y=count.mapped,data=count.comparison) + 
            facet_wrap(~ chrom) + 
            ggtitle(paste(title,"\nCorrelation:", round(read.correlation,3))))
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
