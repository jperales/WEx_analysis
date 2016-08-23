#!/usr/bin/env Rscript

# Inspired from : http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html
require("optparse")

if(length(strsplit(system2("whereis","bedtools",stdout=TRUE),split=" ")[[1]])==1) {
  stop("ERROR #1 : bedtools is not installed in your computer.");
}

option_list = list(
  make_option(c("--bed"), type="character", default=NULL, 
              help="BED file with the regions or intervals from the capture design", metavar="character"),
  make_option(c("--bams"), type="character", default=NULL, 
              help="BAM files separated by commas", metavar="character"),
  make_option(c("--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

bam.files <- strsplit(opt$bams,split=",")[[1]]

fls <- vector("character") # The BedTools output files
for(bam in bam.files) {
  prefix <- gsub(".bam","",basename(bam))
  

  file.tmp <- tempfile(paste0(prefix,".hist_allTR"),fileext = ".tsv");
  cmd <- paste0("bedtools coverage -hist"," -abam ",bam," -b ",opt$bed,
                " | grep ^all 1> ",
                # " /tmp/",prefix,".hist_allTR.txt");
                file.tmp)
  cat(paste0("[sample : ",prefix,"]"," Calculating coverage by bedtools.\n"),file=stdout())
  cat(paste0("[command]\t",cmd,"\n"),file=stdout())
  try(system(cmd))
  
  if(file.exists(file.tmp)) {
    fls <- c(fls,file.tmp)
    names(fls)[length(fls)] <- prefix
  }
}

# Get a list of the bedtools output fls you'd like to read in
# labs <- gsub("^(.*)(_recalibrated|whatever)$","\\1",names(fls))
labs <- gsub("^resultsRS_(.*)_recalibrated$","\\1",names(fls))

# Create lists to hold coverage and cumulative coverage for each alignment,
# and read the data into these lists.
cov <- vector("list")
cov_cumul <- vector("list")
for (i in 1:length(fls)) {
  cov[[i]] <- read.table(fls[i])
  cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
}

# Pick some colors
# Ugly:
# cols <- 1:length(cov)
# Prettier:
# ?colorRampPalette
# display.brewer.all()
library(RColorBrewer)
cols <- brewer.pal(ifelse(length(cov)<3,3,length(cov)), "Dark2")

# Save the graph to a file
cat(paste0("[Figure]","Generating figure.\n"),file=stdout())
png(opt$out, h=600*3, w=800*3, res=280)

# Create plot area, but do not plot anything. Add gridlines and axis labels.
plot(cov[[1]][2:401, 2], cov_cumul[[1]][1:400],
     type='n', xlab="Depth",
     ylab="Fraction of capture target bases \u2265 depth", ylim=c(0,1.0), main="Target Region Coverage",
     las=1)
abline(v = 20, col = "gray60")
abline(v = 50, col = "gray60")
abline(v = 80, col = "gray60")
abline(v = 100, col = "gray60")
abline(h = 0.50, col = "gray60")
abline(h = 0.90, col = "gray60")
axis(1, at=c(20,50,80), labels=c(20,50,80))
axis(2, at=c(0.90), labels=c(0.90),las=1)
axis(2, at=c(0.50), labels=c(0.50),las=1)

# Actually plot the data for each of the alignments (stored in the lists).
for (i in 1:length(cov)) points(cov[[i]][2:401, 2], cov_cumul[[i]][1:400], type='l', lwd=3, col=cols[i])

# Add a legend using the nice sample labeles rather than the full filenames.
legend("topright", legend=labs, col=cols, lty=1, lwd=4)

dev.off()

# awk -F"\t" 'BEGIN{c=0;l=0;b=0;print "Avg_cov\tPerc_basesCov";}{c=c+$6;l=l+1;if($6!=0){b=b+1;}}END{print c/l,"\t",b/l*100;}' resultsRS_normal_recalibrated.bam.hist.TR.txt 