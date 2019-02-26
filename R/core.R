# Qiwne Hu 2018
# core functions - scATAC-seq pipeline

library(Rsamtools)
library(Rcpp)
library(dplyr)
sourceCpp("../src/ccount.ntags.barcode.cpp")
dyn.load("../src/cdensum.so")
dyn.load("../src/peaks.so")

# Function to read bam files and generate vectors of read positions and barcode informaiton
read.bam.tags <- function(filename,read.tag.names=F,fix.chromosome.names=F) {
  if(!is.element("Rsamtools", installed.packages()[, 1])) {
    stop("Rsamtools Bioconductor package is now required for BAM file support. Please install")
  }
  
  ww <- c("flag","rname","pos","isize","strand","mapq","qwidth"); if(read.tag.names) { ww <- c(ww,"qname") };
  bam <- Rsamtools::scanBam(filename,param=Rsamtools::ScanBamParam(what=ww,tag="XC",
                                                                   flag=Rsamtools::scanBamFlag(isUnmappedQuery=FALSE)))[[1]]
  # index barcodes
  barcode.table <- index.barcodes(unique(bam$tag$XC))
  barcode.index <- dplyr::inner_join(data.frame(barcode = bam$tag$XC), barcode.table)
  bam$bindex <- barcode.index$ID
  
  if(is.null(bam$pos) || length(bam$pos)==0) { return(list(tags=c(),quality=c())) }
  strm <- as.integer(bam$strand=="+")
  if(any(bitwAnd(bam$flag,0x1))) {
    # paired-end data
    # use only positive strand mappings
    posvi <- which(strm==1);
    rl <- list(tags=tapply(posvi,bam$rname[posvi],function(ii) as.numeric(bam$pos[ii])),
               flen=tapply(posvi,bam$rname[posvi],function(ii) as.numeric(abs(bam$isize[ii]))))
    rl <- c(rl,list(quality=tapply(1:length(bam$pos[posvi]),bam$rname[posvi],function(ii) bam$mapq[ii])))
    
    # add barcode information
    rl <- c(rl, list(barcode=tapply(1:length(bam$pos[posvi]),bam$rname[posvi],function(ii) bam$tag$XC[ii])))
    rl <- c(rl, list(bindex=tapply(1:length(bam$pos[posvi]),bam$rname[posvi],function(ii) bam$bindex[ii])))
    
    # alternatively, handle reads with NA isize (unpaired?) just like single-ended reads
    #pos <- tapply(1:length(bam$pos),bam$rname,function(ii) ifelse(is.na(bam$isize[ii]), bam$pos[ii]*strm[ii]  - (1-strm[ii])*(bam$pos[ii]+bam$qwidth[ii]), strm[ii]*bam$pos[ii] - (1-strm[ii])*(bam$pos[ii]+bam$isize[ii])))
  } else {
    rl <- list(tags=tapply(1:length(bam$pos),bam$rname,function(ii) bam$pos[ii]*strm[ii]  - (1-strm[ii])*(bam$pos[ii]+bam$qwidth[ii])))
    rl <- c(rl,list(quality=tapply(1:length(bam$pos),bam$rname,function(ii) bam$mapq[ii])))
    
    # add barcode information
    rl <- c(rl, list(barcode=tapply(1:length(bam$pos),bam$rname,function(ii) bam$tag$XC[ii])))
    rl <- c(rl, list(bindex=tapply(1:length(bam$pos[posvi]),bam$rname[posvi],function(ii) bam$bindex[ii])))
  }
  
  if(read.tag.names) {
    rl <- c(rl,list(names=tapply(1:length(bam$pos),bam$rname,function(ii) bam$qname[ii])))
  }
  if(fix.chromosome.names) {
    # remove ".fa"
    names(rl) <- gsub("\\.fa","",names(rl))
  }
  return(list(rl = rl, index = barcode.table))
}

#function to index barcodes
index.barcodes <- function(barcodes){
  #Input: a vector of barcode
  #Output: a dataframe list barcodes and there correpodent idexes
  df <- data.frame(barcode = barcodes)
  df <- df %>% dplyr::mutate(ID = group_indices_(df, .dots=c("barcode")))
  return(df)
}

# calculate cumulative density based on sum of scaled gaussian curves
densum <- function(vin,bw=5,dw=3,match.wt.f=NULL,return.x=T,from=min(vin),to=max(vin),step=1,new.code=T){
  # (by Michael Tolstorukov)
  #
  # vin - input vector; bw -- standard deviation, dw-gaussina cutoff in stdev; dout - output "density")
  # output - if return.x=F vector of cumulative density values corresponding to integer positions described by range(vin)
  # output - if return.x=T a data structure with $x and $y corresponding to the cumulative density
  # optional match.wt.f is a function that will return weights for a tag vector
  
  # construct vector of unique tags and their counts
  tc <- table(vin[vin>=from & vin<=to])
  pos <- as.numeric(names(tc)); storage.mode(pos) <- "double"
  tc <- as.numeric(tc); storage.mode(tc) <- "double"
  n <- length(pos)
  # weight counts
  if(!is.null(match.wt.f)) {
    tc <- tc*match.wt.f(pos);
  }
  
  rng <- c(from,to);
  if(rng[1]<0) { stop("range extends into negative values") }
  if(range(pos)[1]<0) { stop("position vector contains negative values") }
  
  storage.mode(n) <- storage.mode(rng) <- storage.mode(bw) <- storage.mode(dw) <- storage.mode(step) <- "integer";
  
  spos <- rng[1]; storage.mode(spos) <- "double";
  
  dlength <- floor((rng[2] - rng[1])/step) + 1; # length of output array
  if(dlength<1) { stop("zero data range") }
  if(new.code) {
    storage.mode(step) <- storage.mode(dlength) <- storage.mode(bw) <- storage.mode(dw) <-"integer";
    dout <- .Call("ccdensum",pos,tc,spos,bw,dw,dlength,step);
  } else {
    stop("Please set new.code=T to use the new ccdensum function. The old cdensum is deprecated")
  }
  
  
  if(return.x) {
    return(list(x=c(rng[1],rng[1]+step*(dlength-1)),y=dout,step=step))
  } else {
    return(dout)
  }
}

# find positions of all local maxima above some threshold
peaks.c <- function(x,thr=min(x),max.span=1) {
  # by Peter Kharchenko
  storage.mode(x) <- storage.mode(thr) <- "double";
  storage.mode(max.span) <- "integer";
  results <- .Call("find_peaks",x,thr,max.span);
  return(results);
}

# Get peak locations and its correspondent smoothed density values
t.get.density.peaks <- function(signal.tags, bandwidth=150, step=round(bandwidth/3), span=3, thr=3) {
  # Modified from Peter Kharchenko's function
  rng <- range(signal.tags)
  ds <- densum(signal.tags,bw=bandwidth,from=rng[1],to=rng[2],return.x=T,step=step);
  pi <- peaks.c(ds$y,thr,span)
  #return(cbind(x=seq(ds$x[1],ds$x[2],by=step)[pi],y=ds$y[pi]))
  return(cbind(seq(ds$x[1],ds$x[2],by=step)[pi]))
}

# cpp function to generate count matrix
cppFunction('NumericVector cwindow_n_tags_barcodes(NumericVector pos, NumericVector barcode, int nbarcode,
                  NumericVector wpos, double window_half_size) {
  int nw = wpos.size();
  int n = pos.size();
  double whs = window_half_size;
  
  // output count matrix
  NumericMatrix count(nw, nbarcode);
  
  // current array start/end indecies
  int cs = 0; int ce = 0;

  for(int i = 0; i < nw; i++){
    // increment ce to add windows that are overlapped
    double ep = wpos[i] + whs;
    
    while(ep > pos[ce] & ce<n) { ce++ ; } 
    // increment cs to drop windows that have already ended
    
    double sp = wpos[i] - whs;
    while(pos[cs] < sp) { cs++ ; }

    for(int j=cs; j<ce; j++) { 
     // printf("%d", j);
      count(i, barcode[j] - 1) += 1;
    }
  }
  return count;
 }')
