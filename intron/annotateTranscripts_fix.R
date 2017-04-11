annotateTranscripts_fix <-function(txdb, annotationPackage=NULL, by=c("tx","gene"),codingOnly=FALSE,verbose=TRUE,requireAnnotation=FALSE, mappingInfo = NULL){

    if( class(txdb)!="TxDb") stop("txdb must be of class TxDb")

    if(is.null(annotationPackage)){
        organism <- organism(txdb)
        organism <- strsplit(organism," ")[[1]]
        organism <- paste0(substr(organism[1],1,1),
                           tolower(substr(organism[2],1,1)))
        annotationPackage <-  paste("org",organism,"eg.db",sep=".")
        if(verbose)
            message(sprintf("No annotationPackage supplied. Trying %s.",annotationPackage))
    }  

    if(!require(annotationPackage,character.only=TRUE)){
        if(requireAnnotation){
            stop("Can't load ",annotationPackage,".\nMake sure library is installed.\nAnd make sure the organism argument follows the convention here http://www.bioconductor.org/packages/release/data/annotation/.\nFor example for human use Hs")}else{
                message("Could not load ",annotationPackage,". Will continue without annotation")
                annotationPackage <- NULL
            }
    }
    
##################################################
### Get TSS and TSE
#######################################################

    if(verbose) message("Getting TSS and TSE.")
    
    by <- match.arg(by)
        
    if(by=="tx"){
        tt <- transcriptsBy(txdb, by="gene")
        seqinfo <- seqinfo(tt)
        ## the transcript names
        LL <- tt@unlistData@elementMetadata@listData
        txx <- unlist(tt)$tx_name
    } else{
        tt <- genes(txdb)
        seqinfo <- seqinfo(tt)
        txx <- names(tt)
    }
    
    
    chr <- unlist(CharacterList(seqnames(tt)))
    strand <- unlist(CharacterList(strand(tt)))
    
    RR <- ranges(tt)
    geneid <- names(RR)
    ## starts / ends (left/right endpoints):
    TSS <- unlist(start(RR))
    TSE <- unlist(end(RR))
    
    
######################################################
### Get CSS and CSE
#######################################################
    
    if(verbose) message("Getting CSS and CSE.")

    cds <- cdsBy(txdb, by=by, use.names=(by=="tx"))
    with_coding <- which(txx %in% names(cds))
    
    ## take only the left- and right-most endpts:
    cds_ranges <- range(cds[txx[with_coding]])
    
    ## starts and ends (really left / right endpts):
    CSS <- CSE <- integer(length(txx))
    CSS[with_coding] <- unlist(start(cds_ranges))
    CSE[with_coding] <- unlist(end(cds_ranges))
    CSS[-with_coding] <- NA
    CSE[-with_coding] <- NA
    
    ##################################################
    ### Get Exons
    #######################################################

    if(verbose) message("Getting exons.")

    ee <- exonsBy(txdb, by=by, use.names=(by=="tx"))
    ## exon groups of interest in transcript order,
    ## and sorted "left to right" within each group
    Exons <- reduce(ranges(ee)[txx])
    Nexons <- elementNROWS(Exons)

    ##now annotate genes
    if(by=="tx") TT <- elementNROWS(tt) else TT <- rep(1,length(tt))
    Geneid <- Rle(geneid, TT)

    Tx <- txx

    if(!is.null(annotationPackage)){
        if(verbose) message("Annotating genes.")
            
        if(!is.null(mappingInfo)) {
            stopifnot(all(c('column', 'keytype', 'multiVals') %in% names(mappingInfo)))
            geneid <- mapIds(annotationPackage, keys = geneid, column= mappingInfo$column, keytype = mappingInfo$keytype, multiVals = mappingInfo$multiVals)
        }
        

        ##Annotate transcrtipt
        ##cant use select cause map is not 1-1 
        map <- get(gsub("\\.db", "SYMBOL", annotationPackage))
        which <- mappedkeys(map)
        symbols <- sapply(as.list(map[which]), paste, collapse=" ")
        genes <- symbols[geneid]	
        
        map <- get(gsub("\\.db", "REFSEQ", annotationPackage))
        which <- mappedkeys(map)
        symbols <- sapply(as.list(map[which]), paste, collapse=" ")
        refseq <- symbols[geneid]	# a handful of NA
        
        Gene <- Rle(genes, TT)
        Refseq <- Rle(refseq, TT)
    } else {
        Gene <- Rle(NA, sum(TT))
        Refseq <- Rle(NA, sum(TT))
    }

    transcripts=GRanges(ranges=IRanges(start=TSS, end=TSE),
        seqnames=chr, strand=strand,
        seqinfo=seqinfo,
        CSS, CSE, Tx, Geneid, Gene, Refseq, Nexons, Exons)
    if(codingOnly) transcripts <- transcripts[!is.na(values(transcripts)$CSS),]
    attributes(transcripts)$description <- "annotatedTranscripts"
    return(transcripts)
}
