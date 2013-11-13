## read.GenBank.R (2012-02-17)

##   Read DNA Sequences from GenBank via Internet

## Copyright 2002-2012 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

## Read only a specific gene. Added by Tim Lucas 13/02/13
## If argument 'gene' is given a string, only the sequence for that gene will be downloaded.
## If argument CDS is true and a gene is named only the CDS region will be downloaded
## All samples must contain the gene/CDS region asked for otherwise an error will be thrown.

read.GenBankGene <-
    function(access.nb, seq.names = access.nb, species.names = TRUE, gene=FALSE,
             gene.names = FALSE, as.character = FALSE, CDS=FALSE)
{
    if(!is.logical(CDS)) stop("CDS must be TRUE/FALSE")
    N <- length(access.nb)
    ## If there are more than 400 sequences, we need to break down the
    ## requests, otherwise there is a segmentation fault.
    nrequest <- N %/% 400 + as.logical(N %% 400)
    X <- character(0)
    for (i in 1:nrequest) {
        a <- (i - 1) * 400 + 1
        b <- 400 * i
        if (i == nrequest) b <- N
        URL <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                     paste(access.nb[a:b], collapse = ","),
                     "&rettype=gb&retmode=text", sep = "")
        X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
    }
    FI <- grep("^ {0,}ORIGIN", X) + 1
    LA <- which(X == "//") - 1
    obj <- vector("list", N)
    for (i in 1:N) {
        ## remove all spaces and digits
        tmp <- gsub("[[:digit:] ]", "", X[FI[i]:LA[i]])
        obj[[i]] <- unlist(strsplit(tmp, NULL))
    }
    names(obj) <- seq.names
    # Start TL contribution
    if ( is.character( gene ) ) { # Has the user requested a specific gene?
        posLines <- grep( paste('/gene="', gene,'"', sep='') , X )-1 # which lines have gene positions
        geneCDS <- if( CDS ) posLines[grep( 'CDS', X[posLines] )] else posLines[grep( 'gene', X[posLines] )] # Do we want the coding region or the entire gene?
        if( CDS & length( geneCDS )!=N ){
            stop( paste(access.nb[which(!(x[seq(1,N*2,2)]==x[seq(2,N*2,2)]-2))[1]], 'and possibly other samples do not define the CDS region' ) ) }
        if( !CDS & length( geneCDS )!=N ) stop( "Not all entries contain this gene" )
        pos <- lapply( strsplit( gsub( '[^0-9\\.]', '', X[geneCDS] ), '\\.\\.' ), as.numeric ) # extract positions.
        obj <- lapply( 1:length(obj), function(i, obj) { obj[[i]] <- obj[[i]][ pos[[i]][1]:pos[[i]][2] ]}, obj=obj ) # trim gene regions
    } #END TL contribution.
 
    if (!as.character) obj <- as.DNAbin(obj)
    if (species.names) {
        tmp <- character(N)
        sp <- grep("ORGANISM", X)
        for (i in 1:N)
            tmp[i] <- unlist(strsplit(X[sp[i]], " +ORGANISM +"))[2]
        attr(obj, "species") <- gsub(" ", "_", tmp)
    }
    if (gene.names) {
        tmp <- character(N)
        sp <- grep(" +gene +<", X)
        for (i in 1:N)
            tmp[i] <- unlist(strsplit(X[sp[i + 1L]], " +/gene=\""))[2]
        attr(obj, "gene") <- gsub("\"$", "", tmp)
    }
    obj
}
