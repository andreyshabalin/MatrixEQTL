# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

modelLINEAR = 117348L;
modelANOVA  = 47074L;
modelLINEAR_CROSS = 1113461L;

.seq = function(a,b){if(a<=b){a:b}else{NULL}}

formatnum = function(x, digits = 0){
    formatC(
        x = x,
        big.mark = ",",
        format = "f",
        digits = digits);
}

SlicedData = setRefClass("SlicedData",
    fields = list( 
        dataEnv = "environment",
        nSlices1 = "numeric",
        rowNameSlices = "list",
        columnNames = "character",
        fileDelimiter = "character",
        fileSkipColumns = "numeric",
        fileSkipRows = "numeric",
        fileSliceSize = "numeric",
        fileOmitCharacters = "character"
    ),
    methods = list(
        initialize = function( mat = NULL ){
            dataEnv <<- new.env(hash = TRUE, size = 29L);
            nSlices1 <<- 0L;
            if(!is.null(mat)){
                CreateFromMatrix(mat);
            }
            fileSliceSize <<- 1000;
            fileDelimiter <<- "\t";
            fileSkipColumns <<- 1L;
            fileSkipRows <<- 1L;
            fileOmitCharacters <<- "NA"
            return(invisible(.self));
        },
        CreateFromMatrix = function( mat ){
            stopifnot( is(mat, "matrix") );
            setSliceRaw( 1L ,mat );
            rns = rownames( mat, do.NULL = FALSE);
            rowNameSlices <<- list(rns);
            cns = colnames( mat, do.NULL = FALSE );

            columnNames <<- cns;
            return(invisible(.self));
        },
        getSlice = function(sl){
            value = get(paste0(sl), dataEnv);
            if( is.raw(value) ){
                storage.mode(value) = "double";
                value[value == 255] = NA;
            }
            return( value  )    
        },
        getSliceRaw = function(sl){
            return( get(paste0(sl), dataEnv) )    
        },
        setSliceRaw = function(sl, value){
            assign( paste0(sl), value, dataEnv )
            if( nSlices1 < sl ){
                nSlices1 <<- sl;
            }
        },
        setSlice = function(sl, value){
            if( length(value) > 0 ){
                if( all(as.integer(value) == value, na.rm = TRUE) ){
                    if( (min(value, na.rm = TRUE) >= 0 ) && 
                                (max(value, na.rm = TRUE) < 255) )
                    {
                        nv = value;
                        suppressWarnings({storage.mode(nv) = "raw"});
                        nv[ is.na(value)] = as.raw(255);
                        value = nv;
                    } else {
                        storage.mode(value) = "integer";
                    }
                }
            }
            setSliceRaw(sl, value);
        },    
        nSlices = function(){
            return( nSlices1 );
        },
        LoadFile = function(
                        filename,
                        skipRows = NULL,
                        skipColumns = NULL,
                        sliceSize = NULL,
                        omitCharacters = NULL,
                        delimiter = NULL,
                        rowNamesColumn = 1){
            if( !is.null(skipRows) )
                fileSkipRows <<- skipRows;
            if( !is.null(skipColumns) )
                fileSkipColumns <<- skipColumns;
            if( !is.null(omitCharacters) )
                fileOmitCharacters <<- omitCharacters;
            if( !is.null(sliceSize) )
                fileSliceSize <<- sliceSize;
            if( !is.null(delimiter) )
                fileDelimiter <<- delimiter;
            stopifnot(  (fileSkipColumns == 0) || 
                        (rowNamesColumn <= fileSkipColumns) );
            stopifnot(  (fileSkipColumns == 0) || 
                        (rowNamesColumn >= 1) );
    
            fid = file(
                        description = filename,
                        open = "rt",
                        blocking = FALSE,
                        raw = FALSE);
            # clean object if file is open
            Clear(); 
            lines = readLines(
                        con = fid,
                        n = max(fileSkipRows,1L),
                        ok = TRUE,
                        warn = TRUE);
            line1 = tail(lines,1);
            splt = strsplit(line1, split = fileDelimiter, fixed = TRUE);
            if( fileSkipRows > 0L ){
                columnNames <<- splt[[1]]; # [ -(1:fileSkipColumns) ];
            } else {
                seek(fid, 0)
            }        
            
            rm( lines, line1, splt );
            
            rowNameSlices <<- vector("list", 15);
    
            curSliceId = 0L;
            repeat
            {
                # preallocate names and data
                if(length(rowNameSlices) < curSliceId){
                    rowNameSlices[[2L*curSliceId]] <<- NULL;
                }
                curSliceId = curSliceId + 1L;
                
                # read sliceSize rows
                rowtag = vector("character",fileSliceSize);
                rowvals = vector("list",fileSliceSize);
                for(i in seq_len(fileSliceSize)){
                    temp = "";
                    if( fileSkipColumns > 0L ){
                        temp = scan(
                                file = fid,
                                what = character(),
                                n = fileSkipColumns,
                                quiet = TRUE,
                                sep = fileDelimiter);
                    }
                    rowtag[i] = temp[rowNamesColumn];#paste(temp,collapse=" ");
                    rowvals[[i]] = scan(
                                file = fid,
                                what = double(),
                                nlines = 1,
                                quiet = TRUE,
                                sep = fileDelimiter,
                                na.strings = fileOmitCharacters);
                    if( length(rowvals[[i]]) == 0L ){
                        if(i==1L){
                            rowtag = matrix(0, 0, 0);
                            rowvals = character(0);
                        } else {
                            rowtag  = rowtag[  1:(i-1) ];
                            rowvals = rowvals[ 1:(i-1) ];
                        }
                        break;            
                    }
                }
                if( length(rowtag) == 0L ){
                    curSliceId = curSliceId - 1L;
                    break;
                }
                rowNameSlices[[curSliceId]] <<- rowtag;
                data = c(rowvals, recursive = TRUE);
                dim(data) = c(length(rowvals[[1]]), length(rowvals));
                data = t(data);
                setSlice(curSliceId, data);
                if( length(rowtag) < fileSliceSize ){
                    break;
                }
                numtxt = formatnum(curSliceId*fileSliceSize);
                message( "Rows read: ", numtxt);
            }
            close(fid)
            if( fileSkipRows == 0 ){
                columnNames <<- paste0("Col_", seq_len(nCols()));
            } else {
                columnNames <<- tail(columnNames, ncol(getSliceRaw(1)));
            }
            if( fileSkipColumns == 0 ){
                cnt = 0L;
                for( sl in seq_len(nSlices()) ){
                    nr = length(getSliceRaw(sl));
                    rowNameSlices[[sl]] <<- paste0("Row_",cnt + (1:nr));
                    cnt = cnt + nr;
                }
            }
            rowNameSlices <<- rowNameSlices[seq_len(curSliceId)];
            message("Rows read: ", nRows(), " done.");
            return(invisible(.self));
        },
        SaveFile = function(filename){
            if( nSlices() == 0 ){
                message("No data to save");
                return();
            }
            fid = file(filename,"wt");
            for( sl in seq_len(nSlices()) ){
                z = getSlice(sl);
                rownames(z) = rowNameSlices[[sl]];
                colnames(z) = columnNames;
                write.table(
                        x = z,
                        file = fid,
                        sep = "\t", 
                        col.names = (if(sl == 1){NA}else{FALSE}));
            }
            close(fid);
        },
        nRows = function(){
            s = 0L;
            for(sl in seq_len(nSlices())){
                s = s + nrow(getSliceRaw(sl));
            }
            return( s );
        },
        nCols = function(){
            if( nSlices() == 0L ){
                return(0L);
            } else {
                return( ncol(getSliceRaw(1L)) );
            }
        },
        Clear = function(){
            for( sl in seq_len(nSlices()) ){
                rm(list = paste0(sl), envir = dataEnv);
            }
            nSlices1 <<- 0L;
            rowNameSlices <<- list();
            columnNames <<- character();
            return(invisible(.self));
        },
        IsCombined = function(){
            return( nSlices() <= 1L );
        },
        GetAllRowNames = function(){
            return( c(rowNameSlices, recursive=TRUE) );
        },
        GetNRowsInSlice = function(sl){
            return( length( rowNameSlices[[sl]] ) );
        },
        SetNanRowMean = function(){
            if( (nCols() == 0L) ){
                return(invisible(.self));
            }
            for( sl in seq_len(nSlices()) ){
                slice = getSlice(sl);
                if( any(is.na(slice)) ){
                    rowmean = rowMeans(slice, na.rm = TRUE);
                    rowmean[is.na(rowmean)] = 0L;
                    for( j in which(!complete.cases(slice)) ){
                        where1 = is.na(slice[j, ]);
                        slice[j, where1] = rowmean[j];
                    }
                    setSlice(sl, slice);
                }
            }
            return(invisible(.self));
        },
        RowStandardizeCentered = function(){
            for(sl in seq_len(nSlices()) ){
                slice = getSlice(sl);
                div = sqrt( rowSums(slice^2) );
                div[ div == 0 ] = 1;
                setSlice(sl, slice/div);
            }
            return(invisible(.self));
        },
        CombineInOneSlice = function(){
            if( nSlices() <= 1L ){
                return(invisible(.self));            
            }
            nc = nCols();
            nr = nRows();
            datatypes = c("raw","integer","double");
            datafuns = c(as.raw, as.integer, as.double);
            datatype = character(nSlices());
            for(sl in seq_len(nSlices())){
                datatype[sl] = typeof(getSliceRaw(sl));
            }
            mch = max(match(datatype, datatypes, nomatch = length(datatypes)));
            datafun = datafuns[[mch]];
            newData = matrix(datafun(0), nrow = nr, ncol = nc);
            offset = 0;
            for(sl in seq_len(nSlices())){
                if( mch == 1 ){
                    slice = getSliceRaw(sl);
                } else {
                    slice = getSlice(sl);
                }
                newData[ offset + (1:nrow(slice)),] = datafun(slice);
                setSlice(sl, numeric());
                offset = offset + nrow(slice);
            }
            
            nSlices1 <<- 1L;
            setSliceRaw(1L, newData);
            rm(newData);
            
            newrowNameSlices = GetAllRowNames();
            rowNameSlices <<- list(newrowNameSlices)
            return(invisible(.self));
        },
        ResliceCombined = function(sliceSize = -1){
            if( sliceSize > 0L ){
                fileSliceSize <<- sliceSize;
            }
            if( fileSliceSize <= 0 ){
                fileSliceSize <<- 1000;
            }
            if( IsCombined() ){
                nRows1 = nRows();
                if(nRows1 == 0L){
                    return(invisible(.self));
                }
                newNSlices = floor((nRows1 + fileSliceSize - 1)/fileSliceSize);
                oldData = getSliceRaw(1L);
                #oldNames = rowNameSlices[[1]];
                newNameslices = vector("list",newNSlices)
                for( sl in seq_len(newNSlices) ){
                    fr = (sl-1)*fileSliceSize + 1;
                    to = min(sl*fileSliceSize, nRows1);
                    range = fr:to;
                    newpart = oldData[range, , drop = FALSE];
                    if( is.raw(oldData) ){
                        setSliceRaw( sl, newpart);
                    } else {
                        setSlice( sl, newpart);
                    }
                    newNameslices[[sl]] = rowNameSlices[[1]][range];
                }
                rowNameSlices <<- newNameslices ;
            } else {
                stop("Reslice of a sliced matrix is not supported.",
                    "Use CombineInOneSlice first.");
            }
            return(invisible(.self));
        },
        Clone = function(){
            clone = SlicedData$new();
            for(sl in seq_len(nSlices()) ){
                clone$setSliceRaw(sl,getSliceRaw(sl));
            }
            clone$rowNameSlices = rowNameSlices;
            clone$columnNames = columnNames;
            clone$fileDelimiter = fileDelimiter;
            clone$fileSkipColumns = fileSkipColumns;
            clone$fileSkipRows = fileSkipRows;
            clone$fileSliceSize = fileSliceSize;
            clone$fileOmitCharacters = fileOmitCharacters;
            return( clone );        
        },
        RowMatrixMultiply = function(multiplier){
            for(sl in seq_len(nSlices()) ){
                setSlice(sl, getSlice(sl) %*% multiplier);
            }
            return(invisible(.self));
        },
        ColumnSubsample = function(subset){
            for(sl in seq_len(nSlices()) ){
                setSliceRaw(sl, getSliceRaw(sl)[, subset, drop = FALSE]);
            }
            columnNames <<- columnNames[subset];
            return(invisible(.self));
        },
        RowReorderSimple = function(ordr){
            # had to use an inefficient and dirty method
            # due to horrible memory management in R
            if( (typeof(ordr) == "logical") && all(ordr) ){
                return(invisible(.self));
            }
            if( (length(ordr) == nRows()) && 
                all(ordr == (1:length(ordr))) ){
                return(invisible(.self));
            }
            CombineInOneSlice();
            gc();
            setSliceRaw( 1L, getSliceRaw(1L)[ordr, ] );
            rowNameSlices[[1]] <<- rowNameSlices[[1]][ordr];
            gc();
            ResliceCombined();
            gc();
            return(invisible(.self));
        },
        RowReorder = function(ordr){
            # transform logical into indices 
            if( typeof(ordr) == "logical" ){
                if( length(ordr) == nRows() ){
                    ordr = which(ordr);
                } else {
                    stop("Parameter \"ordr\" has wrong length")
                }
            }
            ## first, check that anything has to be done at all
            if( (length(ordr) == nRows()) && 
                all(ordr == (1:length(ordr))) ){
                return(invisible(.self));
            }
            ## check bounds
            #if( (min(ordr) < 1) || (max(ordr) > nRows()) ){
            #    stop("Parameter \"ordr\" is out of bounds");
            #}
            ## slice the data into individual rows
            all_rows = vector("list", nSlices())
            for( i in seq_len(nSlices()) ){
                slice = getSliceRaw(i)
                all_rows[[i]] = split(slice, 1:nrow(slice))
                setSliceRaw(i,numeric())
            }
            gc();
            all_rows = unlist(all_rows, recursive=FALSE, use.names = FALSE);
            ## Reorder the rows
            all_rows = all_rows[ordr];
            ## get row names
            all_names = GetAllRowNames();
            ## erase the set
            rowNameSlices <<- list();
            ## sort names
            all_names = all_names[ordr];
            ##
            ## Make slices back
            nrows = length(all_rows);
            nSlices1 <<- as.integer((nrows+fileSliceSize-1)/fileSliceSize);
            
            rowNameSlices1 = vector("list", nSlices1);
            for( i in seq_len(nSlices1) ){
                fr = 1 + fileSliceSize*(i-1);
                to = min( fileSliceSize*i, nrows);
    
                subset = all_rows[fr:to];
                types = unlist(lapply(subset,typeof));
                israw = (types == "raw")
                if(!all(israw == israw[1])){
                    # some raw and some are not
                    subset = lapply(
                        X = subset, 
                        FUN = function(x){
                            if(is.raw(x)){
                                x=as.integer(x);
                                x[x==255] = NA;
                                return(x)
                            } else {
                                return(x)
                            }
                        });
                }
                subset = unlist(subset);
                dim(subset) = c( length(all_rows[[fr]]) , to - fr + 1)
                #subset = matrix(subset, ncol = (to-fr+1));
                if( is.raw(subset) ){
                    setSliceRaw(i, t(subset)); 
                } else {
                    setSlice(i, t(subset)); 
                }
                rowNameSlices1[[i]] = all_names[fr:to];
                all_rows[fr:to] = 0;
                all_names[fr:to] = 0;
            }
            rowNameSlices <<- rowNameSlices1;
            gc();
            return(invisible(.self));
        },
        RowRemoveZeroEps = function(){
            for(sl in seq_len(nSlices()) ){
                slice = getSlice(sl);
                amean = rowMeans(abs(slice));
                remove = (amean < .Machine$double.eps*nCols());
                if(any(remove)){
                    rowNameSlices[[sl]] <<- rowNameSlices[[sl]][!remove];
                    setSlice(sl, slice[!remove, , drop = FALSE]);
                }
            }
            return(invisible(.self));
        },
        FindRow = function(rowname){
            for(sl in seq_len(nSlices()) ){
                mch = match(rowname,rowNameSlices[[sl]], nomatch = 0);
                if( mch > 0 )
                {
                    row = getSlice(sl)[mch[1], , drop=FALSE];
                    rownames(row) = rowname;
                    colnames(row) = columnNames;
                    return( list(slice = sl, item = mch, row = row) );
                }
            }
            return( NULL );
        },
        show = function(){
            cat("SlicedData object. For more information type: ?SlicedData\n");
            cat("Number of columns:", nCols(), "\n");
            cat("Number of rows:", nRows(), "\n");
            cat("Data is stored in", nSlices(), "slices\n");
            if(nCols()>0){
                z = getSlice(1L);
                if(nrow(z)>0){
                    z = z[1:min(nrow(z),10L), 1:min(ncol(z),10L), drop = FALSE];
                    rownames(z) = rowNameSlices[[1]][1:nrow(z)];
                    colnames(z) = columnNames[1:ncol(z)];
                    cat("Top left corner of the first slice (up to 10x10):\n");
                    methods::show(z)
                }
            }        
        }
    ))

setGeneric("nrow");
setMethod(
    f = "nrow",
    signature = "SlicedData",
    function(x){
        return( x$nRows() );
    });
setGeneric("NROW");
setMethod(
    f = "NROW",
    signature = "SlicedData",
    function(x){
        return( x$nRows() );
    });
setGeneric("ncol");
setMethod(
    f = "ncol",
    signature = "SlicedData",
    function(x){
        return( x$nCols() );
    });
setGeneric("NCOL");
setMethod(
    f = "NCOL",
    signature = "SlicedData",
    function(x){
        return( x$nCols() );
    });
setGeneric("dim");
setMethod(
    f = "dim",
    signature = "SlicedData",
    function(x){
        return( c(x$nRows(),x$nCols()) );
    });
setGeneric("colnames");
setMethod(
    f = "colnames",
    signature = "SlicedData",
    function(x){
        return( x$columnNames );
    });
setGeneric("rownames");
setMethod(
    f = "rownames",
    signature = "SlicedData",
    function(x){
        return( x$GetAllRowNames() );
    });
setMethod(
    f = "[[",
    signature = "SlicedData",
    function(x,i){
        return( x$getSlice(i) );
    });
setGeneric("length");
setMethod(
    f = "length",
    signature = "SlicedData",
    function(x){
        return( x$nSlices() );
    });
setMethod(
    f = "[[<-",
    signature = "SlicedData",
    function(x,i,value){
        x$setSlice(i, value);
        return(x);
    });
summary.SlicedData = function(object, ...){
    z = c(
        nCols = object$nCols(),
        nRows = object$nRows(),
        nSlices = object$nSlices());
    return(z);
}

setGeneric("as.matrix");
setMethod(
    f = "as.matrix", 
    signature = "SlicedData",
    function(x){
        if(x$nSlices() == 0){
            return( matrix(0,0,0) );
        }
        if(x$nSlices() > 1){
            copy = x$Clone();
            copy$CombineInOneSlice();
        } else {
            copy = x;
        }
        mat = copy$getSlice(1L);
        rownames(mat) = rownames(copy);
        colnames(mat) = colnames(copy);
        return( mat );
    })
setGeneric("colnames<-")
setMethod(
    f = "colnames<-",
    signature = "SlicedData",
    function(x,value){
        stopifnot( is.character(value) );
        stopifnot( length(value) == x$nCols() );
        x$columnNames = value;
        return(x);
    })
setGeneric("rownames<-")
setMethod(
    f = "rownames<-",
    signature = "SlicedData",
    function(x,value){
        stopifnot( is.character(value) );
        stopifnot( length(value) == x$nRows() );
        start = 1;
        newNameSlices = vector("list", x$nSlices());
        for( i in seq_len(x$nSlices()) ){
            nr = nrow(x$getSliceRaw(i));
            newNameSlices[[i]] = value[ start:(start+nr-1) ];
            start = start + nr;
        }
        x$rowNameSlices = newNameSlices; 
        return(x);
    })
setGeneric("rowSums")
setMethod(
    f = "rowSums",
    signature = "SlicedData",
    function(x, na.rm = FALSE, dims = 1L){
        if(x$nSlices() == 0){
            return( numeric() );
        }
        stopifnot( dims == 1 );
        thesum = vector("list", x$nSlices());
        for( i in seq_len(x$nSlices()) ){
            thesum[[i]] = rowSums(x$getSlice(i), na.rm);
        }
        return(unlist(thesum, recursive = FALSE, use.names = FALSE));
    })
setGeneric("rowMeans")
setMethod(
    f = "rowMeans",
    signature = "SlicedData",
    function(x, na.rm = FALSE, dims = 1L){
        if(x$nSlices() == 0){
            return( numeric() );
        }
        stopifnot( dims == 1 );
        thesum = vector("list", x$nSlices());
        for( i in seq_len(x$nSlices()) ){
            thesum[[i]] = rowMeans(x$getSlice(i), na.rm);
        }
        return(unlist(thesum, recursive = FALSE, use.names = FALSE));
    })
setGeneric("colSums")
setMethod(
    f = "colSums",
    signature = "SlicedData",
    function(x, na.rm = FALSE, dims = 1L){
        if(x$nCols() == 0){
            return( numeric() );
        }
        stopifnot( dims == 1 );
        thesum = 0;
        for( i in seq_len(x$nSlices()) ){
            thesum = thesum + colSums(x$getSlice(i), na.rm);
        }
        return(thesum);
    })
setGeneric("colMeans")
setMethod(
    f = "colMeans",
    signature = "SlicedData",
    function(x, na.rm = FALSE, dims = 1L){
        if(x$nCols() == 0){
            return( numeric() );
        }
        stopifnot( dims == 1 );
        thesum = 0;
        thecounts = x$nRows();
        for( i in seq_len(x$nSlices()) ){
            slice = x$getSlice(i);
            thesum = thesum + colSums(slice, na.rm);
            if( na.rm ){
                thecounts = thecounts - colSums(is.na(slice));
            }
        }
        return(thesum/thecounts);
    })

.listBuilder = setRefClass(".listBuilder",
    fields = list(
        dataEnv = "environment",
        n = "integer"
    ),
    methods = list(
        initialize = function(){
            dataEnv <<- new.env(hash = TRUE);
            n <<- 0L;
            return(.self);
        },
        add = function(x){
            if(length(x) > 0){
                n <<- n + 1L;
                assign(paste0(n), x, dataEnv );
            }
            return(.self);
        },
        set = function(i,x){
            i = as.integer(i);
            if(length(x) > 0){
                if(i>n)
                    n <<- i;
                assign(paste0(i), x, dataEnv );
            }
            return(.self);
        },
        get = function(i){
            return(base::get(paste0(i),dataEnv));
        },
        list = function(){
            result = vector("list",n);
            for( i in seq_len(n)){
                result[[i]] = .self$get(i);
            }
            return(result);
        },
        unlist = function(){
            return(base::unlist(
                        x = .self$list(),
                        recursive=FALSE,
                        use.names = FALSE));
        },
        show = function(){
            cat(".listBuilder object.\n");
            cat("Iternal object in MatrixEQTL package.\n");
            cat("Number of elements:", .self$n, "\n");
        }
    ))

.histogrammer = setRefClass(".histogrammer",
    fields = list(
        pvbins1 = "numeric",
        statbins1 = "numeric",
        hist.count = "numeric"
    ),
    methods = list(
        initialize = function (pvbins, statbins){
            if(length(pvbins)){
                ord = order(statbins);
                pvbins1 <<- pvbins[ord];
                statbins1 <<- statbins[ord];
                statbins1[length(statbins1)] <<- .Machine$double.xmax;
                hist.count <<- double(length(pvbins)-1);
            }
            return(.self);
        },
        update = function(stats.for.hist){
            t = tabulate(
                    findInterval(
                        x = stats.for.hist,
                        vec = statbins1), 
                    nbins = length(statbins1) - 1L);
            hist.count <<- hist.count + t;
        },
        getResults = function(){
            if(!is.unsorted(pvbins1)){
                return(list(
                        hist.bins = pvbins1,
                        hist.counts = hist.count ));
            } else {
                return(list(
                        hist.bins = rev(pvbins1),
                        hist.counts = rev(hist.count)));
            }
        }
    ))

.minpvalue = setRefClass(".minpvalue",
    fields = list(
        sdata = ".listBuilder",
        gdata = ".listBuilder"
    ),
    methods = list(
        initialize = function(snps, gene){
            sdata <<- .listBuilder$new();
            for( ss in seq_len(snps$nSlices()) ){
                sdata$set( ss, double(snps$GetNRowsInSlice(ss)));
            }
            gdata <<- .listBuilder$new();
            for( gg in seq_len(gene$nSlices()) ){
                gdata$set( gg, double(gene$GetNRowsInSlice(gg)));
            }
            return(.self);
        },
        update = function(ss, gg, astat){
            gmax = gdata$get(gg)
            z1 = max.col(astat,ties.method="first");
            z11 = astat[1:nrow(astat) + nrow(astat) * (z1 - 1)];
            gmax = pmax(gmax, z11);
            gdata$set(gg, gmax);
            
            smax = sdata$get(ss)
            z22 = apply(astat,2,max);
            smax = pmax(smax, z22);
            sdata$set(ss, smax);
            return(.self);
        },
        updatecis = function(ss, gg, select.cis, astat){
            if(length(astat)>0)
            {
                byrows = aggregate(
                                x = astat,
                                by = list(row = select.cis[,1]),
                                FUN = max);
                bycols = aggregate(
                                x = astat,
                                by = list(col=select.cis[,2]),
                                FUN = max);
    
                gmax = gdata$get(gg);
                gmax[byrows$row] = pmax(gmax[byrows$row], byrows$x)
                gdata$set(gg, gmax);
                
                smax = sdata$get(ss)
                smax[bycols$col] = pmax(smax[bycols$col], bycols$x)
                sdata$set(ss, smax);
            }            
            return(.self);
        },
        getResults = function(snps, gene, pvfun){
            min.pv.snps = pvfun(sdata$unlist());
            names(min.pv.snps) = rownames(snps);
            min.pv.gene = pvfun(gdata$unlist());
            names(min.pv.gene) = rownames(gene);
            return(list(min.pv.snps = min.pv.snps, min.pv.gene = min.pv.gene));
        }
    ))

.OutputSaver_FRD = setRefClass(".OutputSaver_FRD",
    fields = list(
        sdata = ".listBuilder",
        gdata = ".listBuilder",
        cdata = ".listBuilder",
        bdata = ".listBuilder",
        fid = "list",
        testfun1 = "list",
        pvfun1 = "list"
    ),
    methods = list(
        initialize = function (){
            sdata <<- .listBuilder$new();
            gdata <<- .listBuilder$new();
            cdata <<- .listBuilder$new();
            bdata <<- .listBuilder$new();
            fid <<- list(0);
            testfun1 <<- list(0);
            pvfun1 <<- list(0);
            return(.self);
        },
        start = function(filename, statistic_name,
                         unused1, unused2,
                         testfun, pvfun){
            testfun1 <<- list(testfun);
            pvfun1 <<- list(pvfun);
            if(length(filename) > 0){
                if( is.character(filename) ){
                    fid <<- list(file(
                                    description = filename,
                                    open = "wt",
                                    blocking = FALSE,
                                    raw = FALSE),
                                TRUE);
                } else {
                    fid <<- list(filename, FALSE)
                }
                writeLines(
                    text = paste0(
                                "SNP\tgene\t",
                                statistic_name,
                                "\tp-value\tFDR"),
                    con = fid[[1]]);
            } else {
                fid <<- list();
            }
        },
        update = function(spos, gpos, sta, beta = NULL){
            if(length(sta)>0){
                sdata$add(spos);
                gdata$add(gpos);
                cdata$add(sta );
                if(!is.null(beta ))
                    bdata$add(beta );
            }
            return(.self);
        },
        getResults = function( gene, snps, FDR_total_count){
            pvalues = NULL;
            if(cdata$n > 0){
                tests = testfun1[[1]](cdata$unlist());
                cdata <<- .listBuilder$new();
                
                ord = sort.list(abs(tests), decreasing = TRUE);
                
                tests = tests[ord];
                pvalues = pvfun1[[1]](tests);
                
                FDR = pvalues * FDR_total_count / (1:length(pvalues));
                FDR[length(FDR)] = min(FDR[length(FDR)], 1);
                FDR = rev(cummin(rev(FDR)));
                
                snps_names  = rownames(snps)[sdata$unlist()[ord]];
                sdata <<- .listBuilder$new();
                gene_names  = rownames(gene)[gdata$unlist()[ord]];
                gdata <<- .listBuilder$new();
                 
                beta = NULL;
                if(bdata$n > 0)
                    beta = bdata$unlist()[ord];
                
                if(length(fid)>0)    {    
                    step = 1000; ########### 100000
                    for( part in seq_len(ceiling(length(FDR)/step)) ){
                        fr = (part-1)*step + 1;
                        to = min(part*step, length(FDR));
                        dump = data.frame(
                                    snps_names[fr:to],
                                    gene_names[fr:to],
                                    if(is.null(beta)){
                                        tests[fr:to] 
                                    } else {
                                        list(
                                            beta[fr:to],
                                            tests[fr:to])
                                    },
                                    pvalues[fr:to],
                                    FDR[fr:to], 
                                    row.names = NULL, 
                                    check.rows = FALSE, 
                                    check.names = FALSE, 
                                    stringsAsFactors = FALSE);
                        write.table(
                                x = dump,
                                file = fid[[1]],
                                quote = FALSE,
                                sep = "\t",
                                row.names = FALSE,
                                col.names = FALSE);
                    }
                 }
            } else {
                message("No significant associations were found.\n",
                    file = if(length(fid)>0){fid[[1]]}else{""});
            }
            if(length(fid)>0){    
                if(fid[[2]]){
                    close(fid[[1]]);
                }
                fid <<- list();
            }
             
            if(!is.null(pvalues)){
                eqtls = list( 
                            snps = snps_names,
                            gene = gene_names,
                            statistic = tests,
                            pvalue = pvalues,
                            FDR = FDR);
                 if(!is.null(beta))
                     eqtls$beta = beta;
            } else {
                eqtls = list( 
                            snps = character(),
                            gene = character(),
                            beta = numeric(),
                            statistic = numeric(),
                            pvalue = numeric(),
                            FDR = numeric());
            }
            return(list(eqtls = data.frame(eqtls)));
        }
    )
)


.OutputSaver_direct = setRefClass(".OutputSaver_direct",
    fields = list(
        gene_names = "character",
        snps_names = "character",
        fid = "list",
        testfun1 = "list",
        pvfun1 = "list"
    ),
    methods = list(
        initialize = function(){
            gene_names <<- character(0);
            snps_names <<- character(0);
            fid <<- list(0);
            testfun1 <<- list(0);
            pvfun1 <<- list(0);
            return(.self);
        },
        start = function(filename, statistic_name, snps, gene, testfun, pvfun){
            # I hope the program stops if it fails to open the file
            if( is.character(filename) ){
                fid <<- list(file(
                            description = filename,
                            open = "wt",
                            blocking = FALSE,
                            raw = FALSE), TRUE);
            } else {
                fid <<- list(filename, FALSE);
            }
            writeLines(
                text = paste0("SNP\tgene\t", statistic_name, "\tp-value"),
                con = fid[[1]]);
            gene_names <<- rownames(gene);
            snps_names <<- rownames(snps);
            testfun1 <<- list(testfun);
            pvfun1 <<- list(pvfun);
        },
        update = function(spos, gpos, sta, beta = NULL){
            if( length(sta) == 0 )
                return();
            sta = testfun1[[1]](sta);
            lst = list(
                    snps = snps_names[spos],
                    gene = gene_names[gpos],
                    beta = beta,
                    statistic = sta,
                    pvalue = pvfun1[[1]](sta));
            lst$beta = lst$beta;
            
            dump2 = data.frame(
                        lst,
                        row.names = NULL,
                        check.rows = FALSE,
                        check.names = FALSE,
                        stringsAsFactors = FALSE);
            write.table(
                        x = dump2,
                        file = fid[[1]],
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = FALSE);
        },
        getResults = function(...){
            if( length(fid)>0 ){    
                if(fid[[2]]){
                    close(fid[[1]]);
                }
                fid <<- list();
            }
            gene_names <<- character(0);
            snps_names <<- character(0);
            return(list());
        }
    )
)

.my.pmin = function(x, val){
    # minimum "pmin" function that can handle empty array
    if(length(x) == 0){
        return(x)
    } else {
        return(pmin.int(x,val));
    }    
}

.my.pmax = function(x, val){
    # minimum "pmin" function that can handle empty array
    if(length(x) == 0){
        return(x)
    } else {
        return(pmax.int(x,val));
    }    
}

.pv.nz = function(x){return( .my.pmax(x, .Machine$double.xmin) )}
 
.SetNanRowMean = function(x){
    if( any(is.na(x)) ){
        rowmean = rowMeans(x, na.rm = TRUE);
        rowmean[ is.na(rowmean) ] = 0;
        for( j in which(!complete.cases(x)) ){
            where1 = is.na( x[j, ] );
            x[j,where1] = rowmean[j];
        }
    }
    return(x);
}

.SNP_process_split_for_ANOVA = function(x,n.groups){
    # split into 2 dummy variables (or more)
    
#     # Get the number of ANOVA groups
#     n.groups = options("MatrixEQTL.ANOVA.categories")[[1]];
#     if( is.null(n.groups))
#         n.groups = 3;
    
    # Unique values in x (make sure it has length of n.groups);
    uniq = unique(as.vector(x));
    uniq = uniq[!is.na(uniq)];
    if( length(uniq) > n.groups ){
        stop("More than declared number of genotype categories ",
                "is detected in ANOVA");
    } else if ( length(uniq) < n.groups ){
        uniq = c(uniq, rep(min(uniq)-1, n.groups-length(uniq)));
    }
    
    # Table of frequencies for each variable (row)
    freq = matrix(0, nrow(x), n.groups);
    for(i in seq_len(n.groups) ){
        freq[ ,i] = rowSums(x==uniq[i], na.rm = TRUE);
    }
    # remove NA-s from x for convenience
    x[is.na(x)] = min(uniq) - 2;
    
    # Output list of matrices
    rez = vector("list", n.groups - 1);

    # Skip the most frequent value
    md = apply(freq, 1, which.max); # most frequent value for each variable
    freq[ cbind(1:nrow(x),md) ] = -1;
    
    # The rest form dumm
    for(j in seq_len(n.groups-1) ){
        md = apply(freq, 1, which.max); 
        freq[ cbind(1:nrow(x),md) ] = -1;
        rez[[j]] = (x == uniq[md]);
    }            
    return( rez );
}

Matrix_eQTL_engine = function(
                        snps, 
                        gene, 
                        cvrt = SlicedData$new(), 
                        output_file_name, 
                        pvOutputThreshold = 1e-5, 
                        useModel = modelLINEAR, 
                        errorCovariance = numeric(), 
                        verbose = TRUE,
                        pvalue.hist = FALSE,
                        min.pv.by.genesnp = FALSE,
                        noFDRsaveMemory = FALSE){
    rez = Matrix_eQTL_main(
                snps = snps, 
                gene = gene, 
                cvrt = cvrt, 
                output_file_name = output_file_name, 
                pvOutputThreshold = pvOutputThreshold,
                useModel = useModel, 
                errorCovariance = errorCovariance, 
                verbose = verbose,
                pvalue.hist = pvalue.hist,
                min.pv.by.genesnp = min.pv.by.genesnp,
                noFDRsaveMemory = noFDRsaveMemory);
    return( rez );
}

Matrix_eQTL_main = function(    
                        snps, 
                        gene, 
                        cvrt = SlicedData$new(), 
                        output_file_name = "", 
                        pvOutputThreshold = 1e-5,
                        useModel = modelLINEAR, 
                        errorCovariance = numeric(), 
                        verbose = TRUE, 
                        output_file_name.cis = "", 
                        pvOutputThreshold.cis = 0,
                        snpspos = NULL, 
                        genepos = NULL,
                        cisDist = 1e6,
                        pvalue.hist = FALSE,
                        min.pv.by.genesnp = FALSE,
                        noFDRsaveMemory = FALSE){
    ########################## Basic variable checks ##########################
    {
        if(!is.null(output_file_name.cis))
            if(nchar(output_file_name.cis)==0)
                output_file_name.cis = NULL;
        
        if(!is.null(output_file_name))
            if(nchar(output_file_name)==0)
                output_file_name = NULL;
        
        # status("Performing basic checks of the input variables");
        stopifnot( is(gene, "SlicedData") );
        stopifnot( is(snps, "SlicedData") );
        stopifnot( is(cvrt, "SlicedData") );
        
        # Check dimensions
        if( min(snps$nRows(),snps$nCols()) == 0 )
            stop("Empty genotype dataset");
        if( min(gene$nRows(),gene$nCols()) == 0 )
            stop("Empty expression dataset");
        if( snps$nCols() != gene$nCols() )
            stop("Different number of samples in ",
                "the genotype and gene expression files");
        if( cvrt$nRows()>0 ){
            if( snps$nCols() != cvrt$nCols() )
                stop("Wrong number of samples in the matrix of covariates");
        }

        stopifnot( is(pvOutputThreshold, "numeric") );
        stopifnot( length(pvOutputThreshold) == 1 );
        stopifnot( pvOutputThreshold >= 0 );
        stopifnot( pvOutputThreshold <= 1 );

        stopifnot( is(noFDRsaveMemory, "logical") );
        stopifnot( length(noFDRsaveMemory) == 1 );

        if( pvOutputThreshold > 0 ){
            stopifnot( !((length(output_file_name) == 0) && noFDRsaveMemory) );
            stopifnot( length(output_file_name) <= 1 );
            if( length(output_file_name) == 1 ){
                stopifnot( any( class(output_file_name) %in% 
                            c("character","connection") ) );
            }
        }
        
        stopifnot( is(pvOutputThreshold.cis, "numeric") );
        stopifnot( length(pvOutputThreshold.cis) == 1 );
        stopifnot( pvOutputThreshold.cis >= 0 );
        stopifnot( pvOutputThreshold.cis <= 1 );
        stopifnot(  (pvOutputThreshold == 0) || 
                    (pvOutputThreshold.cis == 0) || 
                    (pvOutputThreshold <= pvOutputThreshold.cis));
        stopifnot( (pvOutputThreshold > 0) | (pvOutputThreshold.cis > 0) );

        stopifnot( is(useModel, class(modelLINEAR)) );
        stopifnot( length(useModel) == 1 );
        stopifnot( useModel %in% c(modelLINEAR, modelANOVA, modelLINEAR_CROSS));
        if( useModel %in%  c(modelLINEAR, modelLINEAR_CROSS) ){
            if( snps$nCols() <= cvrt$nRows() + 1 + 1){
                stop("The number of covariates exceeds the number of samples.",
                    "\nLinear regression can not be fit.")
            }
        }
        if( useModel == modelLINEAR_CROSS ){
            if( cvrt$nRows() == 0 ){
                stop("Model \"modelLINEAR_CROSS\" requires",
                    " at least one covariate");
            }
        }
        if( useModel == modelANOVA ){
            n.anova.groups = getOption("MatrixEQTL.ANOVA.categories", 3);
            stopifnot( n.anova.groups == floor(n.anova.groups) );
            stopifnot( n.anova.groups >= 3 );
            if( snps$nCols() <= cvrt$nRows() + n.anova.groups){
                stop("The number of covariates exceeds the number of samples.",
                    "\nLinear regression (ANOVA) can not be fit.")
            }
        }
        
        stopifnot( is(verbose, "logical") );
        stopifnot( length(verbose) == 1 );

        stopifnot( is(min.pv.by.genesnp, "logical") );
        stopifnot( length(min.pv.by.genesnp) == 1 );
    
        if( pvOutputThreshold.cis > 0 ){
            stopifnot( (length(output_file_name.cis) > 0) || !noFDRsaveMemory );
            stopifnot( length(output_file_name.cis) <= 1 );
            if( length(output_file_name.cis) == 1 ){
                stopifnot( any( class(output_file_name.cis) %in%
                            c("character","connection") ) );
            }

            stopifnot( is(snpspos, "data.frame") );
            stopifnot( ncol(snpspos) == 3 );
            stopifnot( nrow(snpspos) > 0 );
            stopifnot( is.numeric(snpspos[[3]]) );
            stopifnot( !any(is.na(snpspos[,3])) );
            stopifnot( is(genepos, "data.frame") );
            stopifnot( ncol(genepos) == 4 );
            stopifnot( nrow(genepos) > 0 );
            stopifnot( is.numeric(genepos[[3]]) );
            stopifnot( is.numeric(genepos[[4]]) );
            stopifnot( !any(is.na(genepos[[3]])) );
            stopifnot( !any(is.na(genepos[[4]])) );
            stopifnot( nzchar(output_file_name.cis) );
        }
        
        if( pvOutputThreshold > 0 )
            stopifnot( nzchar(output_file_name) );
        
        stopifnot( is.numeric(errorCovariance) );
        errorCovariance = as.matrix(errorCovariance);
        if( length(errorCovariance) > 0 ){
            if( nrow(errorCovariance) != ncol(errorCovariance) ){
                stop("The covariance matrix is not square");
            }
            if( nrow(errorCovariance) != snps$nCols() ){
                stop("The covariance matrix size",
                    " does not match the number of samples");
            }
            if( !all(errorCovariance == t(errorCovariance)) ){
                stop("The covariance matrix is not symmetric");
            }
        }
    }
    ############################## Initial setup ##############################
    {
        gene.std = .listBuilder$new();
        snps.std = .listBuilder$new();
        
        dont.clone.gene = getOption(
                            x = "MatrixEQTL.dont.preserve.gene.object",
                            default = FALSE);
        if(is.null(dont.clone.gene))
            dont.clone.gene = FALSE;
        
        if( !dont.clone.gene )
            gene = gene$Clone();
        # snps = snps$Clone(); # snps is read only
        cvrt = cvrt$Clone();

        params = list(
            output_file_name = output_file_name, 
            pvOutputThreshold = pvOutputThreshold,
            useModel = useModel, 
            errorCovariance = errorCovariance , 
            verbose = verbose, 
            output_file_name.cis = output_file_name.cis, 
            pvOutputThreshold.cis = pvOutputThreshold.cis,
            cisDist = cisDist ,
             pvalue.hist = pvalue.hist,
            min.pv.by.genesnp = min.pv.by.genesnp);

        if( verbose ){
            lastTime = 0;
            status = function(text){
                # gc();
                newTime = proc.time()[3];
                if(lastTime != 0){
                    message("Task finished in ", round(newTime-lastTime,3),
                            " seconds");
                }
                message(text);
                lastTime <<- newTime;
            }
        } else {
            status = function(text){}
        }
        start.time = proc.time()[3];
    }
    #################### Error covariance matrix processing ####################
    {
        if( length(errorCovariance) > 0 ){
            status("Processing the error covariance matrix");
            eig = eigen(errorCovariance, symmetric = TRUE)
            d = eig$values;
            v = eig$vectors;
            #  errorCovariance == v %*% diag(d) %*% t(v)
            #  errorCovariance^0.5 == v*sqrt(d)*v" (no single quotes anymore)
            #  errorCovariance^(-0.5) == v*diag(1./sqrt(diag(d)))*v"
            if( any(d<=0) ){
                stop("The covariance matrix is not positive definite");
            }
            correctionMatrix = v %*% diag(1./sqrt(d)) %*% t(v);
            rm( eig, v, d, errorCovariance )
        } else {
            rm( errorCovariance );
            correctionMatrix = numeric();
        }
    }
    ##################### Matching gene and SNPs locations #####################
    if( pvOutputThreshold.cis > 0 ){
        status("Matching data files and location files");
        
        # names in the input data    
        gene_names = rownames(gene);
        snps_names = rownames(snps);
        
        # gene range, set: left<right
        if( any(genepos[[3]] > genepos[[4]]) ){
            temp3 = genepos[[3]];
            temp4 = genepos[[4]];
            genepos[[3]] = pmin(temp3, temp4);
            genepos[[4]] = pmax(temp3, temp4);
            rm(temp3, temp4);
        }
        
        # match with the location data
        genematch = match(gene_names, genepos[[1]], nomatch = 0L);
        usedgene = matrix(FALSE, nrow(genepos), 1); 
        # genes in "genepos" that are matching  "gene_names"
        usedgene[ genematch ] = TRUE;
        if( max(genematch) == 0 ){
            stop("Gene names do not match those in the gene location file.");
        }
        message(sum(genematch>0), " of ", length(gene_names), 
                " genes matched");
        
        
        snpsmatch = match( snps_names, snpspos[ ,1],  nomatch = 0L);
        usedsnps = matrix(FALSE, nrow(snpspos), 1);
        usedsnps[ snpsmatch ] = TRUE;
        if( max(snpsmatch) == 0 ){
            stop("SNP names do not match those in the SNP location file.");
        }
        message(sum(snpsmatch>0), " of ", length(snps_names), 
                " SNPs matched\n");
        
        # list used chr names
        chrNames = unique(c( 
                    as.character(unique(snpspos[[2]][usedsnps])), 
                    as.character(unique(genepos[[2]][usedgene]))));
        chrNames = chrNames[ sort.list( 
                                x = suppressWarnings(as.integer(chrNames)),
                                method = "radix",
                                na.last = TRUE)];
        
        # match chr names
        genechr = match(genepos[[2]], chrNames);
        snpschr = match(snpspos[[2]], chrNames);
        
        # max length of a chromosome
        chrMax = max(snpspos[[3]][usedsnps], genepos[[4]][usedgene], na.rm = TRUE) + 
                    cisDist;
        
        # Single number location for all rows in "genepos" and "snpspos"
        genepos2 = as.matrix(genepos[3:4]) + (genechr-1)*chrMax;
        snpspos2 = as.matrix(snpspos[[3]]) + (snpschr-1)*chrMax;
        
        # the final location arrays;
        snps_pos = matrix(0, length(snps_names), 1);
        snps_pos[snpsmatch>0, ] = snpspos2[snpsmatch, , drop = FALSE];
        snps_pos[rowSums(is.na(snps_pos))>0, ] = 0;
        snps_pos[snps_pos==0] = (length(chrNames)+1) * (chrMax+cisDist);
        rm(snps_names, snpsmatch, usedsnps, snpschr, snpspos2)
        
        gene_pos = matrix(0, length(gene_names), 2);
        gene_pos[genematch>0, ] = genepos2[genematch, , drop = FALSE];
        gene_pos[rowSums(is.na(gene_pos))>0, ] = 0;
        gene_pos[gene_pos==0] = (length(chrNames)+2) * (chrMax+cisDist);
        rm(gene_names, genematch, usedgene, genechr, genepos2)
        rm(chrNames, chrMax);

        if( is.unsorted(snps_pos) ){
            status("Reordering SNPs");
            ordr = sort.list(snps_pos);
            snps$RowReorder(ordr);
            snps_pos = snps_pos[ordr, , drop = FALSE];
            rm(ordr);
        }
        if( is.unsorted(rowSums(gene_pos)) ){
            status("Reordering genes");
            ordr = sort.list(rowSums(gene_pos));
            gene$RowReorder(ordr);
            gene_pos = gene_pos[ordr, , drop = FALSE];
            rm(ordr);
        }
        
        # Slice it back
        geneloc = vector("list", gene$nSlices())
        gene_offset = 0;
        for(gc in seq_len(gene$nSlices()) ){
            nr = gene$GetNRowsInSlice(gc);
            geneloc[[gc]] = gene_pos[gene_offset + (1:nr), , drop = FALSE];
            gene_offset = gene_offset + nr;    
        }
        rm(gc, gene_offset, gene_pos);
        
        snpsloc = vector("list", snps$nSlices())
        snps_offset = 0;
        for(sc in seq_len(snps$nSlices()) ){
            nr = snps$GetNRowsInSlice(sc);
            snpsloc[[sc]] = snps_pos[snps_offset + (1:nr), , drop = FALSE];
            snps_offset = snps_offset + nr;    
        }
        rm(nr, sc, snps_offset, snps_pos);
    }
    ########################## Covariates processing ##########################
    {
        status("Processing covariates");
        if( useModel == modelLINEAR_CROSS ){
            last.covariate = as.vector(tail( cvrt$getSlice(cvrt$nSlices()), 1));
        }        
        if( cvrt$nRows()>0 ){
            cvrt$SetNanRowMean();
            cvrt$CombineInOneSlice();
            cvrt = rbind(matrix(1,1,snps$nCols()), cvrt$getSlice(1));
        } else {
            cvrt = matrix(1,1,snps$nCols());
        }
        # Correct for the error covariance structure
        if( length(correctionMatrix)>0 ){
            cvrt = cvrt %*% correctionMatrix;
        }
        # Orthonormalize covariates
        # status("Orthonormalizing covariates");
        q = qr(t(cvrt));
        if( min(abs(diag(qr.R(q)))) < .Machine$double.eps * snps$nCols() ){
            stop("Colinear or zero covariates detected");
        }
        cvrt = t( qr.Q(q) );
        rm(q);
    }
    ######################## Gene expression processing ########################
    {
        status("Processing gene expression data (imputation, residualization)");
        # Impute gene expression
        gene$SetNanRowMean();
        # Correct for the error covariance structure
        if( length(correctionMatrix)>0 ){
            gene$RowMatrixMultiply(correctionMatrix);
        }
        # Orthogonolize expression w.r.t. covariates
        # status("Orthogonolizing expression w.r.t. covariates");
        gene_offsets = double(gene$nSlices()+1);
        for( sl in seq_len(gene$nSlices()) ){
            slice = gene$getSlice(sl);
            gene_offsets[sl+1] = gene_offsets[sl] + nrow(slice);
            rowsq1 = rowSums(slice^2);
            slice = slice - tcrossprod(slice,cvrt) %*% cvrt;
            rowsq2 = rowSums(slice^2);
            # kill rows colinear with the covariates
            delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
            slice[delete.rows,] = 0;
            rowsq2[delete.rows] = 1;
            div = sqrt(rowsq2); #sqrt( rowSums(slice^2) );
            gene.std$set(sl, div);
            gene$setSlice(sl, slice / div);
        }
        rm(rowsq1, rowsq2, delete.rows, div);
        rm( sl, slice );
        #gene$RowRemoveZeroEps();
    }
    ############## snps_process, testfun, pvfun, threshfun, afun ##############
    {
        # snps_process - preprocess SNPs slice
        #
        # afun --- abs for signed stats, identity for non-negative
        # threshfun --- internal stat threshold for given p-value
        # testfun --- t or F statistic from the internal one
        # pvfun --- p-value from the t or F statistic
        
        nSamples = snps$nCols();
        nGenes = gene$nRows();
        nSnps  = snps$nRows();
        nCov = nrow(cvrt);
        # nVarTested = length(snps_list); # set in case(useModel)
        # dfNull = nSamples - nCov;
        # d.f. of the full model
        betafun = NULL;
        
        if( useModel == modelLINEAR ){
            snps_process = function(x){
                return( list(.SetNanRowMean(x)) );
            };
            nVarTested = 1;
            dfFull = nSamples - nCov - nVarTested;
            statistic.fun = function(mat_list){
                return( mat_list[[1]] );
            }
            afun = function(x){return(abs(x))};
            threshfun = function(pv){
                thr = qt(pv/2, dfFull, lower.tail = FALSE);
                thr = thr^2;
                thr[is.infinite(thr)] = .Machine$double.xmax;
                thr = sqrt(  thr / (dfFull + thr) );
                thr[pv >= 1] = 0;
                thr[pv <= 0] = 1;
                return( thr );
            }
            testfun = function(x){ 
                return( x * sqrt( dfFull / (1 - .my.pmin(x^2,1))));
            }
            pvfun = function(x){
                return( .pv.nz(pt(-abs(x),dfFull)*2));
            }
            thresh.cis = threshfun(pvOutputThreshold.cis);
            thresh = threshfun(pvOutputThreshold);
            betafun = function(stat, ss, gg, select){
                return(  stat * gene.std$get(gg)[select[,1]] /
                                snps.std$get(ss)[select[,2]]);
            }
        } else 
        if( useModel == modelANOVA ){
            snps_process = function(x){
                .SNP_process_split_for_ANOVA(x,n.anova.groups);
            }
            nVarTested = n.anova.groups - 1;
            dfFull = nSamples - nCov - nVarTested;
            statistic.fun = function(mat_list){
                x = mat_list[[1]]^2;
                for( j in 2:length(mat_list) )
                    x = x + mat_list[[j]]^2;
                return(x);
            }
            afun = identity;
            threshfun = function(pv){
                thr = qf(pv, nVarTested, dfFull, lower.tail = FALSE);
                thr = thr / (dfFull/nVarTested + thr);
                thr[pv >= 1] = 0;
                thr[pv <= 0] = 1;
                return(thr);
            }
            testfun = function(x){
                return( x / (1 - .my.pmin(x,1)) * (dfFull/nVarTested) );
            }
            pvfun = function(x){
                return( .pv.nz(pf(x, nVarTested, dfFull, lower.tail = FALSE)) );
            }
            thresh.cis = threshfun(pvOutputThreshold.cis);
            thresh = threshfun(pvOutputThreshold);
        } else 
        if( useModel == modelLINEAR_CROSS ){
            last.covariate = as.vector( last.covariate );
            snps_process = .SNP_process_split_for_LINEAR_CROSS = function(x){
                out = vector("list", 2);
                out[[1]] = .SetNanRowMean(x);
                out[[2]] = t( t(out[[1]]) * last.covariate );
                return( out );
            };
            nVarTested = 1;
            dfFull = nSamples - nCov - nVarTested - 1;
            statistic.fun = function(mat_list){
                return( mat_list[[2]] / sqrt(1 - mat_list[[1]]^2) );
            }
            afun = function(x){return(abs(x))};
            threshfun = function(pv){
                thr = qt(pv/2, dfFull, lower.tail = FALSE);
                thr = thr^2;
                thr = sqrt(  thr / (dfFull + thr) );
                thr[pv >= 1] = 0;
                thr[pv <= 0] = 1;
                return( thr );
            }
            testfun = function(x){
                return( x * sqrt( dfFull / (1 - .my.pmin(x^2,1))));
            }
            pvfun = function(x){
                return( .pv.nz(pt(-abs(x),dfFull)*2 ));
            }        
            thresh.cis = threshfun(pvOutputThreshold.cis);
            thresh = threshfun(pvOutputThreshold);                
            betafun = function(stat, ss, gg, select){
                return(  stat * gene.std$get(gg)[select[,1]] / 
                                snps.std$get(ss)[select[,2]]);
            }
        }
        params$dfFull = dfFull;
    }
    ######################### Saver class(es) creation #########################
    {
        status("Creating output file(s)");
        if(noFDRsaveMemory){
            if( pvOutputThreshold > 0 ){
                saver.tra = .OutputSaver_direct$new();
            }
            if( pvOutputThreshold.cis > 0 ){
                saver.cis = .OutputSaver_direct$new();
            }
        } else {
            if( pvOutputThreshold > 0 ){
                saver.tra = .OutputSaver_FRD$new();
            }
            if( pvOutputThreshold.cis > 0 ){
                saver.cis = .OutputSaver_FRD$new();
            }
        }
        if( pvOutputThreshold > 0 )
            if( pvOutputThreshold * gene$nRows() * snps$nRows() > 1000000 )
                if(!noFDRsaveMemory)
                    warning("Warning: pvOutputThreshold may be too large.\n",
                            "Expected number of findings > ", 
                            pvOutputThreshold * gene$nRows() * snps$nRows());
        if( (useModel == modelLINEAR) || (useModel == modelLINEAR_CROSS) ){
            statistic_name = "t-stat";
        } else if( useModel == modelANOVA ){
            statistic_name = "F-test";
        }
        if(!is.null(betafun))
            statistic_name = paste0("beta\t", statistic_name);
        if( pvOutputThreshold > 0 )
            saver.tra$start(output_file_name,
                            statistic_name, snps, gene, testfun, pvfun);
        if( pvOutputThreshold.cis > 0 )
            saver.cis$start(output_file_name.cis,
                            statistic_name, snps, gene, testfun, pvfun);
        rm( statistic_name );
    }
    ########################## Some useful functions ##########################
    {
        orthonormalize.snps = function(cursnps, ss){
            for(p in seq_along(cursnps)){
                if(length(correctionMatrix)>0){
                    cursnps[[p]] = cursnps[[p]] %*% correctionMatrix;
                }
                rowsq1 = rowSums(cursnps[[p]]^2);
                cursnps[[p]] = cursnps[[p]] - 
                            tcrossprod(cursnps[[p]],cvrt) %*% cvrt;
                for(w in seq_len(p-1L))
                    cursnps[[p]] = cursnps[[p]] - 
                            rowSums(cursnps[[p]]*cursnps[[w]]) * cursnps[[w]];
                rowsq2 = rowSums(cursnps[[p]]^2);
                delete.rows = (rowsq2 <= rowsq1 * .Machine$double.eps );
                cursnps[[p]][delete.rows,] = 0;
                div = sqrt( rowsq2 );
                div[ delete.rows ] = 1;
                cursnps[[p]] = cursnps[[p]]/div;
            }
            snps.std$set(ss, div);
            return(cursnps);
        }

        if( pvOutputThreshold.cis > 0 ){
            # sn.l = sapply(snpsloc, function(x)x[1] );
            # sn.r = sapply(snpsloc, function(x)tail(x,1) );
            # ge.l = sapply(geneloc, function(x)x[1,1] );
            # ge.r = sapply(geneloc, function(x)x[nrow(x) , 2] );
            sn.l = sapply(snpsloc, "[", 1 );
            sn.r = sapply(snpsloc, tail, 1 );
            ge.l = sapply(geneloc, min);
            ge.r = sapply(geneloc, max);
            
            # For find interval
            ge.l = rev(cummin(rev(ge.l)));
            ge.r = cummax(ge.r);
            
            gg.1 = findInterval( sn.l, ge.r + cisDist + 1) + 1;
            gg.2 = findInterval( sn.r, ge.l - cisDist );
            rm(sn.l, sn.r, ge.l, ge.r);
         }
    }    
    ################### Prepare counters and histogram bins ###################
    {
        pvbins = NULL; # bin edges for p-values
        statbins = 0;  # bin edges for the test statistic (|t| or F)
        do.hist = FALSE;
        if( length(pvalue.hist) == 1 ){
            if(pvalue.hist == "qqplot"){
                pvbins = c(0, 
                        10^rev(seq(0, log10(.Machine$double.xmin)-1, -0.05)));
            } else
            if( is.numeric(pvalue.hist) ){
                pvbins = seq(from = 0, to = 1, length.out = pvalue.hist+1);
            } else
            if( pvalue.hist == TRUE ){
                pvbins = seq(from = 0, to = 1, length.out = 100+1);
            }
        } else
        if( is.numeric(pvalue.hist) && (length(pvalue.hist) > 1) ){
            pvbins = pvalue.hist;
        }
        if( is.null(pvbins) && (pvalue.hist != FALSE) ){
            stop("Wrong value of pvalue.hist.",
                " Must be FALSE, TRUE, \"qqplot\", or numerical");
        }
        do.hist = !is.null(pvbins);
        if( do.hist ){
            pvbins = sort(pvbins);
            statbins = threshfun(pvbins);
            if( pvOutputThreshold > 0){
                hist.all = .histogrammer$new(pvbins, statbins);
            }
            if( pvOutputThreshold.cis > 0){
                hist.cis = .histogrammer$new(pvbins, statbins);
            }
        }
        rm( pvbins, statbins);
        if(min.pv.by.genesnp){
            if( pvOutputThreshold > 0){
                minpv.tra = .minpvalue$new(snps,gene);
            }
            if( pvOutputThreshold.cis > 0){
                minpv.cis = .minpvalue$new(snps,gene);
            }
        }
    }
    ################################ Main loop ################################
    {
        beta = NULL;
        n.tests.all = 0;
        n.tests.cis = 0;
        n.eqtls.tra = 0;
        n.eqtls.cis = 0;
        
        status("Performing eQTL analysis");
        # ss = 1; gg = 1;
        # ss = snps$nSlices(); gg = gene$nSlices();
        
        snps_offset = 0;
        for(ss in seq_len(snps$nSlices())){
        # for(ss in 1:min(2,snps$nSlices())){ #for debug
            cursnps = NULL;
            nrcs = snps$GetNRowsInSlice(ss);
            
            # loop only through the useful stuff
            if(pvOutputThreshold>0){
                loopset = seq_len(gene$nSlices());
            } else {
                loopset = .seq(gg.1[ss],gg.2[ss]);
            }
            
            for( gg in loopset ){
                gene_offset = gene_offsets[gg];
                curgene = gene$getSlice(gg);
                nrcg = nrow(curgene);
                if(nrcg == 0) next;
                
                rp = "";
                
                statistic = NULL;
                select.cis.raw = NULL;
                ## do cis analysis
                if( (pvOutputThreshold.cis > 0) && 
                    (gg >= gg.1[ss]) && 
                    (gg <= gg.2[ss]) ){
                    
                    if( is.null( statistic ) ){
                        if( is.null( cursnps ) ){
                            cursnps = orthonormalize.snps( 
                                    cursnps = snps_process(snps$getSlice(ss)), 
                                    ss = ss);
                        }                    
                        mat = vector("list", length(cursnps));
                        for(d in seq_along(cursnps)){
                            mat[[d]] = tcrossprod(curgene, cursnps[[d]]);
                        }
                        statistic = statistic.fun( mat );
                        astatistic = afun(statistic);
                    }
                    
                    sn.l = findInterval(
                                x = geneloc[[gg]][ ,1] - cisDist - 1,
                                vec = snpsloc[[ss]]);
                    sn.r = findInterval(
                                x = geneloc[[gg]][ ,2] + cisDist,
                                vec = snpsloc[[ss]]);
                    xx = unlist(lapply(
                                X = which(sn.r>sn.l),
                                FUN = function(x){
                                    (sn.l[x]:(sn.r[x]-1)) * 
                                    as.numeric(nrow(statistic)) + x
                                }));
                    select.cis.raw = xx[ astatistic[xx] >= thresh.cis ];
                    select.cis = arrayInd(select.cis.raw, dim(statistic));
                    
                    n.tests.cis = n.tests.cis + length(xx);
                    n.eqtls.cis = n.eqtls.cis + length(select.cis.raw);
                    
                    if( do.hist )    
                        hist.cis$update(astatistic[xx]);
                    
                    if( min.pv.by.genesnp ){
                        temp = double(length(astatistic));
                        dim(temp) = dim(astatistic);
                        temp[xx] = astatistic[xx];
                        minpv.cis$update(ss, gg, temp);
                    }
                    
                    if(!is.null(betafun))
                        beta = betafun(
                                stat = mat[[length(mat)]][select.cis.raw],
                                ss = ss,
                                gg = gg,
                                select = select.cis);
                                        
                    saver.cis$update( 
                                snps_offset + select.cis[, 2],
                                gene_offset + select.cis[, 1],
                                statistic[select.cis.raw],
                                beta);
                    
                    rp = paste0(rp, ", ", formatnum(n.eqtls.cis), " cis-eQTLs");
                }
                ## do trans/all analysis
                if(pvOutputThreshold>0){
                    if( is.null( statistic ) ){
                        if( is.null( cursnps ) ){
                            cursnps = orthonormalize.snps(
                                cursnps = snps_process( snps$getSlice(ss) ),
                                ss = ss);
                        }
                        mat = vector("list", length(cursnps));
                        for(d in seq_along(cursnps)){
                            mat[[d]] = tcrossprod(curgene, cursnps[[d]]);
                        }
                        statistic = statistic.fun( mat );
                        astatistic = afun(statistic);
                    }
    
                    if( do.hist )    
                        hist.all$update(astatistic);
    
                    if(!is.null(select.cis.raw)) 
                        astatistic[xx] = -1;
                    
                    select.tra.raw = which( astatistic >= thresh);
                    select.tra = arrayInd(select.tra.raw, dim(statistic))
                    
                    n.eqtls.tra = n.eqtls.tra + length(select.tra.raw);
                    n.tests.all = n.tests.all + length(statistic);

                    if(!is.null(betafun))
                        beta = betafun(
                                stat = mat[[length(mat)]][select.tra.raw],
                                ss = ss,
                                gg = gg,
                                select = select.tra);
                                        
                    saver.tra$update( 
                                snps_offset + select.tra[, 2],
                                gene_offset + select.tra[, 1],
                                statistic[select.tra.raw],
                                beta);
                    
                    if( min.pv.by.genesnp ) 
                        minpv.tra$update(ss, gg, astatistic)
                    
                    rp = paste0(rp, ", ", formatnum(n.eqtls.tra), 
                            if(pvOutputThreshold.cis > 0){" trans-"}else{" "},
                            "eQTLs")
                }
    
                #gene_offset = gene_offset + nrcg;
                if( !is.null(statistic) ){
                    per = 100*(gg/gene$nSlices() + ss-1) / snps$nSlices();
                    per = formatC(
                            x = floor(per*100)/100,
                            format = "f",
                            width = 5,
                            digits = 2);
                    
                    message(per, "% done", rp);
                }
            } # gg in 1:gene$nSlices()
            snps_offset = snps_offset + nrcs;
        } # ss in 1:snps$nSlices()
    }
    ############################ Results collection ############################
    {
        rez = list(time.in.sec = proc.time()[3] - start.time);
        rez$param = params;
        
        if(pvOutputThreshold.cis > 0){
            rez.cis = list(ntests = n.tests.cis, neqtls = n.eqtls.cis);
            rez.cis = c(rez.cis, saver.cis$getResults(gene, snps, n.tests.cis));
            if(do.hist)
                rez.cis = c(rez.cis, hist.cis$getResults() );
            if(min.pv.by.genesnp)
                rez.cis = c(rez.cis, minpv.cis$getResults(snps, gene, 
                                pvfun = function(x){pvfun(testfun(x))}) );
        }
        
        if(pvOutputThreshold>0){
            rez.all = list(
                        ntests = n.tests.all,
                        neqtls = n.eqtls.tra + n.eqtls.cis);
            if(pvOutputThreshold.cis > 0){
                rez.tra = list(
                        ntests = n.tests.all - n.tests.cis,
                        neqtls = n.eqtls.tra);
                rez.tra = c(rez.tra, 
                    saver.tra$getResults(gene, snps, n.tests.all - n.tests.cis))
            } else {
                rez.all = c(rez.all, 
                    saver.tra$getResults(gene, snps, n.tests.all));
            }
            if(do.hist){
                rez.all = c(rez.all, hist.all$getResults() );
                if(pvOutputThreshold.cis > 0){
                    rez.tra$hist.bins = rez.all$hist.bins;
                    rez.tra$hist.counts =   rez.all$hist.counts - 
                                            rez.cis$hist.counts;
                }
            }
            if(min.pv.by.genesnp){
                if(pvOutputThreshold.cis > 0){
                    rez.tra = c(rez.tra, minpv.tra$getResults(snps, gene, 
                                pvfun = function(x){pvfun(testfun(x))}) );
                } else {
                    rez.all = c(rez.all, minpv.tra$getResults(snps, gene, 
                                pvfun = function(x){pvfun(testfun(x))}) );
                }
            }
        }
    
        if(exists("rez.all")>0)
            rez$all = rez.all;
        if(exists("rez.tra")>0)
            rez$trans = rez.tra;
        if(exists("rez.cis")>0)
            rez$cis = rez.cis;    
        
        class(rez) = c(class(rez), "MatrixEQTL");
        status("");
    }
    ################################### Done ###################################
    return(rez);
}

.histme = function(m, name1, name2, ...){
    cnts = m$hist.counts;
    bins = m$hist.bins;
    ntst = m$ntests;
    centers = 0.5 * (bins[-1L] + bins[-length(bins)]);
    density = 0.5 / (bins[-1L] - centers) * cnts / ntst;
    ntext = paste0("Histogram for ", name1, formatnum(ntst), name2, " p-values")
    r = structure(list(
                breaks = bins,
                counts = cnts,
                density = density,
                mids = centers,
                equidist = FALSE),
            class = "histogram");
    plot(r, main = ntext, ylab = "Density", xlab = "P-values", ...);
    abline(h = 1, col = "blue");
    return(invisible());
}

.qqme = function(m, lcol, cex, pch, ...){
    cnts = m$hist.counts;
    bins = m$hist.bins;
    ntst = m$ntests;
    
    cusu = cumsum(cnts) / ntst;
    ypos = bins[-1][is.finite(cusu)];
    xpos = cusu[is.finite(cusu)];
    lines(-log10(xpos), -log10(ypos), col = lcol, ...);
    if(length(m$eqtls$pvalue)==0)
        return();
    ypvs = -log10(m$eqtls$pvalue);
    xpvs = -log10(1:length(ypvs) / ntst);
    if(length(ypvs) > 1000){
        # need to filter a bit, make the plotting faster
        levels = as.integer( xpvs/xpvs[1] * 1e3);
        keep = c(TRUE, diff(levels)!=0);
        levels = as.integer( ypvs/ypvs[1] * 1e3);
        keep = keep | c(TRUE, diff(levels)!=0);
        ypvs = ypvs[keep];
        xpvs = xpvs[keep];
        rm(keep)
    }
    points(xpvs, ypvs, col = lcol, pch = pch, cex = cex, ...);
}

plot.MatrixEQTL = function(
            x, 
            cex = 0.5,
            pch = 19,
            xlim = NULL,
            ylim = NULL,
            main = NULL, ...){
    if( x$param$pvalue.hist == FALSE ){
        stop("Cannot plot p-value distribution:",
            " the information was not recorded.\n",
            "Use non FALSE pvalue.hist.");
        return(invisible());
    }
    if( x$param$pvalue.hist == "qqplot" ){
        xmin = 1/max(x$cis$ntests, x$all$ntests);
        ymax = NULL;
        if(!is.null(ylim)){
            ymax = ylim[2];
        } else {
            ymax = -log10(min( 
                    x$cis$eqtls$pvalue[1],
                    x$cis$hist.bins[  c(FALSE,x$cis$hist.counts>0)][1],
                    x$all$eqtls$pvalue[1],
                    x$all$hist.bins[  c(FALSE,x$all$hist.counts>0)][1],
                    x$trans$eqtls$pvalue[1],
                    x$trans$hist.bins[c(FALSE,x$trans$hist.counts>0)][1],
                    na.rm = TRUE))+0.1;
        }
        if(ymax == 0){
            ymax = -log10(.Machine$double.xmin);
        }
        if(!is.null(ymax))
            ylim = c(0,ymax);
        
        if(is.null(xlim))
            xlim =  c(0, -log10(xmin/1.5));
        
        plot(numeric(),numeric(), 
            xlab = expression('\u2013 log'[10]*'(p-value), expected'),
            ylab = expression('\u2013 log'[10]*'(p-value), observed'),
            # xlab = "-Log10(p-value), theoretical",
            # ylab = "-Log10(p-value), observed",
            xlim = c(0, -log10(xmin/1.5)),
            ylim = ylim,
            xaxs="i", yaxs="i", ...);
        lines(c(0,1e3), c(0,1e3), col = "gray");
        if( (x$param$pvOutputThreshold > 0) && 
            (x$param$pvOutputThreshold.cis > 0) ){
            .qqme( x$cis, "red", cex, pch, ...);
            .qqme( x$trans, "blue", cex, pch, ...);
            if(is.null(main)) {
                main = paste("QQ-plot for",
                            formatnum(x$cis$ntests), "local and",
                            formatnum(x$trans$ntests), "distant p-values");
            }
            lset = c(1,2,4);
        } else
        if(x$param$pvOutputThreshold.cis > 0){
            .qqme(x$cis, "red", cex, pch, ...);
            if(is.null(main)){
                main = paste("QQ-plot for",
                            formatnum(x$cis$ntests), "local p-values");
            }
            lset = c(1,4);
        } else {
            .qqme(x$all, "blue", cex, pch, ...);
            if(is.null(main)){
                main = paste("QQ-plot for all",
                            formatnum(x$all$ntests), "p-values");
            }
            lset = c(3,4);
        }
        title(main);

        legend("topleft",
            c(  "Local p-values", 
                "Distant p-values",
                "All p-values",
                "diagonal")[lset],
            col =      c("red", "blue", "blue", "gray")[lset],
            text.col = c("red", "blue", "blue", "gray")[lset],
            pch = 20, lwd = 1, pt.cex = c(1,1,1,0)[lset])
    } else {
        if( (x$param$pvOutputThreshold > 0) && 
            (x$param$pvOutputThreshold.cis > 0) ){
            
            par(mfrow = c(2,1));
            .histme(x$cis, "", " local", ...);
            tran = list(hist.counts = x$all$hist.counts - x$cis$hist.counts,
                    hist.bins = x$all$hist.bins,
                    ntests =  x$all$ntests - x$cis$ntests);
            .histme(x$trans,""," distant", ...);
            par(mfrow=c(1,1));
        } else
        if(x$param$pvOutputThreshold.cis > 0){
            .histme(x$cis, "", " local", ...);
        } else {
            .histme(x$all, "all ", ""  , ...);
        }
    }
    return(invisible());
}

