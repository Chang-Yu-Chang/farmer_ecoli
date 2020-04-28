
# R function to read data from plate reader

"
Description
===========
This function uses the data.table package to read and format the data output
from the two platereaders in the lab. In the platereader 1, data should be
saved as list, not table. In platereader 2, data should be saved in txt. 

Arguments
=========

- f        path to the file containing the data (from the current directory)

- reader   which reader does the data come from. 1 for the old and 2 for the
           new one.
" 

read.spec.data <- function(f, reader)    
    {
        require(data.table)
        options(stringsAsFactors = FALSE)
        if (reader == 1) {
            x <- read.delim(f, row.names=NULL)
            x <- data.table('plate' = x$Plate,
                            'well' = x$Well,
                            'row' = match(substr(x$Well, 2, 2), LETTERS),
                            'col' = ordered(as.numeric(substr(x$Well, 3, 4))),
                            'wl' = x$Wavelength,
                            'sample' = x$Sample,
                            't' = round((x$Meas..Time..s.)/60, 0),
                            'abs' = x$Abs)
            x[,well := gsub(' ', '', well)]  
        }

        if (reader == 2) {

            nread <- 1
            aux <- readLines(f, n=nread)
            
            while (!any(grep('avg.time', aux))) {
                nread <- nread+1
                aux <- readLines(f, n=nread)
            }

            # wavelength
            wl <- grep('Wavelength: ', aux)
            wl <- as.numeric(gsub('\\D+','', aux[wl]))

            # plate name, usually in line 2 
            plate <- gsub('\t', '', aux[2])

            # .. col names 
            cn <- unlist(strsplit(aux[nread], '\t'))[-(1:2)]
            wells <- gsub('.*[(]([^.]+)[)].*', '\\1', cn)

            # .. read data            
            x <- read.delim(f, row.names=NULL, skip=nread, header=FALSE)
            names(x) <- c('read', 't', wells)
            x <- x[,-which(is.na(colnames(x)))]
            x <- x[-which(is.na(x$read)),]
            x <- data.table(x)
            
            # .. melt the data to list format 
            x <- melt(x, id.vars = c('read', 't'),
                      variable.name = "well", value.name = "abs",
                      variable.factor = FALSE)
            x[, t := round(t/60)]
            x[, read :=NULL]
            
            # .. update dataframe
            x[, wl := wl]
            x[, plate := plate]
            x[, row := match(substr(well, 1, 1), LETTERS)]
            x[, col := as.numeric(substr(x$well, 2, 3))]
        }
        return(x)
    }



