#' This script reads and merge the farmer E coli raw data
library(tidyverse)
library(data.table)

# Function for reading the spec output format
read.spec.data.single <- function(f, reader, NUNC96){
    options(stringsAsFactors = FALSE)

    # .. NUNC384 plates
    if(NUNC96 == F){
        if (reader == 1) {
            x <- read.delim(f, skipNul = T, fileEncoding="latin1", header = F)
            colnames(x) <- x[10,]
            x <- x[-c(1:10, nrow(x)),]
            x <- data.table('well' = x$Well,
                            'row' = match(substr(x$Well, 1, 1), LETTERS),
                            'col' = ordered(as.numeric(substr(x$Well, 3, 4))),
                            'wl' = 620,
                            'abs' = x$Abs)
            x <- x[order(x$well),]
            index = NULL
            for(i in 1:24){
                index <- c(index, seq(i, i + 16*23, by = 24))
            }
            x <- x[index,]
        }

        if (reader == 2) {

            # wavelength
            wl <- 620

            # .. col names
            wells <- paste0(rep(LETTERS[1:16], 24),
                            sprintf("%02d", rep(1:24, each=16)))

            # .. read data
            x <- read.delim(f, row.names=NULL, header=F)$V7
            x <- x[-1]
            x <- data.table(x)
            x[,well := wells]
            names(x) <- c("abs", "well")

            # .. update dataframe
            x[, wl := wl]
            x[, row := match(substr(well, 1, 1), LETTERS)]
            x[, col := as.numeric(substr(x$well, 2, 3))]
            setcolorder(x, c('well', 'row', 'col', 'wl', 'abs'))
        }
    }


    # .. NUNC96 plates
    if(NUNC96 == T){
        if (reader == 1) {
            x <- read.delim(f, skipNul = T, fileEncoding="latin1", header = F)
            colnames(x) <- x[10,]
            x <- x[-c(1:10, nrow(x)),]
            x <- data.table('well' = x$Well,
                            'row' = match(substr(x$Well, 1, 1), LETTERS),
                            'col' = ordered(as.numeric(substr(x$Well, 3, 4))),
                            'wl' = 620,
                            'abs' = x$Abs)
            x <- x[order(x$well),]
            index = NULL
            for(i in 1:12){
                index <- c(index, seq(i, i + 8*11, by = 12))
            }
            x <- x[index,]
        }

        if (reader == 2) {

            # .. wavelength
            wl <- 620

            # .. col names
            wells <- paste0(rep(LETTERS[1:8], 12),
                            sprintf("%02d", rep(1:12, each=8)))

            # .. read data
            x <- read.delim(f, row.names=NULL, header=F)$V7
            x <- x[-1]
            x <- data.table(x)
            x[,well := wells]
            names(x) <- c("abs", "well")

            # .. update dataframe
            x[, wl := wl]
            x[, row := match(substr(well, 1, 1), LETTERS)]
            x[, col := as.numeric(substr(x$well, 2, 3))]
            setcolorder(x, c('well', 'row', 'col', 'wl', 'abs'))
        }
    }

    return(x)
}

# Read data
folder <- "../data/raw/farmer_OD/"
file_names <- list.files(folder)
temp_list <- rep(list(NA), length(file_names))

for(i in 1:length(temp_list)){
    # T2 uses different spec
    if (grepl("T2", file_names[i])) {
        d <- read.spec.data.single(paste(folder, file_names[i], sep=''), reader = 1, NUNC96 = F)
    } else {
        d <- read.spec.data.single(paste(folder, file_names[i], sep=''), reader = 2, NUNC96 = F)
    }

    # Map wells from original community DW96 to DW384 plates
    map <- read.table(paste0("../data/raw/farmer_mapping_files/", substr(file_names[i], 8, 9),"_screen_384_layout.txt"))
    map <- data.table(map)[,.(well, media,inoc)]
    d <- merge(d, map, by='well', all=TRUE)

    # Subtract the OD by the median value of blanks
    d[, abs := as.numeric(abs)]
    subtract.blanks <- function(x){
        x$abs <- x$abs - median(x$abs[!x$inoc])
        x$abs[x$abs<0] <- 0
        x
    }

    d <- subtract.blanks(d)
    temp_list[[i]] <- cbind(transfer = sub("Farmer_", "", file_names[i]) %>% sub(".txt", "", .), d)
}

temp_list <- rbindlist(temp_list)

# Clean up the data frame
temp_list %>%
    setNames(c("Treatment", "Well", "Row", "Col", "Wavelength", "Abs", "Media", "Inoculation")) %>%
    separate(col = Treatment, into = c("Transfer", "temp", "Experiment")) %>%
    select(-temp) %>%
    extract(col = Transfer, into = c("temp", "Transfer"), regex = "(\\w)(\\d)") %>%
    select(-temp) %>%
    ungroup() %>%
    # Move one generation forward. Old transfer 1 was cycloheximide treatment. Old transfer 2 becomes communtiy generation 1
    mutate(Transfer = as.numeric(as.character(Transfer)) - 1,
           Transfer = factor(Transfer)) %>%
    fwrite('../data/temp/df_farmer_OD.txt')

