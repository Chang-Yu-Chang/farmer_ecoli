#' This script processes the raw data for the second experiment: farmer E coli
#' 1. Generate parent-offspring origin mapping files for muller diagram used in Figure 2F
#' 2. Remove OD data of dubious communities that originates from a cross-contamination
#' 3. Remove replicate with bubble-caused high OD measures

library(tidyverse)
library(data.table)

# Muller
# Read the raw map files
wells <- c(paste0(rep(LETTERS[1:8], 12), sprintf("%02d", rep(1:12, each = 8))))
wells <- wells[! wells %in% paste0("H", sprintf("%02d", 9:12))]
replicate_well <- fread("../data/raw/farmer_mapping_files/T3_ctrl_wells.csv", header = T) %>% pull(destination)

df_wells <- as.data.frame(matrix(NA, 92, 7)) %>% setNames(paste0("T", 1:7))
df_wells$T1 <- wells

# T1-T2; old T2-T3
temp_wells <- fread(paste0("../data/raw/farmer_mapping_files/T3_expt_wells.csv"), header = T) %>%
    slice(1:92) %>%
    select(Parent = source) %>%
    mutate(Offspring = c(replicate_well)) %>%
    as_tibble %>%
    arrange(factor(Offspring, wells)) %>%
    pull(Parent)

df_wells$T2 <- temp_wells

for (i in 2:6) {
    temp_map <- tibble(Origin = temp_wells, Plate = wells)

    temp_wells <-
        fread(paste0("../data/raw/farmer_mapping_files/T", i+2, "_expt_wells.csv"), header = T) %>% as_tibble %>%
        slice(1:92) %>%
        select(Plate = source) %>%
        left_join(temp_map) %>%
        select(Parent = Origin) %>%
        mutate(Offspring = c(replicate_well)) %>%
        as_tibble %>%
        arrange(factor(Offspring, wells)) %>%
        pull(Parent)

    df_wells[paste0("T", i+1)] <- temp_wells
}

df_wells2 <-
    df_wells %>%
    mutate(Well = T1) %>%
    pivot_longer(cols = -Well, names_to = "Transfer", values_to = "Community") %>%
    arrange(Transfer, Well) %>%
    # Remove NA that accidently selects the contanminated blanks
    filter(!is.na(Community)) %>%
    select(Transfer, Community) %>%
    group_by(Transfer, Community) %>%
    summarize(Count = n()) %>%
    ungroup() %>%
    group_by(Transfer) %>%
    mutate(TotalCount = sum(Count),
           Frequency = Count / TotalCount) %>%
    # Add count = 0 back
    right_join(tibble(Transfer = paste0("T", rep(1:7, each = 92)), Community = rep(wells, 7))) %>%
    replace_na(replace = list(Frequency = 0, Count = 0))

fwrite(df_wells, "../data/temp/df_well_map_origin.txt")
fwrite(df_wells2, "../data/temp/df_muller.txt")


# OD data
# Read farmer data
df_farmer <- fread("../data/raw/farmer_OD.txt") %>% as_tibble
df_well_map_origin <- fread("../data/temp/df_well_map_origin.txt", na.strings = "") %>% as_tibble

# Remove the dubious origins communities
temp_list <- df_farmer %>% split.data.frame(f = df_farmer$Transfer)

for (i in 2:length(temp_list)) {
    # Find dubious communities
    temp <- df_well_map_origin %>%
        select(T1, paste0("T", i)) %>%
        setNames(c("Media", "temp")) %>%
        filter(is.na(temp)) %>%
        select(Media) %>%
        mutate(Transfer = i, Experiment = "expt")
    # Filter out the dubious communtieis from the OD data
    temp_list[[i]] <- df_farmer %>%
        filter(Transfer == i) %>%
        anti_join(temp)
}

df_farmer <- rbindlist(temp_list) %>%
    mutate(Transfer = factor(Transfer)) %>%
    as_tibble

# Compute the community function
df_farmer_positive <- df_farmer %>%
    filter(Media == "M9g", Inoculation == T) %>%
    select(Transfer, Experiment, Media, Abs) %>%
    group_by(Transfer, Experiment) %>%
    summarize(MeanAbs_positive = mean(Abs))

# Remove bubble
df_bubble <- df_farmer %>%
    filter(!(Media %in% c("B", "M9g")), Inoculation == T) %>%
    select(Transfer, Experiment, Media, Abs) %>%
    group_by(Transfer, Experiment, Media) %>%
    arrange(Transfer, Experiment, Media) %>%
    summarise(MeanAbs = mean(Abs), SdAbs = sd(Abs)) %>%
    # Triplicates with high SE potentially have bubble when measured
    filter(SdAbs >= 0.05) %>%
    ungroup %>%
    select(Transfer, Experiment, Media) %>%
    mutate(Bubble = T)

df_farmer2 <- df_farmer %>%
    filter(!(Media %in% c("H09", "H10", "H11", "H12"))) %>% # Blank control on DW96
    filter(!(Media %in% c("B", "M9g")), Inoculation == T) # Controls on DW384

df_removed_bubble <- df_farmer2 %>%
    left_join(df_bubble) %>%
    filter(Bubble == T) %>%
    group_by(Transfer, Experiment, Media) %>%
    filter(Abs != max(Abs))

df_farmer_func <-
    df_farmer2 %>%
    # Filter out possible bubble triplicate
    left_join(df_bubble) %>% replace_na(replace = list(Bubble = FALSE)) %>%
    filter(Bubble == F) %>% select(-Bubble) %>%
    # Add back the removed bubble duplicates
    bind_rows(df_removed_bubble) %>%
    # Function = OD normalized by positive control
    left_join(df_farmer_positive) %>%
    mutate(Function = Abs / MeanAbs_positive) %>%
    # Mean functon of triplicates of each community
    select(Transfer, Experiment, Media, Function) %>%
    group_by(Transfer, Experiment, Media) %>%
    arrange(Transfer, Experiment, Media) %>%
    summarise(MeanFunction = mean(Function), SdFunction = sd(Function)) %>%
    {.}


fwrite(df_farmer_func, file = "../data/temp/df_farmer_func.txt")


