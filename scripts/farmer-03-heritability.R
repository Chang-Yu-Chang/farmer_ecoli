#' Heritability for figure 2E
#' 1. Identify the lineages of the two top dominant communiteis
#' 2. Map the parent and offspring functions within these two lineages
#' 3. Compute the heritability per gernation using boostrapped linear regression

library(tidyverse)
library(data.table)
library(broom)

df_farmer_func <- fread("../data/temp/df_farmer_func.txt") %>% as_tibble() %>% mutate(Transfer = factor(Transfer))
df_farmer_func_raw <- fread("../data/temp/df_farmer_func_raw.txt") %>% as_tibble() %>% mutate(Transfer = factor(Transfer))
df_well_map_origin <- fread("../data/temp/df_well_map_origin.txt", na.strings = "") %>% as_tibble

# Identify the parent-offspring pairs within the lineages of the top two communities
df_well_map_top <- df_well_map_origin %>%
    mutate(Media = T1) %>%
    pivot_longer(cols = -Media, names_to = "Transfer", names_pattern = "T(.)", values_to = "Origin") %>%
    filter(Origin %in% c("B01", "D01")) %>% # Top two communities
    mutate(Transfer = factor(Transfer))

# Descendants of the two top communities
df_farmer_func_top <- df_farmer_func %>%
    filter(Experiment == "expt") %>%
    right_join(df_well_map_top) %>%
    select(Transfer, Media, MeanFunction, SdFunction, SeFunction) %>%
    arrange(Transfer, Media)

df_farmer_func_raw_top <- df_farmer_func_raw %>%
    filter(Experiment == "expt") %>%
    right_join(df_well_map_top) %>%
    select(Transfer, Media, Function, MeasurementReplicate) %>%
    arrange(Transfer, Media)


# Map the parent function to the offspring function
wells <- paste0(rep(LETTERS[1:8], each = 4), sprintf("%02d", c(rep(1:4, 8), c(rep(5:8, 8)), c(rep(9:12, 8)))))

df_well_map <- list.files(path = "../data/raw/farmer_mapping_files", pattern = "expt_wells.csv", full.names = T) %>%
    lapply(function(x){
        fread(x, header = T) %>% select(-V1) %>%
            mutate(temp = wells[1:nrow(.)]) %>%
            setNames(c("Parent", "Offspring")) %>%
            select(Parent, Offspring)
    }) %>%
    rbindlist(idcol = "Transfer") %>%
    as_tibble()

temp_list <- rep(list(NA), 6); names(temp_list) <- 1:6
temp_list_raw <- rep(list(NA), 6); names(temp_list_raw) <- 1:6

for (i in 1:length(temp_list)) {
    temp_map <- df_well_map %>%
        filter(Transfer == i) %>%
        select(Parent, Offspring)
    # Mean of measurement
    parent <- df_farmer_func %>%
        filter(Transfer == i, Experiment == "expt") %>%
        select(Media, MeanFunction, SdFunction, SeFunction) %>%
        setNames(c("Parent", "ParentFunction", "ParentFunctionSd", "ParentFunctionSe"))

    offspring <- df_farmer_func %>%
        filter(Transfer == i+1, Experiment == "expt") %>%
        select(Media, MeanFunction, SdFunction, SeFunction) %>%
        setNames(c("Offspring", "OffspringFunction", "OffspringFunctionSd", "OffspringFunctionSe"))

    temp_list[[i]] <-
        offspring %>%
        left_join(temp_map, by = "Offspring") %>%
        filter(!is.na(Parent)) %>%
        left_join(parent, by = "Parent")

    # Raw measurements
    parent <- df_farmer_func_raw %>%
        filter(Transfer == i, Experiment == "expt") %>%
        select(Media, Function, MeasurementReplicate) %>%
        setNames(c("Parent", "ParentFunction", "MeasurementReplicate"))

    offspring <- df_farmer_func_raw %>%
        filter(Transfer == i+1, Experiment == "expt") %>%
        select(Media, Function, MeasurementReplicate) %>%
        setNames(c("Offspring", "OffspringFunction", "MeasurementReplicate"))

    temp_list_raw[[i]] <-
        offspring %>%
        left_join(temp_map, by = "Offspring") %>%
        filter(!is.na(Parent)) %>%
        left_join(parent, by = c("Parent", "MeasurementReplicate"))
}

# Heritablity, aggregated triplicate measurement error
df_farmer_heritability <- rbindlist(temp_list, idcol = "Transfer") %>% as_tibble %>%
    filter(!is.na(ParentFunction)) %>% # Remove dubious communities
    left_join(df_well_map_top %>% setNames(c("Parent", "Transfer", "Origin")), Joining, by = c("Transfer", "Parent"))

# Heritability, with raw replicate mearuements. Used for boostrapping
df_farmer_heritability_raw <- rbindlist(temp_list_raw, idcol = "Transfer") %>% as_tibble %>%
    filter(!is.na(ParentFunction)) %>% # Remove dubious communities
    left_join(df_well_map_top %>% setNames(c("Parent", "Transfer", "Origin")), Joining, by = c("Transfer", "Parent")) %>%
    arrange(Transfer, Parent, Offspring)



# Boostrapping
## Generate bootstrapped data; this section may take a couple minutes
days <-  unique(df_farmer_heritability_raw$Transfer) # Generatations
temp_list <- rep(list(NA), length(days)); names(temp_list) <- 1:length(days)
b = 1000

for (i in 1:length(days)) {
    temp_list_boot <- rep(list(NA), b)
    for (k in 1:b) {
        set.seed(k)
        temp_list_boot[[k]] <-
            df_farmer_heritability_raw %>%
            filter(Transfer == i) %>%
            group_by(Parent, Offspring) %>%
            group_split() %>%
            # For each parent-offspring pair, randomly draw a match from 3x3 measurement matches
            lapply(function(x) slice(x, sample(1:nrow(x), size = 1))) %>%
            bind_rows() %>%
            # For each bootstrap, do a simple linear regression
            lm(formula = OffspringFunction ~ ParentFunction, data = .) %>%
            tidy() %>%
            mutate(term = c("Intercept", "Slope")) %>%
            select(term, estimate, p.value) %>%
            pivot_wider(names_from = term, values_from = c(estimate, p.value))
        cat(k, " ")
    }
    temp_list[[i]] <- rbindlist(temp_list_boot)
}

df_farmer_heritability_boot <- rbindlist(temp_list, idcol = "Transfer")

fwrite(df_farmer_heritability_boot, file = "../data/temp/df_farmer_heritability_boot.txt")

## Estimate the bootstrap p value
df_farmer_heritability_boot <- fread("../data/temp/df_farmer_heritability_boot.txt")

df_farmer_heritability_boot_stat <- df_farmer_heritability_boot %>%
    mutate(PositiveSlope = estimate_Slope > 0, Significant = p.value_Slope < 0.05) %>%
    group_by(Transfer, PositiveSlope, Significant) %>%
    # Count the number of bootstraps that have significant lm fit
    summarize(Count = n()) %>%
    right_join(data.frame(Transfer = rep(1:6, each = 4), PositiveSlope = c(T, F), Significant = c(T,T,F,F))) %>%
    replace_na(list(Count = 0)) %>%
    # Compute the frequency
    mutate(Frequency = Count /1000) %>%
    filter(Significant) %>%
    select(Transfer, Count) %>%
    pivot_wider(names_from = PositiveSlope, values_from = Count) %>%
    setNames(c("Transfer", "PositiveSlope", "NegativeSlope")) %>%
    # Compute p-value of bootstrap by the number of significant lm fits of bootstrap
    mutate(p_boot = 1 - (PositiveSlope + NegativeSlope) / 1000) %>%
    # Mean heritability
    left_join(df_farmer_heritability_boot %>% group_by(Transfer) %>%
                  summarize(MeanSlope = mean(estimate_Slope), MeanIntercept = mean(estimate_Intercept)))

# p value
df_farmer_heritability_boot_stat$p.value_plot <- paste0("p = ", round(df_farmer_heritability_boot_stat$p_boot, 3))
df_farmer_heritability_boot_stat$p.value_plot[df_farmer_heritability_boot_stat$p_boot < 0.001] <- "p < 0.001"

# set non-significant slope to 0
df_farmer_heritability_boot_stat$MeanSlope[df_farmer_heritability_boot_stat$p_boot>0.05] <- 0

fwrite(df_farmer_heritability_boot_stat, file = "../data/temp/df_farmer_heritability_boot_stat.txt")








