#' Heritability for figure 2E
#' 1. Identify the lineages of the two top dominant communiteis
#' 2. Map the parent and offspring functions within these two lineages
#' 3. Compute the heritability from linear regression

library(tidyverse)
library(data.table)
library(broom)

df_farmer_func <- fread("../data/temp/df_farmer_func.txt") %>% as_tibble() %>% mutate(Transfer = factor(Transfer))
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

temp_list <- rep(list(NA), 6)
names(temp_list) <- 1:6

for (i in 1:length(temp_list)) {
    temp_map <- df_well_map %>%
        filter(Transfer == i) %>%
        select(Parent, Offspring)

    parent <- df_farmer_func %>%
        filter(Transfer == i, Experiment == "expt") %>%
        select(Media, MeanFunction, SeFunction) %>%
        setNames(c("Parent", "ParentFunction", "ParentFunctionSe"))

    offspring <- df_farmer_func %>%
        filter(Transfer == i+1, Experiment == "expt") %>%
        select(Media, MeanFunction, SeFunction) %>%
        setNames(c("Offspring", "OffspringFunction", "OffspringFunctionSe"))

    temp_list[[i]] <-
        offspring %>%
        left_join(temp_map, by = "Offspring") %>%
        filter(!is.na(Parent)) %>%
        left_join(parent, by = "Parent")
#        select(ParentFunction, ParentFunctionSd, OffspringFunction, OffspringFunctionSd) %>%
}

df_farmer_heritability <- rbindlist(temp_list, idcol = "Transfer") %>% as_tibble %>%
    filter(!is.na(ParentFunction)) %>% # Remove dubious communities
    left_join(df_well_map_top %>% setNames(c("Parent", "Transfer", "Origin")), Joining, by = c("Transfer", "Parent"))


# Compute heritability from slope of parent-offspring pairs
## Overall
df_farmer_heritability_stat_collapsed <- df_farmer_heritability %>%
    lm(data = ., formula = OffspringFunction ~ ParentFunction) %>%
    tidy() %>%
    {.}

df_farmer_heritability_stat_collapsed$p.value_plot <- paste0("p-value = ", round(df_farmer_heritability_stat_collapsed$p.value, 3))
df_farmer_heritability_stat_collapsed$p.value_plot[df_farmer_heritability_stat_collapsed$p.value < 0.001] <- "p-value < 0.001"

## Each generation
df_farmer_heritability_stat <-
    df_farmer_heritability %>%
    group_by(Transfer) %>%
    nest() %>%
    mutate(test = map(data, ~lm(.x$OffspringFunction ~ .x$ParentFunction)),
           tidied = map(test, tidy)) %>%
    unnest(tidied) %>%
    select(Transfer, estimate, p.value) %>%
    group_by(Transfer) %>%
    mutate(temp = c("Intercept", "Slope")) %>%
    pivot_wider(names_from = temp, values_from = c(estimate, p.value))


df_farmer_heritability_stat$p.value_plot <- paste0("p-value = ", round(df_farmer_heritability_stat$p.value_Slope, 3))
df_farmer_heritability_stat$p.value_plot[df_farmer_heritability_stat$p.value_Slope < 0.001] <- "p-value < 0.001"

fwrite(df_farmer_heritability, "../data/temp/df_farmer_heritability.txt")
fwrite(df_farmer_heritability_stat, "../data/temp/df_farmer_heritability_stat.txt")
fwrite(df_farmer_heritability_stat_collapsed, "../data/temp/df_farmer_heritability_stat_collapsed.txt")



# Compute heritability from breeder's equation
df_farmer_func <- fread("../data/temp/df_farmer_func.txt") %>% as_tibble() %>% mutate(Transfer = factor(Transfer))
df_selected_parent <- df_well_map %>% distinct(Transfer, Parent) %>% mutate(Transfer = factor(Transfer))


df_selected_parent_func <- df_farmer_func %>%
    filter(Experiment == "expt") %>%
    select(Transfer, Media, MeanFunction) %>%
    setNames(c("Transfer", "Parent", "Function")) %>%
    right_join(df_selected_parent)


parents <- df_farmer_func %>%
    filter(Experiment == "expt") %>%
    group_by(Transfer) %>%
    summarize(MeanParent = mean(MeanFunction, na.rm = T))

selected_parents <- df_selected_parent_func %>%
    group_by(Transfer) %>%
    summarize(MeanSelectedParent = mean(Function, na.rm = T))

selected_parents %>%
    mutate(MeanParent = parents$MeanParent[1:6],
           MeanOffspring = parents$MeanParent[2:7]) %>%
    mutate(RealizedHeritability = (MeanParent-MeanOffspring)/(MeanParent-MeanSelectedParent))
















