#' Statistics for the second expeirment: farmer E coli
#' 1. Mean function difference using t test
#' 2. Maximum function differnece using t test
#' 3. Mann-Kendall test between time and maximam function

library(tidyverse)
library(data.table)
library(broom)

df_farmer_func <- fread("../data/temp/df_farmer_func.txt") %>% as_tibble()

# Mean function
## Compute pairwise t test
df_farmer_ttest_mean <-
    df_farmer_func %>%
    select(Transfer, Experiment, Media, MeanFunction) %>%
    group_by(Transfer) %>%
    pivot_wider(names_from = Experiment, values_from = MeanFunction) %>%
    nest() %>%
    mutate(test = map(data, ~t.test(.x$ctrl, .x$expt)),
           tidied = map(test, glance)) %>%
    unnest(tidied) %>%
    select(Transfer, estimate1, estimate2, p.value, conf.low, conf.high, method, alternative)

df_farmer_ttest_mean$Asterisk <- rep(NA, nrow(df_farmer_ttest_mean))
df_farmer_ttest_mean$Asterisk[df_farmer_ttest_mean$p.value >= 0.05] <- ""
df_farmer_ttest_mean$Asterisk[df_farmer_ttest_mean$p.value < 0.05] <- "*"
df_farmer_ttest_mean$Asterisk[df_farmer_ttest_mean$p.value < 0.01] <- "**"
df_farmer_ttest_mean$Asterisk[df_farmer_ttest_mean$p.value < 0.001] <- "***"

fwrite(df_farmer_ttest_mean, "../data/temp/df_farmer_ttest_mean.txt")

## Compute mean
df_farmer_func_mean <- df_farmer_func %>%
    select(Transfer, Experiment, Media, CommunityFunction = MeanFunction) %>%
    group_by(Transfer, Experiment) %>%
    summarize(MeanCommunityFunction = mean(CommunityFunction),
              SdCommunityFunction = sd(CommunityFunction),
              SeCommunityFunction = SdCommunityFunction/sqrt(n())) %>%
    arrange(Experiment, Transfer) %>%
    left_join(df_farmer_ttest_mean)

## Add asterisks
df_farmer_func_mean_star <- df_farmer_func_mean %>%
    mutate(temp_max = MeanCommunityFunction + SeCommunityFunction) %>%
    select(Transfer, Experiment, temp_max, Asterisk) %>%
    pivot_wider(names_from = Experiment, values_from = temp_max) %>%
    mutate(AsteriskLocation = max(ctrl, expt) + 0.01) %>%
    select(Transfer, AsteriskLocation, Asterisk)

fwrite(df_farmer_func_mean, file = "../data/temp/df_farmer_func_mean.txt")
fwrite(df_farmer_func_mean_star, file = "../data/temp/df_farmer_func_mean_star.txt")

# Maximum function
## Compute t test from mean, se, and sample size
t.test2 <- function(m1, m2, s1, s2, n1, n2){
    se <- sqrt((s1^2/n1) + (s2^2/n2))
    degree_freedom <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    t_stat <- (m1-m2)/se
    c(m1-m2, se, t_stat, 2*pt(-abs(t_stat),degree_freedom)) %>%
        matrix(nrow = 1) %>%
        as.data.frame() %>%
        setnames(c("Difference of means", "Std Error", "t", "p.value")) %>%
        return()
}

df_farmer_ttest_max <-
    df_farmer_func %>%
    group_by(Transfer, Experiment) %>%
    filter(MeanFunction == max(MeanFunction)) %>%
    select(Transfer, Experiment, Mean = MeanFunction, Sd = SdFunction) %>%
    group_by(Transfer) %>%
    pivot_wider(names_from = Experiment, values_from = c(Mean, Sd)) %>%
    mutate(N_expt = 3, N_ctrl = 3) %>%
    nest() %>%
    mutate(test = map(data, ~t.test2(m1 = .x$Mean_ctrl, m2 = .x$Mean_expt,
                                     s1 = .x$Sd_ctrl, s2 = .x$Sd_expt,
                                     n1 = .x$N_ctrl, n2 = .x$N_expt))) %>%
    unnest(test) %>%
    select(Transfer, t, p.value)


df_farmer_ttest_max$Asterisk <- rep(NA, nrow(df_farmer_ttest_max))
df_farmer_ttest_max$Asterisk[df_farmer_ttest_max$p.value >= 0.05] <- ""
df_farmer_ttest_max$Asterisk[df_farmer_ttest_max$p.value < 0.05] <- "*"
df_farmer_ttest_max$Asterisk[df_farmer_ttest_max$p.value < 0.01] <- "**"
df_farmer_ttest_max$Asterisk[df_farmer_ttest_max$p.value < 0.001] <- "***"

fwrite(df_farmer_ttest_max, "../data/temp/df_farmer_ttest_max.txt")


## Subset for max
df_farmer_func_max <- df_farmer_func %>%
    group_by(Transfer, Experiment) %>%
    filter(MeanFunction == max(MeanFunction, na.rm = T)) %>%
    arrange(Experiment, Transfer) %>%
    left_join(df_farmer_ttest_max)

## Add asterisks
df_farmer_func_max_star <- df_farmer_func_max %>%
    mutate(temp_max = MeanFunction + SeFunction) %>%
    select(Transfer, Experiment, temp_max, Asterisk) %>%
    pivot_wider(names_from = Experiment, values_from = temp_max) %>%
    mutate(AsteriskLocation = max(ctrl, expt) + 0.01) %>%
    select(Transfer, AsteriskLocation, Asterisk)

fwrite(df_farmer_func_max, file = "../data/temp/df_farmer_func_max.txt")
fwrite(df_farmer_func_max_star, file = "../data/temp/df_farmer_func_max_star.txt")


## Mann-Kendall non-parametric test to test the trend significantly increases or descreases
x <- df_farmer_func_max %>%
    filter(Transfer, Experiment == "expt") %>%
    select(Transfer, Experiment, MeanFunction) %>%
    pull(MeanFunction)

cor.test(x=1:length(x),y=x, meth="kendall", continuity = TRUE)





