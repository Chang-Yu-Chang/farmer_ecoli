#' Read the
library(tidyverse)
library(data.table)
library(broom)

df_amylase_func <- fread("../data/raw/amylase_summary.csv") %>%
    setNames(c("Transfer", paste0(rep(c("Max", "Error", "Mean", "Sd"), 2), "_", rep(c("expt", "ctrl"), each = 4)))) %>%
    pivot_longer(cols = -Transfer, names_to = c("Stat", "Experiment"), names_pattern = "(.*)_(.*)", values_to = "CommunityFunction")

fwrite(df_amylase_func, "../data/temp/df_amylase_func.txt")

# Mean function
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
    if( equal.variance==FALSE )
    {
        se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        # welch-satterthwaite df
        df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    } else
    {
        # pooled standard deviation, scaled by the sample sizes
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
        df <- n1+n2-2
    }
    t <- (m1-m2-m0)/se
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))
    dat <- as.data.frame(matrix(dat, nrow = 1))
    colnames(dat) <- c("Difference of means", "Std Error", "t", "p.value")
    return(dat)
}

# t test
df_amylase_ttest_mean <- df_amylase_func %>%
    filter(Stat %in% c("Mean", "Sd")) %>%
    group_by(Transfer) %>%
    pivot_wider(names_from = c(Stat, Experiment), values_from = CommunityFunction) %>%
    mutate(N_expt = 24, N_ctrl = 24) %>%
    nest() %>%
    mutate(test = map(data, ~t.test2(m1 = .x$Mean_ctrl, m2 = .x$Mean_expt,
                                     s1 = .x$Sd_ctrl, s2 = .x$Sd_expt,
                                     n1 = .x$N_ctrl, n2 = .x$N_expt))) %>%
    unnest(test) %>%
    select(Transfer, t, p.value)

df_amylase_ttest_mean$Asterisk <- rep(NA, nrow(df_amylase_ttest_mean))
df_amylase_ttest_mean$Asterisk[df_amylase_ttest_mean$p.value >= 0.05] <- ""
df_amylase_ttest_mean$Asterisk[df_amylase_ttest_mean$p.value < 0.05] <- "*"
df_amylase_ttest_mean$Asterisk[df_amylase_ttest_mean$p.value < 0.01] <- "**"
df_amylase_ttest_mean$Asterisk[df_amylase_ttest_mean$p.value < 0.001] <- "***"

fwrite(df_amylase_ttest_mean, "../data/temp/df_amylase_ttest_mean.txt")

## Reshape data
df_amylase_func_mean <- df_amylase_func %>%
    filter(Stat %in% c("Mean", "Sd")) %>%
    pivot_wider(names_from = Stat, values_from = CommunityFunction) %>%
    setNames(c("Transfer", "Experiment", "MeanCommunityFunction", "SdCommunityFunction")) %>%
    arrange(Experiment, Transfer) %>%
    left_join(df_amylase_ttest_mean)

## Asterisk
df_amylase_func_mean_star <- df_amylase_func_mean %>%
    mutate(temp_max = MeanCommunityFunction + SdCommunityFunction) %>%
    select(Transfer, Experiment, temp_max, Asterisk) %>%
    pivot_wider(names_from = Experiment, values_from = temp_max) %>%
    group_by(Transfer) %>%
    mutate(AsteriskLocation = max(ctrl, expt) + 0.02) %>%
    select(Transfer, AsteriskLocation, Asterisk)

fwrite(df_amylase_func_mean, "../data/temp/df_amylase_func_mean.txt")
fwrite(df_amylase_func_mean_star, "../data/temp/df_amylase_func_mean_star.txt")


# Maximum function
## t test
df_amylase_ttest_max <-
    df_amylase_func %>%
    filter(Stat %in% c("Max", "Error")) %>%
    group_by(Transfer) %>%
    pivot_wider(names_from = c(Stat, Experiment), values_from = CommunityFunction) %>%
    mutate(N_expt = 4, N_ctrl = 4) %>%
    nest() %>%
    mutate(test = map(data, ~t.test2(m1 = .x$Max_ctrl, m2 = .x$Max_expt,
                                     s1 = .x$Error_ctrl, s2 = .x$Error_expt,
                                     n1 = .x$N_ctrl, n2 = .x$N_expt))) %>%
    unnest(test) %>%
    select(Transfer, t, p.value)


df_amylase_ttest_max$Asterisk <- rep(NA, nrow(df_amylase_ttest_max))
df_amylase_ttest_max$Asterisk[df_amylase_ttest_max$p.value >= 0.05] <- ""
df_amylase_ttest_max$Asterisk[df_amylase_ttest_max$p.value < 0.05] <- "*"
df_amylase_ttest_max$Asterisk[df_amylase_ttest_max$p.value < 0.01] <- "**"
df_amylase_ttest_max$Asterisk[df_amylase_ttest_max$p.value < 0.001] <- "***"

fwrite(df_amylase_ttest_max, "../data/temp/df_amylase_ttest_max.txt")

## Join ttest result
df_amylase_func_max <- df_amylase_func %>%
    filter(Stat %in% c("Max", "Error")) %>%
    pivot_wider(names_from = Stat, values_from = CommunityFunction) %>%
    setNames(c("Transfer", "Experiment", "MaxFunction", "ErrorFunction")) %>%
    left_join(df_amylase_ttest_max)

# Asterisk
df_amylase_func_max_star <- df_amylase_func_max %>%
    mutate(temp_max = MaxFunction + ErrorFunction) %>%
    select(Transfer, Experiment, temp_max, Asterisk) %>%
    pivot_wider(names_from = Experiment, values_from = temp_max) %>%
    group_by(Transfer) %>%
    mutate(AsteriskLocation = max(ctrl, expt) + 0.02) %>%
    select(Transfer, AsteriskLocation, Asterisk)


fwrite(df_amylase_func_max, "../data/temp/df_amylase_func_max.txt")
fwrite(df_amylase_func_max_star, "../data/temp/df_amylase_func_max_star.txt")



