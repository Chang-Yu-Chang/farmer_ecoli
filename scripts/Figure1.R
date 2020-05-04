#' Code to make the figure plotting the mean and maximum functions of amylolytic communities

library(tidyverse)
library(data.table)
library(cowplot)

df_amylase_func <- fread("../data/temp/df_amylase_func.txt")
df_amylase_func_mean <- fread("../data/temp/df_amylase_func_mean.txt")
df_amylase_func_mean_star <- fread("../data/temp/df_amylase_func_mean_star.txt")
df_amylase_func_max <- fread("../data/temp/df_amylase_func_max.txt")
df_amylase_func_max_star <- fread("../data/temp/df_amylase_func_max_star.txt")

# Panel A. Cartoon
p_cartoon <- ggdraw() + draw_image("../data/experimental_scheme/Figure1A.png")

# Panel B. Mean
p_B <- df_amylase_func_mean %>%
    ggplot(aes(x = Transfer, y = MeanCommunityFunction, color = Experiment, group = Experiment)) +
    geom_errorbar(aes(x=Transfer, ymin = MeanCommunityFunction - SdCommunityFunction,
                      ymax = MeanCommunityFunction + SdCommunityFunction),
                  width = .2, position = position_dodge(0.3)) +
    geom_line(linetype = 2, position = position_dodge(0.3)) +
    geom_point(shape = 21, fill = "white", size = 3, position = position_dodge(0.3)) +
    # Add asterisk
    geom_text(data = df_amylase_func_mean_star, inherit.aes = F, aes(x = Transfer, y = AsteriskLocation, label = Asterisk),
              size = 4, color = "black") +
    scale_x_continuous(breaks = 1:17) +
    scale_color_discrete(labels = c("control", "selection")) +
    coord_cartesian(ylim=c(0, 0.45)) +
    theme_cowplot() +
    theme(legend.position = c(0.1, 0.95), legend.direction = "vertical", legend.title = element_blank()) +
    labs(x = "Generation", y = "Mean Function")


p_B_inset <- df_amylase_func %>%
    filter(Stat %in% c("Mean", "Sd")) %>%
    pivot_wider(names_from = Stat, values_from = CommunityFunction) %>%
    setNames(c("Transfer", "Experiment", "MeanCommunityFunction", "SdCommunityFunction")) %>%
    filter(Transfer >= 4) %>%
    ggplot(aes(x = Transfer, y = MeanCommunityFunction, color = Experiment, group = Experiment)) +
    geom_smooth(method = "lm", lwd = .5, se = T, aes(fill = Experiment)) +
    geom_errorbar(aes(x=Transfer, ymin = MeanCommunityFunction-SdCommunityFunction, ymax = MeanCommunityFunction+SdCommunityFunction),
                  width = .2, position = position_dodge(0.3)) +
    geom_point(shape = 21, fill = "white", size = 1, position = position_dodge(0.3)) +
    scale_x_continuous(breaks = 4:17, labels = c("", 5, rep("", 4), 10, rep("", 4), 15, rep("", 2))) +
    coord_cartesian(ylim=c(0.08, 0.25)) +
    theme_cowplot() +
    theme(legend.position = "none", legend.title = element_blank(), axis.title = element_blank()) +
    labs(x = "Generation", y = "Mean Function")


p_C <- df_amylase_func_max %>%
    ggplot(aes(x = Transfer, y = MaxFunction, color = Experiment, group = Experiment)) +
    geom_errorbar(aes(x=Transfer, ymin = MaxFunction-ErrorFunction, ymax = MaxFunction+ErrorFunction),
                  width = .2, position = position_dodge(0.3)) +
    geom_line(linetype = 2, position = position_dodge(0.3)) +
    geom_point(shape = 24, fill = "white", size = 3, position = position_dodge(0.3)) +
    # Add asterisk
    geom_text(data = df_amylase_func_max_star, inherit.aes = F, aes(x = Transfer, y = AsteriskLocation, label = Asterisk),
              size = 4, color = "black") +
    scale_x_continuous(breaks = 1:17) +
    coord_cartesian(ylim=c(0, 0.45)) +
    theme_cowplot() +
    theme(legend.position = "none", legend.title = element_blank()) +
    labs(x = "Generation", y = "Maximum Function")


# Put together the panels
p_bottom_row <- plot_grid(plotlist = list(p_B, p_C), ncol = 2, labels = LETTERS[2:3], axis = "t", align = "hv") +
    draw_plot(p_B_inset, x = 0.25, y = 0.6, width = 0.25, height = 0.4)

p <- plot_grid(p_cartoon, p_bottom_row, ncol = 1, labels = c("A", ""), rel_heights = c(1, 1))

ggsave("../figure/Fig1.png", plot = p, width = 9, height = 7)
ggsave("../figure/Fig1.pdf", plot = p, width = 9, height = 7)
