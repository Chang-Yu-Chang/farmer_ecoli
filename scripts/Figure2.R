#' This script makes figure from the farmer E coli

library(tidyverse)
library(data.table)
library(cowplot)

df_farmer_func <- fread("../data/temp/df_farmer_func.txt") %>% as_tibble
df_farmer_func_mean <- fread("../data/temp/df_farmer_func_mean.txt") %>% as_tibble
df_farmer_func_mean_star <- fread("../data/temp/df_farmer_func_mean_star.txt") %>% as_tibble
df_farmer_func_max <- fread("../data/temp/df_farmer_func_max.txt") %>% as_tibble
df_farmer_func_max_star <- fread("../data/temp/df_farmer_func_max_star.txt") %>% as_tibble
df_muller <- fread("../data/temp/df_muller.txt") %>% as_tibble()
df_farmer_heritability <- fread("../data/temp/df_farmer_heritability.txt")
df_farmer_heritability_boot_stat <- fread("../data/temp/df_farmer_heritability_boot_stat.txt")


# Panel A. Cartoon
p_cartoon <- ggdraw() + draw_image("../data/experimental_scheme/Figure2A.png")


# Panel B. Mean function
p_B <- df_farmer_func_mean %>%
    ggplot(aes(x = Transfer, y = MeanCommunityFunction, color = Experiment, group = Experiment)) +
    geom_errorbar(aes(x=Transfer, ymin = MeanCommunityFunction-SdCommunityFunction,
                      ymax = MeanCommunityFunction+SdCommunityFunction),
                  width = .2, position = position_dodge(0.3)) +
    geom_line(linetype = 2,position = position_dodge(0.3)) +
    geom_point(shape = 21, fill = "white", size = 3, position = position_dodge(0.3)) +
    # Add asterisk
    geom_text(data = df_farmer_func_mean_star, inherit.aes = F, aes(x = Transfer, y = AsteriskLocation, label = Asterisk), size = 5, color = "black") +
    coord_cartesian(ylim=c(0, 0.61)) +
    scale_x_continuous(breaks = 1:7) +
    scale_color_discrete(labels = c("control", "selection")) +
    theme_cowplot() +
    theme(legend.position = c(0.1, 0.95), legend.title = element_blank()) +
    labs(x = "Generation", y = "Mean Function")


# Panel C. Maximum function
p_C <- df_farmer_func_max %>%
    ggplot(aes(x = Transfer, y = MeanFunction, color = Experiment, group = Experiment)) +
    geom_errorbar(aes(x=Transfer, ymin = MeanFunction-SeFunction, ymax = MeanFunction+SeFunction),
                  width = .2, position = position_dodge(0.3)) +
    geom_line(linetype = 2, position = position_dodge(0.3)) +
    geom_point(shape = 24, fill = "white", size = 3, position = position_dodge(0.3)) +
    # Add asterisk
    geom_text(data = df_farmer_func_max_star, inherit.aes = F, aes(x = Transfer, y = AsteriskLocation, label = Asterisk), size = 5, color = "black") +
    scale_x_continuous(breaks = 1:7) +
    coord_cartesian(ylim=c(0, 0.61)) +
    theme_cowplot() +
    theme(legend.position = "none", legend.title = element_blank()) +
    labs(x = "Generation", y = "Maximum Function")

# Panel D. Histogram
p_D <- df_farmer_func %>%
    ggplot() +
    geom_histogram(aes(x = MeanFunction, y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..], fill = Experiment),
                   color = 1, alpha = 0.6, binwidth = 0.02, position = position_identity()) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_discrete(labels = c("control", "selection")) +
    # Facet strips
    geom_text(aes(label = paste0("Generation ", Transfer)), x = Inf, y = Inf, hjust = 1.2, vjust = 1.5) +
    facet_grid(Transfer~.,
               labeller = labeller(Transfer = paste0("Generation ", 1:7) %>% setNames(1:7))) +
    theme_bw() +
    theme(legend.position ="none", legend.title = element_blank(), legend.background = element_blank(),
          panel.grid = element_blank(), panel.border = element_blank(),
          strip.background = element_blank(), strip.text = element_blank()) +
    labs(x = "Function", y = "Frequency") +
    panel_border(color = "black")

# Panel E. Heritability
p_E <- df_farmer_heritability %>%
    ggplot(aes(x = ParentFunction, y = OffspringFunction)) +
    #    geom_abline(slope = 1, intercept = 0, color = "grey", linetype = 2, size = .5) +
    geom_abline(data = df_farmer_heritability_boot_stat %>% filter(Transfer <= 2),
                color = "grey", linetype = 2,
                aes(slope = MeanSlope, intercept = MeanIntercept)) +
    geom_errorbar(aes(color = Origin, ymin = OffspringFunction-OffspringFunctionSe, ymax = OffspringFunction+OffspringFunctionSe)) +
    geom_errorbarh(aes(color = Origin, xmin = ParentFunction-ParentFunctionSe, xmax = ParentFunction+ParentFunctionSe)) +
    geom_point(aes(color = Origin, shape = Origin), fill = "white") +
    scale_color_manual(values = c("black", "#F3949B", "#E4211D")) +
    scale_shape_manual(values = c(21, 24, 22)) +
    scale_x_continuous(breaks = c(0.1, 0.5)) +
    scale_y_continuous(breaks = c(0.1, 0.5)) +
    coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 0.6)) +
    # Strip
    geom_text(aes(label = paste0("Generation ", Transfer, "-", Transfer+1)), x = 0.3, y = Inf, vjust = 1.5)+ #, hjust = 1.5, vjust = 1.5) +
    #    geom_text(data = , x = 0.3, y = 0.5, paste0("Day", 1:6)) +
    # Add stat
    geom_text(data = df_farmer_heritability_boot_stat %>% filter(Transfer == 1),
              aes(label = paste0("h^2==", round(MeanSlope, 3))), size = 3,
              x = -Inf, y = 0.45, hjust = -0.3, parse = T) +
    geom_text(data = df_farmer_heritability_boot_stat %>% filter(Transfer == 1),
              aes(label = paste0(p.value_plot)), size = 3,
              x = -Inf, y = 0.39, hjust = -0.3) +
    geom_text(data = df_farmer_heritability_boot_stat %>% filter(Transfer == 2),
              aes(label = paste0("h^2==", round(MeanSlope, 3))), size = 3,
              x = 0.3, y = 0.25, hjust = -0.1, parse = T) +
    geom_text(data = df_farmer_heritability_boot_stat %>% filter(Transfer == 2),
              aes(label = paste0(p.value_plot)), size = 3,
              x = 0.3, y = 0.19, hjust = -0.1, parse = T) +
    geom_text(data = df_farmer_heritability_boot_stat %>% filter(Transfer >= 3),
              aes(label = paste0(p.value_plot)), size = 3,
              x = -Inf, y = -Inf, hjust = -.3, vjust = -1) +
    facet_wrap(Transfer~., labeller = labeller(Transfer = paste0("Generation ", 1:6, " - ", 2:7) %>% setNames(1:6))) +
    theme_bw() +
    theme(legend.position = "none",
          panel.border = element_blank(), panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()) +
    panel_border(color = "black", size = 1) +
    labs(x = "Parent", y = "Offspring")


# Panel F. Muller diagram
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
p_F <- df_muller %>%
    extract(col = Community, into = c("Replicate", "Sample"), regex = "(\\w)(\\d+)", remove = F) %>%
    extract(col = Transfer, into = "Transfer", regex = "T(\\d+)", remove = T) %>%
    mutate(Transfer = as.numeric(Transfer)) %>%
    ggplot(aes(x = Transfer, y = Frequency, fill = Sample, alpha = Replicate)) +
    geom_area(color = "black", size = .1) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values =  getPalette(12)) +
    scale_alpha_manual(values = seq(0.5, 1, length.out = 8)) +
    theme_cowplot() +
    theme(legend.position = "right") +
    guides(alpha = F) +
    panel_border(color = "black") +
    labs(x = "Generation", y = "Frequency")

# Panel G. Heritability versus mean
df_G <- df_farmer_func_mean %>%
    filter(Experiment == "expt") %>%
    bind_rows(df_farmer_heritability_boot_stat %>% mutate(Transfer = Transfer + 0.5))

ylim.prim <- c(0, 1)   # in this example, heritability
ylim.sec <- c(0, 0.4)    # in this example, mean
b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])
mycolors = c(left = "#E7861B", right = "#7997FF") # #00BFC4

p_G <- df_G %>%
    ggplot(aes(x = Transfer)) +
    # Heritability, primary axis
    geom_line(aes(y = MeanSlope, group = Experiment), linetype = 2, color = mycolors["left"]) +
    geom_point(aes(y = MeanSlope, group = Experiment), shape = 22, fill = "white", size = 3, color = mycolors["left"]) +
    # Mean, secondary axis
    geom_errorbar(aes(x=Transfer, ymin = a + b * (MeanCommunityFunction-SeCommunityFunction),
                      ymax = a + b * (MeanCommunityFunction+SeCommunityFunction)),
                  width = .2, position = position_dodge(0.01), color = mycolors["right"]) +
    geom_line(aes(y = a + b * MeanCommunityFunction, group = Experiment), linetype = 2, color = mycolors["right"]) +
    geom_point(aes(y = a + b * MeanCommunityFunction, group = Experiment), shape = 21, fill = "white", size = 3, color = mycolors["right"]) +
    scale_x_continuous(breaks = 1:7) +
    scale_y_continuous("Heritability", sec.axis = sec_axis(~ (. - a)/b, name = "Mean")) +
    theme_cowplot() +
    theme(legend.position = "top", legend.title = element_blank(),
          axis.title.y = element_text(color = mycolors["left"]),
          axis.text.y = element_text(color = mycolors["left"]),
          axis.line.y = element_line(color = mycolors["left"]),
          axis.ticks.y = element_line(color = mycolors["left"]),
          axis.title.y.right = element_text(color = mycolors["right"]),
          axis.text.y.right = element_text(color = mycolors["right"]),
          axis.line.y.right = element_line(color = mycolors["right"]),
          axis.ticks.y.right = element_line(color = mycolors["right"])
    ) +
    labs(x = "Generation")


# Put together panels
p_BC <- plot_grid(plotlist = list(p_B, p_C), nrow = 1, labels = c("B", "C"), align = "hv", greedy = T)
p_EFG <- plot_grid(plotlist = list(p_E, p_F, p_G), nrow = 1, labels = c("E", "F", "G"), rel_widths = c(4,3,2.5))
p_left_column <- plot_grid(plotlist = list(p_cartoon, NULL, p_BC, NULL, p_EFG), ncol = 1, labels = c("A", "", ""), rel_heights = c(3.5, 0.3, 3, 0.5, 3))
p <- plot_grid(p_left_column, p_D, ncol = 2, labels = c("", "D"), rel_widths = c(6, 1.5))

ggsave("../figure/Fig2.png", plot = p, width = 16, height = 12)
ggsave("../figure/Fig2.pdf", plot = p, width = 16, height = 12)

