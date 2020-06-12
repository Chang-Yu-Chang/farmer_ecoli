#' Suppplementary figures and tables

library(tidyverse)
library(data.table)
library(cowplot)
library(knitr)
library(kableExtra)

# Read data
df_amylase_strain <- read_csv("../data/raw/amylase/single_fits.csv")
df_amylase_lichen <- read_csv("../data/raw/amylase/licheniformis_single.csv")
df_farmer_OD <- read_csv("../data/temp/df_farmer_OD.txt")
df_farmer_func <- read_csv("../data/temp/df_farmer_func.txt")

# Table S1
## The data is adapted from Sanchez-Gorostiaga2019
t1 <- df_amylase_strain %>%
    mutate(Vj.fit = round(Vj.fit, 3), se = round(se, 3)) %>%
    unite("Vi", Vj.fit, se, sep = "±") %>%
    filter(S1 %in% c("Subtilis", "Megaterium", "Polymyxa", "licheniformis")) %>%
    mutate(S1 = c("Bacillus subtilis", "Bacillus megaterium", "Paenibacillus polymyxa", "Bacillus licheniformis")) %>%
    setNames(c(" ", "V_{i}(hr^{-1})")) %>%
    kable(format = "latex", booktabs = T, escape = F) %>%
    column_spec(1, italic = T, )

save_kable(t1, file = "../figure/TabS1.png")


# Fig S1
p1 <- ggplot(df_amylase_lichen, aes(t, f.j, color=ordered(x, sort(unique(df_amylase_lichen$x))), group=x)) +
    geom_point() + geom_line(linetype='dashed') +
    theme_cowplot() +
    theme(panel.grid=element_blank()) +
    labs(x='time (hr.)', y='Fraction of starch degraded',
         color = 'SN\ndilution')
ggsave('../figure/FigS1.png', plot = p1, width = 4, height = 3)

# Fig S2
p2_pre <- df_farmer_OD %>%
    # Remove blank wells on DW384
    filter(Media != "B", Inoculation == T) %>%
    # Remove blank wells on DW96
    filter(!(Media %in% paste0("H", sprintf("%02d", 9:12)))) %>%
    filter(Media != "M9g") %>%
    ggplot(aes(x = Transfer, y = Abs, group = Transfer)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_jitter(width = 0.2, size = 0.5) +
    scale_x_continuous(breaks = 1:7) +
    theme_classic() +
    facet_grid(.~Experiment, labeller = as_labeller(c(ctrl="Control", expt = "Selection"))) +
    theme(legend.position = "top", axis.line = element_blank()) +
    labs(x = "Generation", y = expression(OD[620])) +
    panel_border(color = 1)+
    NULL

p2_inset <- df_farmer_OD %>%
    # Remove blank wells on DW384
    filter(Media != "B", Inoculation == T) %>%
    # Remove blank wells on DW96
    filter(!(Media %in% paste0("H", sprintf("%02d", 9:12)))) %>%
    filter(Media == "M9g") %>%
    ggplot(aes(x = Transfer, y = Abs, group = Transfer)) +
    geom_boxplot(outlier.size = 0.2) +
    geom_jitter(width = 0.2, size = 0.2) +
    scale_x_continuous(breaks = 1:7) +
    theme_classic() +
    theme(legend.position = "top", axis.line = element_blank(), axis.title = element_blank()) +
    panel_border(color = 1)+
    NULL

p2 <- plot_grid(p2_pre) +
    draw_plot(p2_inset, x = 0.75, y = 0.55, width = 0.2, height = 0.3)
ggsave("../figure/FigS2.png", plot = p2, width = 8, height = 3)


# Figure S3
df_var <- df_farmer_func %>%
    group_by(Transfer, Experiment) %>%
    summarise(VarCommunityFunction = var(MeanFunction),
              MeanCommunityFunction = mean(MeanFunction),
              SeCommunityFunction = sd(MeanFunction)/sqrt(n()),
              Count = n(),
              CoefficientVar = VarCommunityFunction/MeanCommunityFunction)

ylim.prim <- c(0, 0.4)   # in this example, coeffiecient of variance
ylim.sec <- c(0, 0.4)    # in this example, mean
b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])
mycolors = c(left = "#E7861B", right = "#7997FF")

p3 <- df_var %>%
    ggplot(aes(x = Transfer)) +
    # Var, primary axis
    geom_line(aes(y = CoefficientVar, linetype = Experiment, group = Experiment), color = mycolors["left"]) +
    geom_point(aes(y = CoefficientVar, group = Experiment), shape = 22, fill = "white", size = 3, color = mycolors["left"]) +
    # Mean, secondary axis
    geom_errorbar(aes(x=Transfer, ymin = a + b * (MeanCommunityFunction-SeCommunityFunction),
                      ymax = a + b * (MeanCommunityFunction+SeCommunityFunction)),
                  width = .2, position = position_dodge(0.01), color = mycolors["right"]) +
    geom_line(aes(y = a + b * MeanCommunityFunction, group = Experiment, linetype = Experiment), color = mycolors["right"]) +
    geom_point(aes(y = a + b * MeanCommunityFunction, group = Experiment), shape = 21, fill = "white", size = 3, color = mycolors["right"]) +
    scale_x_continuous(breaks = 1:7) +
    scale_y_continuous("Coefficient of Variance", sec.axis = sec_axis(~ (. - a)/b, name = "Mean")) +
    scale_linetype_manual(values = c(2,1), labels = c("Control", "Selection")) +
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
    labs(x = "Generation") +
     guides(color = F)

ggsave("../figure/FigS3.png", plot = p3, width = 4, height = 3)


# Muller diagram of control line
df_muller_ctrl <- fread("../data/temp/df_muller_ctrl.txt")
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
p4 <- df_muller_ctrl %>%
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

ggsave("../figure/FigS4.png", plot = p4, width = 5, height = 4)



