#' Suppplementary figures and tables

library(tidyverse)
library(data.table)
library(cowplot)
library(knitr)
library(kableExtra)

# Read data
df_farmer_OD <- read_csv("../data/temp/df_farmer_OD.txt")
df_farmer_func <- read_csv("../data/temp/df_farmer_func.txt")

# Table S1
## The data is adapted from Sanchez-Gorostiaga2019
t1 <- tibble(
    ` ` = "V_{i}(hr^{-1})",
    `Bacillus subtilis` = "V_{S}=1.21±0.13",
    `Bacillus megaterium` = "V_{M}=0.04±0.01",
    `Paenibacillus polymyxa` = "V_{P}=0.19±0.02") %>%
    kable(format = "latex", booktabs = T, escape = F) %>%
    column_spec(1, italic = T) %>%
    row_spec(0, italic = T)

save_kable(t1, file = "../figure/TabS1.png")


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
    theme(legend.position = "top", axis.line = element_blank()) +
    labs(x = "", y = "") +
    panel_border(color = 1)+
    NULL

p2 <- plot_grid(p2_pre) +
    draw_plot(p2_inset, x = 0.75, y = 0.55, width = 0.2, height = 0.35)
ggsave("../figure/FigS2.png", plot = p2, width = 10, height = 4)




# Coefficient of Variance
df_var <- df_farmer_func %>%
    #filter(Experiment == "expt") %>%
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

ggsave("../figure/FigS3.png", plot = p3, width = 5, height = 5)

