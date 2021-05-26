library(ggplot2)
library(readr)
library(dplyr)
library(perm)

data <- read_csv("data/storms_summary.csv") %>%
  mutate(RI = replace(RI, RI == 0, "non-RI")) %>%
  mutate(RI = replace(RI, RI == 1, "RI")) %>%
  mutate(RI = factor(RI, levels = c("RI", "non-RI"), ordered = TRUE))

data %>% ggplot(aes(x = RI, y = shear, fill = RI)) +
  geom_boxplot(outlier.size = 1.1, show.legend = FALSE) +
  theme_classic() +
  labs(y = "Shear Magnitude (m/s)", x = "") +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  theme(axis.text.x = element_text(size = 15, color = "black", vjust = -1),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14))
  

permTS(data %>% filter(RI == "RI") %>% pull(shear),
       data %>% filter(RI == "non-RI") %>% pull(shear),
       alternative = "two.sided",
       method = "exact.mc")

data %>% mutate(RI = factor(RI, levels = c("non-RI", "RI"))) %>%
  ggplot(aes(x = storm, fill = RI)) +
  geom_bar(color = "black", size = 0.25) + scale_fill_manual(values = c("#00BFC4","#F8766D")) +
  labs(x = "Storm #", y = "Number of Observations", fill= "") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 130, 20)) +
  scale_y_continuous(breaks = seq(0, 75, 25)) +
  theme(legend.position = c(0.9, 0.925), legend.direction = "horizontal") 
  

  