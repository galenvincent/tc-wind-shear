library(readr)
library(plyr)
library(tidyr)
library(Hmisc)
library(perm)
library(lubridate)
library(dplyr)
library(ggplot2)


data <- read_csv("data/storms_summary.csv") %>%
  mutate(RI = replace(RI, RI == 0, "non-RI")) %>%
  mutate(RI = replace(RI, RI == 1, "RI")) %>%
  mutate(RI = factor(RI, levels = c("RI", "non-RI"), ordered = TRUE))

data %>% ggplot(aes(x = RI, y = shear, fill = RI)) +
  geom_boxplot(outlier.size = 1.1, show.legend = FALSE) +
  theme_bw() +
  labs(y = "Average Shear Magnitude (m/s)", x = "") +
  scale_fill_manual(values = c("#EE442F", "#63ACBE")) +
  theme(axis.text.x = element_text(size = 15, color = "black", vjust = -1),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14))

data %>% ggplot(aes(x = RI, y = shear, fill = RI)) +
  geom_violin(show.legend = FALSE) +
  #stat_summary(fun.data = "mean_sdl", geom = "pointrange", 
  #             width = 0.05, color = "black",
  #             show.legend = FALSE) +
  geom_boxplot(outlier.size = 1, show.legend = FALSE,
               width = 0.2) +
  theme_bw() +
  labs(y = "Average Shear Magnitude (m/s)", x = "") +
  scale_fill_manual(values = c("#EE442F", "#63ACBE")) +
  theme(axis.text.x = element_text(size = 15, color = "black", vjust = -1),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14))
  

permTS(data %>% filter(RI == "RI") %>% pull(shear),
       data %>% filter(RI == "non-RI") %>% pull(shear),
       alternative = "two.sided",
       method = "exact.mc")

data %>% mutate(RI = factor(RI, levels = c("non-RI", "RI"))) %>%
  ggplot(aes(x = storm, fill = RI)) +
  geom_bar(color = "black", size = 0.25) + scale_fill_manual(values = c("#63ACBE", "#EE442F")) +
  labs(x = "Storm Index", y = "Number of Observations", fill= "") +
  theme_bw(base_size = 13) +
  scale_x_continuous(breaks = seq(0, 130, 20)) +
  scale_y_continuous(breaks = seq(0, 75, 25)) +
  theme(legend.position = c(0.9, 0.925), legend.direction = "horizontal") 

data %>% group_by(storm) %>% 
  summarize(Length = length(storm), Proportion_RI = sum(RI == "RI")/length(RI)) %>%
  ggplot(aes(x = Length, y = Proportion_RI)) +
  geom_point() +
  theme_bw(base_size = 13) +
  labs(x = "Number of Observations in Storm", y = "Proportion RI Observations")


### RI/RW plots
data <- read_csv("data/filtered_storm_list.csv")

ri_plot <- function(id, data){
  
  len = data %>% filter(ID == id) %>% nrow
  
  if (len == 0) {
    print("no data")
    return()
  }
  
  nicole <- data %>%
    filter(ID == id) %>%
    select(WIND, RI) %>%
    bind_cols(data.frame(num = (0:(len-1))*6)) %>%
    filter(num < 39*6)
  
  ggplot(nicole, aes(x = num, y = WIND)) + 
    geom_line() +
    geom_point(
      data = nicole %>% filter(RI == TRUE),
      mapping = aes(x = num, y = WIND),
      shape = 17, color = "red", size = 3
    ) +
    scale_x_continuous(name = "Hours from beginning of observation", 
                       breaks = seq(0, (len-1)*6, 24),
                       minor_breaks = seq(0, (len-1)*6, 6)) +
    scale_y_continuous(name = "Wind Speed (kt)", 
                       breaks = seq(0, 120, 10)) +
    theme_bw(base_size = 14) 
    
}

ri_plot("AL152016", data)

## Storms and observations per month and year
data <- read_csv("data/filtered_storm_list.csv")

# Per year
data %>% group_by(YEAR) %>% summarize(Observations = length(ID)) %>%
  cbind(
    data %>% group_by(YEAR) %>% summarize(Storms = length(unique(ID))) %>% select(Storms)
  ) %>% pivot_longer(c(Observations, Storms), names_to = "type", values_to = "count") %>%
  ggplot(aes(x = YEAR, y = count)) +
  geom_bar(stat = "identity", fill = "#63ACBE") +
  facet_grid(type ~ ., scales = "free_y") +
  theme_bw(base_size = 13) +
  labs(x = "Year", y = "Count")

# Per month
data <- data %>% mutate(month = month(DATETIME))
data$month_str <- factor(mapvalues(data$month, 
                           from = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
                           to = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")),
                         levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

data %>% group_by(month_str, .drop = FALSE) %>% summarize(Observations = length(ID)) %>%
  cbind(
    data %>% group_by(month_str, .drop = FALSE) %>% summarize(Storms = length(unique(ID))) %>% select(Storms)
  ) %>% pivot_longer(c(Observations, Storms), names_to = "type", values_to = "count") %>%
  ggplot(aes(x = month_str, y = count)) +
  geom_bar(stat = "identity", fill = "#63ACBE") +
  facet_grid(type ~ ., scales = "free_y") +
  theme_bw(base_size = 13) +
  labs(x = "Month", y = "Count")

## Integrated circulation EDA
data <- read_csv("data/int_circ_summary.csv") %>%
  mutate(RI = replace(RI, RI == FALSE, "non-RI")) %>%
  mutate(RI = replace(RI, RI == TRUE, "RI")) %>%
  mutate(RI = factor(RI, levels = c("RI", "non-RI"), ordered = TRUE))

data %>% ggplot(aes(x = RI, y = ic, fill = RI)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(outlier.size = 1, show.legend = FALSE,
               width = 0.2) +
  theme_bw() +
  labs(y = "Mean Integrated Circulation", x = "") +
  scale_fill_manual(values = c("#EE442F", "#63ACBE")) +
  theme(axis.text.x = element_text(size = 15, color = "black", vjust = -1),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14))

data %>% ggplot(aes(x = WIND, y = ic)) + geom_point()

# More detailed integrated circulation EDA
data <- read_csv("data/int_circ_summary_detailed.csv") %>%
  mutate(intensity = case_when(
    `25+` ~ '25+',
    `15-25` ~ '[15, 25)',
    `m15-15` ~ '(-15, 15)',
    `m25-m15` ~ '(-25, -15]',
    `25-` ~ '25-'
  )) %>%
  mutate(intensity = factor(intensity, 
                            levels = c('25-', '(-25, -15]', '(-15, 15)', '[15, 25)', '25+'), 
                            ordered = TRUE))

data %>% ggplot(aes(x = intensity, y = ic, fill = intensity)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(outlier.size = 1, show.legend = FALSE,
               width = 0.2) +
  theme_bw() +
  labs(y = "Mean Integrated Circulation", x = "Intensification threshold (kt)") +
  scale_fill_brewer(palette  = "Set3") +
  theme(axis.text.x = element_text(size = 15, color = "black", vjust = -0.25),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14, vjust = -0.75))

data %>%
  filter(intensity == '25+') %>%
  filter(shear < 20) %>%
  ggplot(aes(x = shear, y = ic)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = 'red') +
  labs(y = "Mean Integrated Circulation", x = "Mean Vertical Wind Shear (m/s)", title = "IC vs VWS for rapidly intensifying storms") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14))

data %>%
  filter(intensity == '25+') %>%
  mutate(shear_disc = case_when(
    shear >= 10 ~ '> 10 m/s',
    shear < 10 ~ '< 10 m/s'
  )) %>%
  mutate(shear_disc = factor(shear_disc)) %>%
  ggplot(aes(x = shear_disc, y = ic, fill = shear_disc)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(outlier.size = 1, show.legend = FALSE,
               width = 0.2) +
  theme_bw() +
  labs(y = "Mean Integrated Circulation", x = "Mean Vertical Wind Shear", title = "IC for rapidly intensifying storms in high and low VWS") +
  scale_fill_brewer(palette  = "Set3") +
  theme(axis.text.x = element_text(size = 15, color = "black", vjust = -0.25),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14, vjust = -0.75))

