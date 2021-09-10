library(readr)
library(plyr)
library(tidyr)
library(Hmisc)
library(perm)
library(lubridate)
library(dplyr)
library(ggplot2)
library(magick)
library(scales)
library(naniar)
library(circular)
library(viridis)
library(forcats)
library(metR)
library(np)

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
data <- read_csv("data/int_circ_min_summary_detailed_200hpa.csv") %>%
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


## Filter to RI storms and plot IC maps after ranking by mean shear
ranked_data <- data %>%
  filter(`25+` == TRUE) %>%
  arrange(desc(shear))

n_plot <- 100
figs_file <- "figures/integrated_circulation_vortex_removed/"
for (ii in 1:n_plot) {
  id <- ranked_data[ii,]$ID
  storm_index <- ranked_data[ii,]$storm_index
  storm_name <- ranked_data[ii,]$NAME
  shear_map <- image_read(paste0(figs_file, id, "_", storm_index, ".png"))
  image_annotate(shear_map, paste0(ii, ": ", storm_name, " - ", storm_index), size = 40) %>%
    image_write(path = paste0("figures/ranked_by_shear/", ii, "-", id, "_", storm_index, ".png"), format = "png")
}


# LEVEL SET ANALYSIS (run before anything below) --------------------------
data <- read_csv("data/int_circ_level_set_summary_-0.5_lag-4.csv") %>%
  mutate(intensity = case_when(
    `25+` ~ '25+',
    `15-25` ~ '[15, 25)',
    `m15-15` ~ '(-15, 15)',
    `m25-m15` ~ '(-25, -15]',
    `25-` ~ '25-'
  )) %>%
  mutate(pm20 = case_when(
    `20+` ~ '20+',
    `20-` ~ '20-'
  )) %>%
  mutate(intensity = factor(intensity, 
                            levels = c('25-', '(-25, -15]', '(-15, 15)', '[15, 25)', '25+'), 
                            ordered = TRUE)) %>%
  mutate(velocity_mag = sqrt(velocity_u^2 + velocity_v^2)) %>%
  mutate(basin = substr(ID, 1, 2),
         basin = replace(basin, basin == 'CP', 'EP')) %>% 
  mutate(pm20 = factor(pm20, levels = c('20+', '20-')))


pi_scales <- math_format(.x * pi, format = function(x) x / pi)

df_counts <- data %>%
  filter(ls_size != 0) %>%
  group_by(intensity) %>%
  summarize(n = n(), 
            ls_velocity_angle = median(ls_velocity_angle),
            ls_shear_angle = median(ls_shear_angle),
            ls_north_angle = median(ls_north_angle),
            ls_radius = median(ls_radius),
            ls_size = median(ls_size))

df_counts_pm20_linear_mean <- data %>%
  filter(ls_size != 0, !is.na(pm20)) %>%
  group_by(pm20) %>%
  summarize(n = n(),
            ls_velocity_angle = median(ls_velocity_angle),
            ls_shear_angle = median(ls_shear_angle),
            ls_north_angle = median(ls_north_angle),
            ls_radius = median(ls_radius),
            ls_size = median(ls_size))

# Violin Plots ------------------------------------------------------------

## Velocity angle
data %>% 
  filter(ls_size != 0) %>%
  ggplot(aes(x = intensity, y = ls_velocity_angle, fill = intensity)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(outlier.size = 1, show.legend = FALSE, width = 0.2) +
  stat_summary(fun = circ_mean, geom = "point", color = "black", shape = 17, size = 3, show.legend = FALSE) +
  geom_text(data = df_counts, aes(label=n), nudge_y = 0.25, size = 3) +
  theme_bw() +
  scale_y_continuous(labels = pi_scales, breaks = seq(0, 2*pi, pi / 2)) +
  labs(y = "IC Level Set Angle (TC Velocity Reference)", x = "Intensification threshold (kt)") +
  scale_fill_brewer(palette  = "Set3") +
  theme(axis.text.x = element_text(size = 15, color = "black", vjust = -0.25),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14, vjust = -0.75))

## Shear Angle
data %>% 
  filter(ls_size != 0) %>%
  ggplot(aes(x = intensity, y = ls_shear_angle, fill = intensity)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(outlier.size = 1, show.legend = FALSE, width = 0.2) +
  stat_summary(fun = mean, geom = "point", color = "black", shape = 17, size = 3, show.legend = FALSE) +
  geom_text(data = df_counts, aes(label=n), nudge_y = -0.25, size = 3) +
  theme_bw() +
  scale_y_continuous(labels = pi_scales, breaks = seq(0, 2*pi, pi / 2)) +
  labs(y = "IC Level Set Angle (Shear Reference)", x = "Intensification threshold (kt)") +
  scale_fill_brewer(palette  = "Set3") +
  theme(axis.text.x = element_text(size = 15, color = "black", vjust = -0.25),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14, vjust = -0.75))

## Angle w.r.t. north
data %>% 
  filter(ls_size != 0, ls_radius < 2000) %>%
  ggplot(aes(x = intensity, y = ls_north_angle, fill = intensity)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(outlier.size = 1, show.legend = FALSE, width = 0.2) +
  stat_summary(fun = mean, geom = "point", color = "black", shape = 17, size = 3, show.legend = FALSE) +
  geom_text(data = df_counts, aes(label=n), nudge_y = 0.25, size = 3) +
  theme_bw() +
  scale_y_continuous(labels = pi_scales, breaks = seq(0, 2*pi, pi / 2)) +
  labs(y = "IC Level Set Angle (North Reference)", x = "Intensification threshold (kt)") +
  scale_fill_brewer(palette  = "Set3") +
  theme(axis.text.x = element_text(size = 15, color = "black", vjust = -0.25),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14, vjust = -0.75))

## Radius
data %>% 
  filter(ls_size != 0) %>%
  ggplot(aes(x = intensity, y = ls_radius, fill = intensity)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(outlier.size = 1, show.legend = FALSE, width = 0.2) +
  stat_summary(fun = mean, geom = "point", color = "black", shape = 17, size = 3, show.legend = FALSE) +
  geom_text(data = df_counts, aes(label=n), nudge_y = 50, size = 3) +
  theme_bw() +
  labs(y = "IC Level Set Radius from Storm Center (km)", x = "Intensification threshold (kt)") +
  scale_fill_brewer(palette  = "Set3") +
  theme(axis.text.x = element_text(size = 15, color = "black", vjust = -0.25),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14, vjust = -0.75))

## Size
data %>% 
  filter(ls_size != 0) %>%
  ggplot(aes(x = intensity, y = ls_size, fill = intensity)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(outlier.size = 1, show.legend = FALSE, width = 0.2) +
  stat_summary(fun = mean, geom = "point", color = "black", shape = 17, size = 3, show.legend = FALSE) +
  geom_text(data = df_counts, aes(label=n), nudge_y = 75, size = 3) +
  theme_bw() +
  labs(y = "IC Level Set Size (# grid points)", x = "Intensification threshold (kt)") +
  scale_fill_brewer(palette  = "Set3") +
  theme(axis.text.x = element_text(size = 15, color = "black", vjust = -0.25),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14, vjust = -0.75))

size_data_na <- data %>% 
  mutate(ls_size = ifelse(ls_size == 0, NA, ls_size)) %>%
  bind_shadow() %>%
  mutate(ls_size_NA = relevel((data %>% 
                              mutate(ls_size = ifelse(ls_size == 0, NA, ls_size)) %>%
                              bind_shadow())$ls_size_NA, "NA"
                           )
         ) %>%
  group_by(intensity) %>%
  summarize(missing_prop = round(sum(ls_size_NA == "NA")/n(), 3),
            ls_size = -50)

size_data <- data %>%
  filter(ls_size != 0)

ggplot(size_data, aes(x = intensity, y = ls_size, fill = intensity)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(outlier.size = 1, show.legend = FALSE, width = 0.2) +
  stat_summary(fun = mean, geom = "point", color = "black", shape = 17, size = 3, show.legend = FALSE) +
  geom_text(data = df_counts, aes(label=n), nudge_y = 75, size = 3) +
  geom_text(data = size_data_na, aes(label = missing_prop), size = 4, color = "red") +
  #geom_jitter(data = size_data_na %>% filter(ls_size_NA == "NA"), aes(x = intensity, y = -50), alpha = 0.5, show.legend = FALSE) +
  theme_bw() +
  labs(y = "IC Level Set Size (# grid points)", x = "Intensification threshold (kt)") +
  scale_fill_brewer(palette  = "Set3") +
  theme(axis.text.x = element_text(size = 15, color = "black", vjust = -0.25),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14, vjust = -0.75))  
  
data %>% 
  bind_shadow() %>%
  mutate(Vaccinated_NA = relevel((data %>% bind_shadow())$Vaccinated_NA, "NA")) %>%
  ggplot(aes(x = Hesitant, fill = Vaccinated_NA)) +
  geom_density(alpha = 0.35, bw = 0.01) +
  labs(x = "Vaccine Hesitancy", y = "Density") +
  scale_fill_discrete(name = "", labels = c("Missing", "Not Missing")) +
  theme_bw() +
  theme(text = element_text(size = 22),
        legend.key.height = unit(1, "cm"),
        legend.position = "top")



# Circular Plots/2 Categories ---------------------------------------------
#### setup ####
circ_mean <- function(x){
  x %>%
    circular::circular(units = "radians") %>%
    circular::mean.circular(na.rm = TRUE)
}

df_counts_pm20 <- data %>%
  filter(ls_size != 0, !is.na(pm20)) %>%
  group_by(pm20) %>%
  summarize(n = n(),
            ls_velocity_angle = circ_mean(ls_velocity_angle),
            ls_shear_angle = circ_mean(ls_shear_angle),
            ls_north_angle = circ_mean(ls_north_angle),
            ls_radius = median(ls_radius),
            ls_size = median(ls_size)) %>%
  mutate(
    ls_velocity_angle = case_when(
      ls_velocity_angle < 0 ~ ls_velocity_angle + 2*pi,
      TRUE ~ ls_velocity_angle
    ),
    ls_shear_angle = case_when(
      ls_shear_angle < 0 ~ ls_shear_angle + 2*pi,
      TRUE ~ ls_shear_angle
    ),
    ls_north_angle = case_when(
      ls_north_angle < 0 ~ ls_north_angle + 2*pi,
      TRUE ~ ls_north_angle
    )
  )

df_counts_pm20_vel <- data %>%
  filter(ls_size != 0, !is.na(pm20)) %>%
  filter(velocity_mag * 100/(6*1.852) > 5) %>% # Filter to storms that are approx > 5 knot velocity
  group_by(pm20) %>%
  summarize(n = n(),
            ls_velocity_angle = circ_mean(ls_velocity_angle),
            ls_shear_angle = circ_mean(ls_shear_angle),
            ls_north_angle = circ_mean(ls_north_angle),
            ls_radius = median(ls_radius),
            ls_size = median(ls_size)) %>%
  mutate(
    ls_velocity_angle = case_when(
      ls_velocity_angle < 0 ~ ls_velocity_angle + 2*pi,
      TRUE ~ ls_velocity_angle
    ),
    ls_shear_angle = case_when(
      ls_shear_angle < 0 ~ ls_shear_angle + 2*pi,
      TRUE ~ ls_shear_angle
    ),
    ls_north_angle = case_when(
      ls_north_angle < 0 ~ ls_north_angle + 2*pi,
      TRUE ~ ls_north_angle
    )
  )

df_counts_pm20_shear <- data %>%
  filter(ls_size != 0, !is.na(pm20)) %>%
  filter(shear_mag > 5) %>% # Filter to storms that are approx > 5 m/s shear
  group_by(pm20) %>%
  summarize(n = n(),
            ls_velocity_angle = circ_mean(ls_velocity_angle),
            ls_shear_angle = circ_mean(ls_shear_angle),
            ls_north_angle = circ_mean(ls_north_angle),
            ls_radius = median(ls_radius),
            ls_size = median(ls_size)) %>%
  mutate(
    ls_velocity_angle = case_when(
      ls_velocity_angle < 0 ~ ls_velocity_angle + 2*pi,
      TRUE ~ ls_velocity_angle
    ),
    ls_shear_angle = case_when(
      ls_shear_angle < 0 ~ ls_shear_angle + 2*pi,
      TRUE ~ ls_shear_angle
    ),
    ls_north_angle = case_when(
      ls_north_angle < 0 ~ ls_north_angle + 2*pi,
      TRUE ~ ls_north_angle
    )
  )

#### velocity ####
data %>% 
  filter(ls_size != 0, !is.na(pm20), velocity_mag * 100/(6*1.852) > 5, velocity_mag < 200) %>%
  ggplot(aes(x = ls_velocity_angle, fill = pm20)) + 
  geom_histogram(breaks = seq(0, 2*pi, by = pi/16), position = "stack", color = "black") + 
  geom_vline(data = df_counts_pm20_vel, aes(xintercept = ls_velocity_angle), size = 1.5) +
  coord_polar(direction = -1) +
  scale_x_continuous("", limits = c(0, 2*pi), 
                     breaks = seq(0, 2*pi, by = pi/2), labels = pi_scales)+
  scale_y_continuous("",limits=c(-20,60)) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~ pm20) +
  ggtitle("TC Velocity Referenced IC Level Set Angle Distribution")

# Take a look at the distribution without filtering out low velocity observations
data %>% 
  filter(ls_size != 0, !is.na(pm20)) %>%
  ggplot(aes(x = ls_velocity_angle, fill = pm20)) + 
  geom_histogram(breaks = seq(0, 2*pi, by = pi/16), position = "stack", color = "black") + 
  geom_vline(data = df_counts_pm20, aes(xintercept = ls_velocity_angle), size = 1.5) +
  coord_polar(direction = -1) +
  scale_x_continuous("", limits = c(0, 2*pi), 
                     breaks = seq(0, 2*pi, by = pi/2), labels = pi_scales)+
  scale_y_continuous("",limits=c(-20,60)) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~ pm20) +
  ggtitle("TC Velocity Referenced IC Level Set Angle Distribution")

#### shear ####
data %>% 
  filter(ls_size != 0, !is.na(pm20), shear_mag > 5, velocity_mag < 200) %>%
  ggplot(aes(x = ls_shear_angle, fill = pm20)) + 
  geom_histogram(breaks = seq(0, 2*pi, by = pi/16), position = "stack", color = "black") + 
  geom_vline(data = df_counts_pm20_shear, aes(xintercept = ls_shear_angle), size = 1.5) +
  coord_polar(direction = -1) +
  scale_x_continuous("", limits = c(0, 2*pi), 
                     breaks = seq(0, 2*pi, by = pi/2), labels = pi_scales)+
  scale_y_continuous("",limits=c(-20,60)) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~ pm20) +
  ggtitle("TC Shear Referenced IC Level Set Angle Distribution")

# Take a look at the distribution without filtering out low shear observations
data %>% 
  filter(ls_size != 0, !is.na(pm20)) %>%
  ggplot(aes(x = ls_shear_angle, fill = pm20)) + 
  geom_histogram(breaks = seq(0, 2*pi, by = pi/16), position = "stack", color = "black") + 
  geom_vline(data = df_counts_pm20, aes(xintercept = ls_shear_angle), size = 1.5) +
  coord_polar(direction = -1) +
  scale_x_continuous("", limits = c(0, 2*pi), 
                     breaks = seq(0, 2*pi, by = pi/2), labels = pi_scales)+
  scale_y_continuous("",limits=c(-20,60)) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~ pm20) +
  ggtitle("TC Shear Referenced IC Level Set Angle Distribution")

#### north/south ####
data %>% 
  filter(ls_size != 0, !is.na(pm20), velocity_mag < 200) %>%
  ggplot(aes(x = ls_north_angle, fill = pm20)) + 
  geom_histogram(breaks = seq(0, 2*pi, by = pi/16), position = "stack", color = "black") + 
  geom_vline(data = df_counts_pm20, aes(xintercept = ls_north_angle), size = 1.5) +
  coord_polar(direction = -1) +
  scale_x_continuous("", limits = c(0, 2*pi), 
                     breaks = seq(0, 2*pi, by = pi/2), labels = pi_scales)+
  scale_y_continuous("",limits=c(-20,60)) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~ pm20) +
  ggtitle("North Referenced IC Level Set Angle Distribution")



# Scatter Plots -----------------------------------------------------------
#### velocity ####
data %>%
  filter(ls_size != 0, !is.na(pm20)) %>%
  filter(velocity_mag * 100/(6*1.852) > 5, velocity_mag < 200) %>%
  mutate(ls_vel_x = -1 * ls_radius * sin(ls_velocity_angle), ls_vel_y = ls_radius * cos(ls_velocity_angle)) %>%
  ggplot(aes(x = ls_vel_x, y = ls_vel_y, color = pm20)) +
  geom_point(alpha = 0.6) +
  labs(x = 'Component Perpendicular to Velocity (km)',
       y = 'Component Parallel to Velocity (km)',
       color = "Intensification\nCategory",
       title = "TC Velocity Referenced IC Level Set Locations") +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "none")

#### shear ####
data %>%
  filter(ls_size != 0, !is.na(pm20)) %>%
  mutate(ls_shear_x = -1 * ls_radius * sin(ls_shear_angle), ls_shear_y = ls_radius * cos(ls_shear_angle)) %>%
  ggplot(aes(x = ls_shear_x, y = ls_shear_y, color = pm20)) +
  geom_point(alpha = 0.6) +
  labs(x = 'Component Perpendicular to Shear (km)',
       y = 'Component Parallel to Shear (km)',
       color = "Intensification\nCategory",
       title = "TC Shear Referenced IC Level Set Locations") +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "none")

#### north/south ####
data %>%
  filter(ls_size != 0, !is.na(pm20)) %>%
  filter(velocity_mag < 200) %>%
  filter(ls_u < 6.5, ls_u > -6.5) %>%
  ggplot(aes(x = ls_u, y = ls_v, color = pm20)) +
  geom_point(alpha = 0.6) +
  labs(x = 'Zonal (degrees)',
       y = 'Meridional (degrees)',
       color = "Intensification\nCategory",
       title = "North Referenced IC Level Set Locations") +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "none")


## Proportion hex plots
data %>%
  filter(ls_size != 0, !is.na(pm20)) %>%
  mutate(pm20_ind = ifelse(pm20 == '20+', 1, 0)) %>%
  ggplot(aes(x = ls_u, y = ls_v, z = pm20_ind)) + 
  stat_summary_hex(fun = function(z){(sum(z)/length(z))*(727/762)}, bins = 17) + 
  scale_fill_viridis(name = "Proportion\n20+") +
  theme_bw() +
  theme() +
  coord_fixed()

data %>%
  filter(ls_size != 0, !is.na(pm20)) %>%
  mutate(pm20_ind = ifelse(pm20 == '20+', 1, 0)) %>%
  mutate(ls_vel_x = -1 * ls_radius * sin(ls_velocity_angle), ls_vel_y = ls_radius * cos(ls_velocity_angle)) %>%
  ggplot(aes(x = ls_vel_x, y = ls_vel_y, z = pm20_ind)) + 
  stat_summary_hex(fun = function(z){ifelse(length(z) > 2, (sum(z)/length(z))*(727/762), NA)}, bins = 15) + 
  scale_fill_viridis(name = "Proportion\n20+") +
  theme_bw() +
  theme() +
  coord_fixed()



# Density Estimation ------------------------------------------------------
#### velocity ####
de_vel_data <- data %>%
  filter(ls_size != 0, !is.na(pm20)) %>%
  filter(velocity_mag * 100/(6*1.852) > 5, velocity_mag < 200) %>%
  mutate(ls_vel_x = -1 * ls_radius * sin(ls_velocity_angle), ls_vel_y = ls_radius * cos(ls_velocity_angle))

grid_pts <- 50
ls_vel_x.seq <- seq(min(de_vel_data$ls_vel_x), max(de_vel_data$ls_vel_x), length.out = grid_pts)
ls_vel_y.seq <- seq(min(de_vel_data$ls_vel_y), max(de_vel_data$ls_vel_y), length.out = grid_pts)

de_vel_eval <- expand.grid(ls_vel_x = ls_vel_x.seq, ls_vel_y = ls_vel_y.seq)

bw_create <- function(data, bw_size){
  npudensbw(formula = ~ ls_vel_x + ls_vel_y, 
            data = data,
            bws = c(bw_size, bw_size),
            bandwidth.compute = FALSE)
}

bw_size = 75
de_vel_bw_20p <- filter(de_vel_data, pm20 == '20+') %>% bw_create(bw_size)
de_vel_bw_20m <- filter(de_vel_data, pm20 == '20-') %>% bw_create(bw_size)
de_vel_bw_20p_AL <- filter(de_vel_data, pm20 == '20+', basin == "AL") %>% bw_create(bw_size)
de_vel_bw_20p_EP <- filter(de_vel_data, pm20 == '20+', basin == "EP") %>% bw_create(bw_size)
de_vel_bw_20m_AL <- filter(de_vel_data, pm20 == '20-', basin == "AL") %>% bw_create(bw_size)
de_vel_bw_20m_EP <- filter(de_vel_data, pm20 == '20-', basin == "EP") %>% bw_create(bw_size)

de_vel_20p <- npudens(de_vel_bw_20p, newdata = de_vel_eval)
de_vel_20m <- npudens(de_vel_bw_20m, newdata = de_vel_eval)
de_vel_20p_AL <- npudens(de_vel_bw_20p_AL, newdata = de_vel_eval)
de_vel_20p_EP <- npudens(de_vel_bw_20p_EP, newdata = de_vel_eval)
de_vel_20m_AL <- npudens(de_vel_bw_20m_AL, newdata = de_vel_eval)
de_vel_20m_EP <- npudens(de_vel_bw_20m_EP, newdata = de_vel_eval)

de_vel_toplot_init <- data.frame(
  x = de_vel_20p_AL$eval$ls_vel_x,
  y = de_vel_20p_AL$eval$ls_vel_y,
  density_20p = de_vel_20p$dens,
  density_20m = de_vel_20m$dens
)
de_vel_toplot <- pivot_longer(de_vel_toplot_init,
                              cols = starts_with('density'),
                              names_to = "intensity",
                              names_prefix = "density_",
                              values_to = "density") %>%
  mutate(intensity = factor(intensity, levels = c('20p', '20m'), labels = c('20+ (Intensifying)', '20- (Weakening)')))

de_vel_toplot_init_bybasin <- data.frame(
  x = de_vel_20p_AL$eval$ls_vel_x,
  y = de_vel_20p_AL$eval$ls_vel_y,
  density_20p_AL = de_vel_20p_AL$dens,
  density_20p_EP = de_vel_20p_EP$dens,
  density_20m_AL = de_vel_20m_AL$dens,
  density_20m_EP = de_vel_20m_EP$dens
)
de_vel_toplot_bybasin <- pivot_longer(de_vel_toplot_init_bybasin,
                                      cols = starts_with('density'),
                                      names_to = "intensity_basin",
                                      names_prefix = "density_",
                                      values_to = "density") %>%
  mutate(intensity = substr(intensity_basin, 1, 3),
         basin = substr(intensity_basin, 5, 6)) %>%
  mutate(intensity = factor(intensity, levels = c('20p', '20m'), labels = c('20+ (Intensifying)', '20- (Weakening)')),
         basin = factor(basin, levels = c('AL', 'EP'), labels = c('NAL', 'ENP')))

de_vel_counts <- de_vel_data %>%
  group_by(pm20) %>%
  summarize(n = n(), mean_size = round(mean(ls_size), 1), mean_avg_ic = mean(ls_avg_ic)) %>%
  mutate(x = 0, y = 0)%>%
  rename(intensity = pm20) %>%
  mutate(intensity = factor(intensity, levels = c('20+', '20-'), labels = c('20+ (Intensifying)', '20- (Weakening)')))

de_vel_counts_bybasin <- de_vel_data %>%
  group_by(pm20, basin) %>%
  summarize(n = n(), mean_size = round(mean(ls_size), 1), mean_avg_ic = mean(ls_avg_ic)) %>%
  mutate(x = 0, y = 0) %>%
  rename(intensity = pm20) %>%
  mutate(intensity = factor(intensity, levels = c('20+', '20-'), labels = c('20+ (Intensifying)', '20- (Weakening)')),
         basin = factor(basin, levels = c('AL', 'EP'), labels = c('NAL', 'ENP')))

# by intensity
de_vel_toplot %>%
  ggplot(aes(x = x, y = y)) + 
  geom_raster(aes(fill = density)) + 
  geom_hline(yintercept = 0, color = 'red') +
  geom_vline(xintercept = 0, color = 'red') +
  geom_text(data = de_vel_counts, aes(label=paste0('Samples: ', n)), nudge_y = -576, nudge_x = -370, size = 3, col = 'white') +
  geom_text(data = de_vel_counts, aes(label=paste0('Avg. LS size: ', mean_size)), nudge_y = -475, nudge_x = -275, size = 3, col = 'white') +
  labs(x = "Component Perpendicular to Velocity (km)", 
       y = "Component Parallel to Velocity (km)", 
       fill = "Density") +
  scale_fill_viridis(limits = c(0, max(de_vel_toplot$density))) +
  coord_equal() +
  theme_bw() +
  facet_grid(~ intensity)

# by intensity + basin
de_vel_toplot_bybasin %>%
  ggplot(aes(x = x, y = y)) + 
  geom_raster(aes(fill = density)) + 
  geom_hline(yintercept = 0, color = 'red') +
  geom_vline(xintercept = 0, color = 'red') +
  geom_text(data = de_vel_counts_bybasin, aes(label=paste0('Samples: ', n)), nudge_y = -576, nudge_x = -370, size = 3, col = 'white') +
  geom_text(data = de_vel_counts_bybasin, aes(label=paste0('Avg. LS size: ', mean_size)), nudge_y = -475, nudge_x = -275, size = 3, col = 'white') +
  labs(x = "Component Perpendicular to Velocity (km)", 
       y = "Component Parallel to Velocity (km)", 
       fill = "Density") +
  scale_fill_viridis(limits = c(0, max(de_vel_toplot_bybasin$density))) +
  coord_equal() +
  theme_bw() +
  facet_grid(basin ~ intensity)

#### shear ####
de_she_data <- data %>%
  filter(ls_size != 0, !is.na(pm20)) %>%
  filter(shear_mag > 5, velocity_mag < 200) %>%
  mutate(ls_she_x = -1 * ls_radius * sin(ls_shear_angle), ls_she_y = ls_radius * cos(ls_shear_angle))

grid_pts <- 50
ls_she_x.seq <- seq(min(de_she_data$ls_she_x), max(de_she_data$ls_she_x), length.out = grid_pts)
ls_she_y.seq <- seq(min(de_she_data$ls_she_y), max(de_she_data$ls_she_y), length.out = grid_pts)

de_she_eval <- expand.grid(ls_she_x = ls_she_x.seq, ls_she_y = ls_she_y.seq)

bw_create <- function(data, bw_size){
  npudensbw(formula = ~ ls_she_x + ls_she_y, 
            data = data,
            bws = c(bw_size, bw_size),
            bandwidth.compute = FALSE)
}

bw_size = 75
de_she_bw_20p <- filter(de_she_data, pm20 == '20+') %>% bw_create(bw_size)
de_she_bw_20m <- filter(de_she_data, pm20 == '20-') %>% bw_create(bw_size)
de_she_bw_20p_AL <- filter(de_she_data, pm20 == '20+', basin == "AL") %>% bw_create(bw_size)
de_she_bw_20p_EP <- filter(de_she_data, pm20 == '20+', basin == "EP") %>% bw_create(bw_size)
de_she_bw_20m_AL <- filter(de_she_data, pm20 == '20-', basin == "AL") %>% bw_create(bw_size)
de_she_bw_20m_EP <- filter(de_she_data, pm20 == '20-', basin == "EP") %>% bw_create(bw_size)

de_she_20p <- npudens(de_she_bw_20p, newdata = de_she_eval)
de_she_20m <- npudens(de_she_bw_20m, newdata = de_she_eval)
de_she_20p_AL <- npudens(de_she_bw_20p_AL, newdata = de_she_eval)
de_she_20p_EP <- npudens(de_she_bw_20p_EP, newdata = de_she_eval)
de_she_20m_AL <- npudens(de_she_bw_20m_AL, newdata = de_she_eval)
de_she_20m_EP <- npudens(de_she_bw_20m_EP, newdata = de_she_eval)

de_she_toplot_init <- data.frame(
  x = de_she_20p_AL$eval$ls_she_x,
  y = de_she_20p_AL$eval$ls_she_y,
  density_20p = de_she_20p$dens,
  density_20m = de_she_20m$dens
)
de_she_toplot <- pivot_longer(de_she_toplot_init,
                              cols = starts_with('density'),
                              names_to = "intensity",
                              names_prefix = "density_",
                              values_to = "density") %>%
  mutate(intensity = factor(intensity, levels = c('20p', '20m'), labels = c('20+ (Intensifying)', '20- (Weakening)')))

de_she_toplot_init_bybasin <- data.frame(
  x = de_she_20p_AL$eval$ls_she_x,
  y = de_she_20p_AL$eval$ls_she_y,
  density_20p_AL = de_she_20p_AL$dens,
  density_20p_EP = de_she_20p_EP$dens,
  density_20m_AL = de_she_20m_AL$dens,
  density_20m_EP = de_she_20m_EP$dens
)
de_she_toplot_bybasin <- pivot_longer(de_she_toplot_init_bybasin,
                                     cols = starts_with('density'),
                                     names_to = "intensity_basin",
                                     names_prefix = "density_",
                                     values_to = "density") %>%
  mutate(intensity = substr(intensity_basin, 1, 3),
         basin = substr(intensity_basin, 5, 6)) %>%
  mutate(intensity = factor(intensity, levels = c('20p', '20m'), labels = c('20+ (Intensifying)', '20- (Weakening)')),
         basin = factor(basin, levels = c('AL', 'EP'), labels = c('NAL', 'ENP')))

de_she_counts <- de_she_data %>%
  group_by(pm20) %>%
  summarize(n = n(), mean_size = round(mean(ls_size), 1), mean_avg_ic = mean(ls_avg_ic)) %>%
  mutate(x = 0, y = 0)%>%
  rename(intensity = pm20) %>%
  mutate(intensity = factor(intensity, levels = c('20+', '20-'), labels = c('20+ (Intensifying)', '20- (Weakening)')))

de_she_counts_bybasin <- de_she_data %>%
  group_by(pm20, basin) %>%
  summarize(n = n(), mean_size = round(mean(ls_size), 1), mean_avg_ic = mean(ls_avg_ic)) %>%
  mutate(x = 0, y = 0) %>%
  rename(intensity = pm20) %>%
  mutate(intensity = factor(intensity, levels = c('20+', '20-'), labels = c('20+ (Intensifying)', '20- (Weakening)')),
         basin = factor(basin, levels = c('AL', 'EP'), labels = c('NAL', 'ENP')))

# by intensity
de_she_toplot %>%
  ggplot(aes(x = x, y = y)) + 
  geom_raster(aes(fill = density)) + 
  geom_hline(yintercept = 0, color = 'red') +
  geom_vline(xintercept = 0, color = 'red') +
  geom_text(data = de_she_counts, aes(label=paste0('Samples: ', n)), nudge_y = -576, nudge_x = -370, size = 3, col = 'white') +
  geom_text(data = de_she_counts, aes(label=paste0('Avg. LS size: ', mean_size)), nudge_y = -475, nudge_x = -275, size = 3, col = 'white') +
  labs(x = "Component Perpendicular to Shear (km)", 
       y = "Component Parallel to Shear (km)", 
       fill = "Density") +
  scale_fill_viridis(limits = c(0, max(de_she_toplot$density))) +
  coord_equal() +
  theme_bw() +
  facet_grid(~ intensity)

# by intensity + basin
de_she_toplot_bybasin %>%
  ggplot(aes(x = x, y = y)) + 
  geom_raster(aes(fill = density)) + 
  geom_hline(yintercept = 0, color = 'red') +
  geom_vline(xintercept = 0, color = 'red') +
  geom_text(data = de_she_counts_bybasin, aes(label=paste0('Samples: ', n)), nudge_y = -576, nudge_x = -370, size = 3, col = 'white') +
  geom_text(data = de_she_counts_bybasin, aes(label=paste0('Avg. LS size: ', mean_size)), nudge_y = -475, nudge_x = -275, size = 3, col = 'white') +
  labs(x = "Component Perpendicular to Shear (km)", 
       y = "Component Parallel to Shear (km)", 
       fill = "Density") +
  scale_fill_viridis(limits = c(0, max(de_she_toplot_bybasin$density))) +
  coord_equal() +
  theme_bw() +
  facet_grid(basin ~ intensity)


#### north/south ####
de_ns_data <- data %>%
  filter(ls_size != 0, !is.na(pm20)) %>%
  filter(velocity_mag < 200) %>%
  mutate(ls_ns_x = ls_u, ls_ns_y = ls_v)

grid_pts <- 50
ls_ns_x.seq <- seq(min(de_ns_data$ls_ns_x), max(de_ns_data$ls_ns_x), length.out = grid_pts)
ls_ns_y.seq <- seq(min(de_ns_data$ls_ns_y), max(de_ns_data$ls_ns_y), length.out = grid_pts)

de_ns_eval <- expand.grid(ls_ns_x = ls_ns_x.seq, ls_ns_y = ls_ns_y.seq)

bw_create <- function(data, bw_size){
  npudensbw(formula = ~ ls_ns_x + ls_ns_y, 
            data = data,
            bws = c(bw_size, bw_size),
            bandwidth.compute = FALSE)
}

bw_size = 0.75
de_ns_bw_20p <- filter(de_ns_data, pm20 == '20+') %>% bw_create(bw_size)
de_ns_bw_20m <- filter(de_ns_data, pm20 == '20-') %>% bw_create(bw_size)
de_ns_bw_20p_AL <- filter(de_ns_data, pm20 == '20+', basin == "AL") %>% bw_create(bw_size)
de_ns_bw_20p_EP <- filter(de_ns_data, pm20 == '20+', basin == "EP") %>% bw_create(bw_size)
de_ns_bw_20m_AL <- filter(de_ns_data, pm20 == '20-', basin == "AL") %>% bw_create(bw_size)
de_ns_bw_20m_EP <- filter(de_ns_data, pm20 == '20-', basin == "EP") %>% bw_create(bw_size)

de_ns_20p <- npudens(de_ns_bw_20p, newdata = de_ns_eval)
de_ns_20m <- npudens(de_ns_bw_20m, newdata = de_ns_eval)
de_ns_20p_AL <- npudens(de_ns_bw_20p_AL, newdata = de_ns_eval)
de_ns_20p_EP <- npudens(de_ns_bw_20p_EP, newdata = de_ns_eval)
de_ns_20m_AL <- npudens(de_ns_bw_20m_AL, newdata = de_ns_eval)
de_ns_20m_EP <- npudens(de_ns_bw_20m_EP, newdata = de_ns_eval)

de_ns_toplot_init <- data.frame(
  x = de_ns_20p_AL$eval$ls_ns_x,
  y = de_ns_20p_AL$eval$ls_ns_y,
  density_20p = de_ns_20p$dens,
  density_20m = de_ns_20m$dens
)
de_ns_toplot <- pivot_longer(de_ns_toplot_init,
                             cols = starts_with('density'),
                             names_to = "intensity",
                             names_prefix = "density_",
                             values_to = "density") %>%
  mutate(intensity = factor(intensity, levels = c('20p', '20m'), labels = c('20+ (Intensifying)', '20- (Weakening)')))


de_ns_toplot_init_bybasin <- data.frame(
  x = de_ns_20p_AL$eval$ls_ns_x,
  y = de_ns_20p_AL$eval$ls_ns_y,
  density_20p_AL = de_ns_20p_AL$dens,
  density_20p_EP = de_ns_20p_EP$dens,
  density_20m_AL = de_ns_20m_AL$dens,
  density_20m_EP = de_ns_20m_EP$dens
)
de_ns_toplot_bybasin <- pivot_longer(de_ns_toplot_init_bybasin,
                              cols = starts_with('density'),
                              names_to = "intensity_basin",
                              names_prefix = "density_",
                              values_to = "density") %>%
  mutate(intensity = substr(intensity_basin, 1, 3),
         basin = substr(intensity_basin, 5, 6)) %>%
  mutate(intensity = factor(intensity, levels = c('20p', '20m'), labels = c('20+ (Intensifying)', '20- (Weakening)')),
         basin = factor(basin, levels = c('AL', 'EP'), labels = c('NAL', 'ENP')))

de_ns_counts <- de_ns_data %>%
  group_by(pm20) %>%
  summarize(n = n(), mean_size = round(mean(ls_size), 1), mean_avg_ic = mean(ls_avg_ic)) %>%
  mutate(x = 0, y = 0)%>%
  rename(intensity = pm20) %>%
  mutate(intensity = factor(intensity, levels = c('20+', '20-'), labels = c('20+ (Intensifying)', '20- (Weakening)')))

de_ns_counts_bybasin <- de_ns_data %>%
  group_by(pm20, basin) %>%
  summarize(n = n(), mean_size = round(mean(ls_size), 1), mean_avg_ic = mean(ls_avg_ic)) %>%
  mutate(x = 0, y = 0) %>%
  rename(intensity = pm20) %>%
  mutate(intensity = factor(intensity, levels = c('20+', '20-'), labels = c('20+ (Intensifying)', '20- (Weakening)')),
         basin = factor(basin, levels = c('AL', 'EP'), labels = c('NAL', 'ENP')))

# by intensity
de_ns_toplot %>%
  ggplot(aes(x = x, y = y)) + 
  geom_raster(aes(fill = density)) + 
  geom_hline(yintercept = 0, color = 'red') +
  geom_vline(xintercept = 0, color = 'red') +
  geom_text(data = de_ns_counts, aes(label=paste0('Samples: ', n)), nudge_y = -5, nudge_x = -3.5, size = 3, col = 'white') +
  geom_text(data = de_ns_counts, aes(label=paste0('Avg. LS size: ', mean_size)), nudge_y = -4, nudge_x = -2.5, size = 3, col = 'white') +
  labs(x = "Zonal (degrees)", 
       y = "Meridional (degrees)", 
       fill = "Density") +
  scale_fill_viridis(limits = c(0, max(de_ns_toplot$density))) +
  theme_bw() +
  facet_grid(~ intensity)

# by intensity + basin
de_ns_toplot_bybasin %>%
  ggplot(aes(x = x, y = y)) + 
  geom_raster(aes(fill = density)) + 
  geom_hline(yintercept = 0, color = 'red') +
  geom_vline(xintercept = 0, color = 'red') +
  geom_text(data = de_ns_counts_bybasin, aes(label=paste0('Samples: ', n)), nudge_y = -5, nudge_x = -3.5, size = 3, col = 'white') +
  geom_text(data = de_ns_counts_bybasin, aes(label=paste0('Avg. LS size: ', mean_size)), nudge_y = -4, nudge_x = -2.5, size = 3, col = 'white') +
  labs(x = "Zonal (degrees)", 
       y = "Meridional (degrees)", 
       fill = "Density") +
  scale_fill_viridis(limits = c(0, max(de_ns_toplot_bybasin$density))) +
  theme_bw() +
  facet_grid(basin ~ intensity)




# Density Estimation with Lag ---------------------------------------------
#### velocity ####
#### shear ####
#### north/south ####

