# Load necessary libraries
# install.packages(c("survival", "ggplot2", "survminer", "dplyr", "svglite"))
library(survival)
library(ggplot2)
library(survminer)
library(dplyr)
library(svglite)

# --- 1. Plot for Proportional Hazards (PH) Assumption Met ---

set.seed(123)
n <- 200
# Simulate event times with a constant hazard ratio (HR = 0.5)
event_time_ph_met <- rexp(n, rate = 0.05 * rep(c(1, 0.5), each = n/2))
# Simulate censoring times
censoring_time_ph_met <- rexp(n, rate = 0.02)

sim_data_ph_met <- data.frame(
  time = pmin(event_time_ph_met, censoring_time_ph_met),
  status = as.numeric(event_time_ph_met <= censoring_time_ph_met),
  group = factor(rep(c("Control", "Treatment"), each = n/2))
)

fit_met <- survfit(Surv(time, status) ~ group, data = sim_data_ph_met)

p_met <- ggsurvplot(
  fit_met,
  data = sim_data_ph_met,
  palette = c("#E7B800", "#2E9FDF"),
  title = "PH Assumption Met",
  subtitle = "Hazard ratio is constant over time.",
  legend.title = "Group",
  legend.labs = c("Control", "Treatment"),
  xlab = "Time",
  ylab = "Survival Probability",
  censor.shape="|", censor.size = 4,
  ggtheme = theme_classic()
)

ggsave("images/ph_assumption_met.svg", plot = p_met$plot, width = 6, height = 4)
ggsave("images/ph_assumption_met.png", plot = p_met$plot, width = 6, height = 4)


# --- 2. Plot for Proportional Hazards (PH) Assumption Violated ---

set.seed(456)
n <- 200
# Simulate event times with crossing survival curves (violating PH)
time_control <- rexp(n/2, rate = 0.04)
time_treatment <- c(rexp(n/4, rate = 0.1), rexp(n/4, rate = 0.02)) # Delayed effect
event_time_ph_violated <- c(time_control, time_treatment)
# Simulate censoring times
censoring_time_ph_violated <- rexp(n, rate = 0.015)

sim_data_ph_violated <- data.frame(
  time = pmin(event_time_ph_violated, censoring_time_ph_violated),
  status = as.numeric(event_time_ph_violated <= censoring_time_ph_violated),
  group = factor(rep(c("Control", "Treatment"), each = n/2))
)

fit_violated <- survfit(Surv(time, status) ~ group, data = sim_data_ph_violated)

p_violated <- ggsurvplot(
  fit_violated,
  data = sim_data_ph_violated,
  palette = c("#E7B800", "#2E9FDF"),
  title = "PH Assumption Violated",
  subtitle = "Survival curves cross, hazard ratio is not constant.",
  legend.title = "Group",
  legend.labs = c("Control", "Treatment"),
  xlab = "Time",
  ylab = "Survival Probability",
  censor.shape="|", censor.size = 4,
  ggtheme = theme_classic()
)

ggsave("images/ph_assumption_violated.svg", plot = p_violated$plot, width = 6, height = 4)
ggsave("images/ph_assumption_violated.png", plot = p_violated$plot, width = 6, height = 4)


# --- 3. RMST Plot for PH Assumption Met Case ---

L <- 70 # Truncation time

# Use the data from the PH Met simulation
fit_rmst_met <- survfit(Surv(time, status) ~ group, data = sim_data_ph_met)
p_rmst_base_met <- ggsurvplot(fit_rmst_met, data = sim_data_ph_met, censor.shape="|", censor.size = 4)

# Data for the full survival curve
full_plot_df_met <- p_rmst_base_met$plot$data

# Data for the shaded area (up to L)
shaded_plot_df_met <- full_plot_df_met %>% filter(time <= L)
shaded_plot_df_met <- shaded_plot_df_met %>%
  group_by(strata) %>%
  do({
    last_row <- tail(., 1)
    bind_rows(., data.frame(time = L, surv = last_row$surv, strata = last_row$strata))
  }) %>%
  arrange(strata, time)

p_rmst_met <- ggplot(full_plot_df_met, aes(x = time, y = surv, color = strata)) +
  geom_step(linewidth = 1.2) +
  geom_ribbon(data = filter(shaded_plot_df_met, strata == "group=Treatment"), aes(ymin = 0, ymax = surv), fill = "#2E9FDF", alpha = 0.3) +
  geom_ribbon(data = filter(shaded_plot_df_met, strata == "group=Control"), aes(ymin = 0, ymax = surv), fill = "#E7B800", alpha = 0.3) +
  geom_vline(xintercept = L, linetype = "dashed", color = "red") +
  annotate("text", x = L, y = 0.5, label = paste("L =", L), color = "red", angle = 90, vjust = -0.5, fontface="bold") +
  scale_color_manual(name = "Group", labels = c("Control", "Treatment"), values = c("group=Control" = "#E7B800", "group=Treatment" = "#2E9FDF")) +
  labs(title = "RMST Interpretation (PH Met)", subtitle = "Treatment effect is the difference in shaded areas.", x = "Time", y = "Survival Probability") +
  theme_classic() + theme(legend.position = "bottom") +
  geom_point(data = filter(full_plot_df_met, n.censor > 0), aes(x = time, y = surv), shape = "|", size = 4) +
  coord_cartesian(xlim = c(0, max(full_plot_df_met$time) + 5))

ggsave("images/rmst_causal_plot_ph_met.svg", plot = p_rmst_met, width = 7, height = 5)
ggsave("images/rmst_causal_plot_ph_met.png", plot = p_rmst_met, width = 7, height = 5)


# --- 4. RMST Plot for PH Assumption Violated Case ---

# Use the data from the PH Violated simulation
fit_rmst_violated <- survfit(Surv(time, status) ~ group, data = sim_data_ph_violated)
p_rmst_base_violated <- ggsurvplot(fit_rmst_violated, data = sim_data_ph_violated, censor.shape="|", censor.size = 4)

# Data for the full survival curve
full_plot_df_violated <- p_rmst_base_violated$plot$data

# Data for the shaded area (up to L)
shaded_plot_df_violated <- full_plot_df_violated %>% filter(time <= L)
shaded_plot_df_violated <- shaded_plot_df_violated %>%
  group_by(strata) %>%
  do({
    last_row <- tail(., 1)
    bind_rows(., data.frame(time = L, surv = last_row$surv, strata = last_row$strata))
  }) %>%
  arrange(strata, time)

p_rmst_violated <- ggplot(full_plot_df_violated, aes(x = time, y = surv, color = strata)) +
  geom_step(linewidth = 1.2) +
  geom_ribbon(data = filter(shaded_plot_df_violated, strata == "group=Treatment"), aes(ymin = 0, ymax = surv), fill = "#2E9FDF", alpha = 0.3) +
  geom_ribbon(data = filter(shaded_plot_df_violated, strata == "group=Control"), aes(ymin = 0, ymax = surv), fill = "#E7B800", alpha = 0.3) +
  geom_vline(xintercept = L, linetype = "dashed", color = "red") +
  annotate("text", x = L, y = 0.5, label = paste("L =", L), color = "red", angle = 90, vjust = -0.5, fontface="bold") +
  scale_color_manual(name = "Group", labels = c("Control", "Treatment"), values = c("group=Control" = "#E7B800", "group=Treatment" = "#2E9FDF")) +
  labs(title = "RMST Interpretation (PH Violated)", subtitle = "RMST remains a valid and interpretable measure.", x = "Time", y = "Survival Probability") +
  theme_classic() + theme(legend.position = "bottom") +
  geom_point(data = filter(full_plot_df_violated, n.censor > 0), aes(x = time, y = surv), shape = "|", size = 4) +
  coord_cartesian(xlim = c(0, max(full_plot_df_violated$time) + 5))

ggsave("images/rmst_causal_plot_ph_violated.svg", plot = p_rmst_violated, width = 7, height = 5)
ggsave("images/rmst_causal_plot_ph_violated.png", plot = p_rmst_violated, width = 7, height = 5)


print("All SVG and PNG images have been generated and saved in the 'images' folder.")
