library(dplyr)
library(tidyr)

# reshape to long format
boot_long <- boot_df_wide %>%
  pivot_longer(cols = starts_with("boot_"),
               names_to = "replicate",
               values_to = "value") %>%
  filter(grepl("^boot_[0-9]+$", replicate))  # keep only boot_1 ... boot_10

# summarize variation per quadrat
boot_summary <- boot_long %>%
  group_by(Name) %>%
  summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value),
    cv = sd / mean,
    ci_lower = mean - qt(0.975, df = n - 1) * sd / sqrt(n),
    ci_upper = mean + qt(0.975, df = n - 1) * sd / sqrt(n),
    ci_width = ci_upper - ci_lower
  )

# define thresholds (you can adjust these)
cv_threshold <- 0.1        # 10% variation
ci_width_threshold <- 0.2  # depends on your scale

boot_summary <- boot_summary %>%
  mutate(
    flag_high_cv = cv > cv_threshold,
    flag_wide_ci = ci_width > ci_width_threshold,
    unreliable = flag_high_cv | flag_wide_ci
  )

# view problematic quadrats
problem_quadrats <- boot_summary %>%
  filter(unreliable)

print(problem_quadrats)
