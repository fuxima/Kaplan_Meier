set.seed(123) # For reproducibility

# Generate synthetic lung cancer survival data
n <- 300 # Number of patients
synthetic_lung <- data.frame(
  patient_id = 1:n,
  time_days = round(rweibull(n, shape = 0.8, scale = 300)), # Time in days
  status = rbinom(n, size = 1, prob = 0.7),
  sex = factor(sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4)),
               levels = 1:2, labels = c("Male", "Female")),
  age = round(rnorm(n, mean = 62, sd = 8)),
  ph.ecog = sample(0:4, n, replace = TRUE, prob = c(0.4, 0.3, 0.2, 0.08, 0.02)),
  ph.karno = round(runif(n, min = 30, max = 100)),
  meal.cal = round(rnorm(n, mean = 1000, sd = 300)),
  wt.loss = round(abs(rnorm(n, mean = 5, sd = 7)))
)

# Convert time from days to years (365.25 days/year)
synthetic_lung$time_years <- synthetic_lung$time_days / 365.25

# Ensure survival times are positive
synthetic_lung$time_years <- abs(synthetic_lung$time_years)

# Add censoring (some patients still alive at end of study)
synthetic_lung$status[synthetic_lung$time_years > quantile(synthetic_lung$time_years, 0.8)] <- 0

# Calculate Odds Ratio for sex difference
library(epitools)
or_table <- table(synthetic_lung$sex, synthetic_lung$status)
or_result <- oddsratio(or_table)
cat("\nOdds Ratio (Female vs Male):\n")
print(or_result$measure[2,]) # Shows OR, lower CI, upper CI

# Survival analysis with time in years
library(survival)
library(survminer)

# Fit Cox model to get Hazard Ratio
cox_model <- coxph(Surv(time_years, status) ~ sex, data = synthetic_lung)
cox_summary <- summary(cox_model)
cat("\nCox Model Hazard Ratio:\n")
print(cox_summary$conf.int[, c(1,3,4)])

surv_plot <- ggsurvplot(
  survfit(Surv(time_years, status) ~ sex, data = synthetic_lung),
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  xlab = "Time (years)",
  ylab = "Survival Probability",
  legend.labs = c("Male", "Female"),
  title = "Survival by Sex (Time in Years)",
  subtitle = paste("OR =", round(or_result$measure[2,1], 2),
                   "[", round(or_result$measure[2,2], 2), "-",
                   round(or_result$measure[2,3], 2), "]",
                   " | HR =", round(cox_summary$conf.int[1], 2),
                   "[", round(cox_summary$conf.int[3], 2), "-",
                   round(cox_summary$conf.int[4], 2), "]"),
  palette = c("#1f77b4", "#ff7f0e"),
  break.time.by = 1,
  ggtheme = theme_minimal()
)

print(surv_plot)
