# Advanced Applied Econometrics - Final Essay
# Alexander Lievens & Tobias Vanoverschelde

### Preparation
set.seed(42)

DATA_PATH <- "data/mss_repdata.dta"
dir.create("tables",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

pkgs <- c("haven",         # read .dta (Stata) files
          "tidyverse",     # data wrangling
          "plm",           # panel models + Hausman + System-GMM (pgmm)
          "AER",           # ivreg() with diagnostics
          "lmtest",        # coeftest, bptest
          "sandwich",      # robust / clustered SE
          "fixest",        # high-performance FE + IV + cluster SE
          "car",           # linearHypothesis
          "modelsummary",  # LaTeX tables
          "boot",          # bootstrap helpers
          "ggplot2"
)
invisible(lapply(pkgs, library, character.only = TRUE))

# Read Stata file and adjust to R using haven (+ everything that is not a string to plain numeric)
df <- haven::read_dta(DATA_PATH) %>%
        haven::zap_labels() %>%
        haven::zap_formats() %>%
        as_tibble()

df <- df %>% mutate(across(where(~ !is.character(.)), as.numeric))

# Build modelling data
df <- df %>%
  arrange(ccode, year) %>%
  group_by(ccode) %>%
  mutate(gdp_g_l = dplyr::lag(gdp_g)) %>%
  ungroup() %>%
  mutate(country = factor(ccode),
         year_f  = factor(year))

# Remove missings
modvars <- c("any_prio", "gdp_g", "gdp_g_l",
             "GPCP_g", "GPCP_g_l",
             "y_0", "polity2l", "ethfrac", "relfrac",
             "Oil", "lmtnest", "lpopl1")
df_est <- df %>% tidyr::drop_na(dplyr::all_of(modvars))

# Panel structure for plm
pdata <- pdata.frame(df_est, index = c("ccode", "year"))

### EDA
desc <- df_est %>%
  dplyr::select(any_prio, gdp_g, GPCP_g, y_0, polity2l, ethfrac,
                relfrac, Oil, lmtnest, lpopl1) %>%
  pivot_longer(everything(), names_to = "var", values_to = "x") %>%
  group_by(var) %>%
  summarise(N    = sum(!is.na(x)),
            mean = mean(x, na.rm = TRUE),
            sd   = sd(x,   na.rm = TRUE),
            min  = min(x,  na.rm = TRUE),
            max  = max(x,  na.rm = TRUE), .groups = "drop")
print(desc)
write.csv(desc, "tables/descriptives.csv", row.names = FALSE)

### Analysis
# Pooled OLS (LPM)
f_main <- any_prio ~ gdp_g + gdp_g_l +
                     y_0 + polity2l + ethfrac + relfrac +
                     Oil + lmtnest + lpopl1

# FE formula: drop time-invariant controls and keep only time-varying controls
f_fe_re <- any_prio ~ gdp_g + gdp_g_l + polity2l

m_pols    <- plm(f_main, data = pdata, model = "pooling")
v_pols_cl <- vcovHC(m_pols, type = "HC1", cluster = "group")
print(coeftest(m_pols, vcov = v_pols_cl))

# Heteroscedasticity diagnostics
m_pols_lm <- lm(f_main, data = df_est)
bp_test <- bptest(m_pols_lm)
wh_test <- bptest(m_pols_lm,
                  ~ fitted(m_pols_lm) + I(fitted(m_pols_lm)^2))
cat("\nBreusch-Pagan:  chi2 =", round(bp_test$statistic, 3),
    " p =", format.pval(bp_test$p.value, digits = 3), "\n")
cat("White (simplified):  chi2 =", round(wh_test$statistic, 3),
    " p =", format.pval(wh_test$p.value, digits = 3), "\n")
# LPM is inherently heteroscedastic: Var(e) = p(1-p) -> Robust SE essential

# Panel: FE and RE
# f_fe_re drops time-invariant covariates
m_fe <- plm(f_fe_re, data = pdata, model = "within", effect = "twoways")
m_re <- plm(f_fe_re, data = pdata, model = "random", effect = "twoways")

v_fe_cl <- vcovHC(m_fe, type = "HC1", cluster = "group", method = "arellano")
v_re_cl <- vcovHC(m_re, type = "HC1", cluster = "group", method = "arellano")

# Pooled OLS on the SAME trimmed specification, used as benchmark for F-test
m_pols_trim <- plm(f_fe_re, data = pdata, model = "pooling")

# (a) F test for country fixed effects (H0: c_i = 0)
F_FE  <- pFtest(m_fe, m_pols_trim)
cat("\nF-test country FE:  F =", round(F_FE$statistic, 3),
    " p =", format.pval(F_FE$p.value, digits = 3), "\n")

# (b) Hausman (FE vs RE)
H_FE_RE <- phtest(m_fe, m_re)
cat("Hausman FE vs RE:  chi2 =", round(H_FE_RE$statistic, 3),
    " p =", format.pval(H_FE_RE$p.value, digits = 3), "\n")

# Reduced form: rainfall directly on conflict
f_rf <- any_prio ~ GPCP_g + GPCP_g_l + polity2l
m_rf    <- plm(f_rf, data = pdata, model = "within", effect = "twoways")
v_rf_cl <- vcovHC(m_rf, type = "HC1", cluster = "group", method = "arellano")
print(coeftest(m_rf, vcov = v_rf_cl))

# First stage
f_fs <- gdp_g ~ GPCP_g + GPCP_g_l + polity2l
m_fs    <- plm(f_fs, data = pdata, model = "within", effect = "twoways")
v_fs_cl <- vcovHC(m_fs, type = "HC1", cluster = "group", method = "arellano")
print(coeftest(m_fs, vcov = v_fs_cl))

# F-stat on excluded instruments (cluster-robust)
F_excl <- linearHypothesis(m_fs,
                           c("GPCP_g = 0", "GPCP_g_l = 0"),
                           vcov = v_fs_cl, test = "F")
cat("\nFirst-stage F on excluded IVs (cluster-robust):  F =",
    round(F_excl$F[2], 3),
    " p =", format.pval(F_excl$`Pr(>F)`[2], digits = 3), "\n")
# Stock-Yogo 10% relative-bias critical value for 2 IVs, 1 endog = 19.93

# 2SLS
# Implementation 1: fixest with FE absorbed + cluster SE + KP F
m_iv_fx <- feols(any_prio ~ gdp_g_l + polity2l |
                            ccode + year |
                            gdp_g ~ GPCP_g + GPCP_g_l,
                  cluster = ~ccode, data = df_est)
print(summary(m_iv_fx))
fs_kp <- tryCatch(fitstat(m_iv_fx, ~ ivf1 + ivwald + sargan + wh,
                          simplify = TRUE),
                  error = function(e) NULL)
print(fs_kp)

# Implementation 2: AER::ivreg with explicit country & year dummies
m_iv_aer <- ivreg(any_prio ~ gdp_g + gdp_g_l + polity2l +
                              factor(ccode) + factor(year) |
                  GPCP_g + GPCP_g_l + gdp_g_l + polity2l +
                              factor(ccode) + factor(year),
                  data = df_est)
diag_iv <- summary(m_iv_aer, diagnostics = TRUE, vcov. = sandwich)
cat("\nAER ivreg diagnostics (HC0):\n")
print(diag_iv$diagnostics)

# Like MSS (2004): country FE + country-specific linear time trends (robustness check)
df_est <- df_est %>%
  mutate(t_lin = as.numeric(year) - min(as.numeric(year)))    # 0..17

# Build country-trend regressors as factor(ccode):t_lin interactions
m_iv_fx_tr <- feols(any_prio ~ gdp_g_l + polity2l +
                                factor(ccode):t_lin |
                                ccode |
                                gdp_g ~ GPCP_g + GPCP_g_l,
                     cluster = ~ccode, data = df_est)
cat("\n--- 2SLS with country-specific linear time trends (MSS-style) ---\n")
print(summary(m_iv_fx_tr))
fs_kp_tr <- tryCatch(fitstat(m_iv_fx_tr, ~ ivf1 + ivwald + sargan + wh,
                              simplify = TRUE),
                      error = function(e) NULL)
print(fs_kp_tr)

# First-stage with country-trends for the F-statistic
m_fs_tr <- feols(gdp_g ~ GPCP_g + GPCP_g_l + polity2l +
                          factor(ccode):t_lin |
                          ccode,
                  cluster = ~ccode, data = df_est)
F_excl_tr <- tryCatch(
  car::linearHypothesis(m_fs_tr,
                        c("GPCP_g = 0", "GPCP_g_l = 0"),
                        vcov = vcov(m_fs_tr), test = "F"),
  error = function(e) NULL)
if (!is.null(F_excl_tr)) {
  cat("\nMSS-spec first-stage F on excluded IVs:  F =",
      round(F_excl_tr$F[2], 3),
      " p =", format.pval(F_excl_tr$`Pr(>F)`[2], digits = 3), "\n")
}

# Sensitivity: alternative DV (war_prio >= 1000 deaths)
if ("war_prio" %in% names(df_est)) {
  m_iv_war <- tryCatch(
    feols(war_prio ~ gdp_g_l + polity2l |
                     ccode + year |
                     gdp_g ~ GPCP_g + GPCP_g_l,
           cluster = ~ccode, data = df_est),
    error = function(e) NULL)
  if (!is.null(m_iv_war)) {
    cat("\n--- 2SLS with alternative DV war_prio (>=1000 deaths) ---\n")
    print(summary(m_iv_war))
  }
}

# Panel balance & missingness diagnostic
balance_tab <- df_est %>%
  group_by(ccode, country_name) %>%
  summarise(N_obs = n(), year_min = min(year), year_max = max(year),
            .groups = "drop")
cat("\n--- Panel structure (obs per country) ---\n")
cat("Range of obs per country:", range(balance_tab$N_obs), "\n")
cat("Countries with full panel (T=18):",
    sum(balance_tab$N_obs == max(balance_tab$N_obs)),
    "out of", nrow(balance_tab), "\n")
print(table(balance_tab$N_obs))

# Anderson-Rubin weak-IV-robust test
Z <- cbind(df_est$GPCP_g, df_est$GPCP_g_l)
X <- model.matrix(~ gdp_g_l + polity2l + factor(ccode) + factor(year) - 1,
                  data = df_est)

ivm <- ivmodel::ivmodel(Y = df_est$any_prio, D = df_est$gdp_g, Z = Z, X = X,
                        intercept = FALSE)
ar  <- ivmodel::AR.test(ivm, beta0 = 0)

cat("Anderson-Rubin test (H0: beta_gdpg = 0):\n")
cat("  F  =", round(ar$Fstat, 3), "\n")
cat("  p  =", format.pval(ar$p.value, digits = 3), "\n")
cat("  CI = [", round(min(ar$ci), 3), ", ", round(max(ar$ci), 3), "]\n", sep = "")

# Bootstrap SE for the IV coefficient on gdp_g
B <- 500
country_ids <- unique(df_est$ccode)
n_c <- length(country_ids)
boot_b <- replicate(B, {
  picked  <- sample(country_ids, n_c, replace = TRUE)
  d_b     <- do.call(rbind, lapply(picked,
                                   function(c) df_est[df_est$ccode == c, ]))
  # Re-index resampled clusters: a country sampled twice would otherwise
  # be one collapsed group. Give every resampled country a unique id.
  d_b$rep_id <- rep(seq_along(picked),
                    times = sapply(picked,
                                   function(c) sum(df_est$ccode == c)))
  fit <- tryCatch(
    feols(any_prio ~ gdp_g_l + polity2l |
                     rep_id + year |
                     gdp_g ~ GPCP_g + GPCP_g_l,
           data = d_b, cluster = ~rep_id, warn = FALSE,
           notes = FALSE),
    error = function(e) NULL)
  if (is.null(fit)) NA_real_
  else coef(fit)["fit_gdp_g"]
})

boot_b <- boot_b[is.finite(boot_b)]
cat("\nCluster-bootstrap SE for beta_{gdp_g}:", round(sd(boot_b), 4), "\n")
cat("95% percentile CI:", round(quantile(boot_b, c(0.025, 0.975)), 4), "\n")

# IV-Probit via control function (Rivers-Vuong)
# Step 1: first-stage residuals (uses fixest with absorbed FE)
m_fs_fx <- feols(gdp_g ~ GPCP_g + GPCP_g_l + polity2l |
                          ccode + year,
                  cluster = ~ccode, data = df_est)
df_est$v_hat <- residuals(m_fs_fx)

# Step 2: probit with the first-stage residual as control function term
probit_cf <- glm(any_prio ~ gdp_g + v_hat + gdp_g_l + polity2l +
                              factor(ccode) + factor(year),
                  family = binomial(link = "probit"),
                  data = df_est)
v_hat_z  <- coef(summary(probit_cf))["v_hat", "z value"]
v_hat_p  <- coef(summary(probit_cf))["v_hat", "Pr(>|z|)"]
cat("\nControl-function endogeneity test (probit):\n",
    " z(v_hat) =", round(v_hat_z, 3),
    " p =", format.pval(v_hat_p, digits = 3),
    "  (rejects exogeneity if small)\n")

# Average marginal effect of gdp_g in probit
ame_gdp <- mean(dnorm(predict(probit_cf, type = "link"))) *
           coef(probit_cf)["gdp_g"]
cat(" Probit AME of gdp_g:", round(ame_gdp, 4),
    " (vs. LPM coef:", round(coef(m_pols)["gdp_g"], 4), ")\n")

# Dynamic panel data: System-GMM (Blundell-Bond)
# Conflict is persistent (any_prio_{t-1} likely matters). Static FE is then
# misspecified. We try Arellano-Bond and Blundell-Bond as robustness checks.
df_est <- df_est %>% arrange(ccode, year)
pdata2 <- pdata.frame(df_est, index = c("ccode", "year"))

# Difference GMM (Arellano-Bond)
m_abond <- tryCatch(
  pgmm(any_prio ~ lag(any_prio, 1) + gdp_g + gdp_g_l + polity2l |
                  lag(any_prio, 2:3) + lag(gdp_g, 2:3),
        data           = pdata2,
        effect         = "individual",
        model          = "twosteps",
        transformation = "d",
        collapse       = TRUE),
  error = function(e) { message("Arellano-Bond failed: ", e$message); NULL })

# System GMM (Blundell-Bond)
m_bbond <- tryCatch(
  pgmm(any_prio ~ lag(any_prio, 1) + gdp_g + gdp_g_l + polity2l |
                  lag(any_prio, 2:3) + lag(gdp_g, 2:3),
        data           = pdata2,
        effect         = "individual",
        model          = "twosteps",
        transformation = "ld",
        collapse       = TRUE),
  error = function(e) { message("Blundell-Bond failed: ", e$message); NULL })

cat("\n--- Arellano-Bond (difference GMM) ---\n")
print(summary(m_abond, robust = TRUE))

cat("\n--- Blundell-Bond (system GMM) ---\n")
print(summary(m_bbond, robust = TRUE))

# Generate LaTex tables
models_main <- list(
  "POLS (LPM)"     = m_pols,
  "FE 2-way"       = m_fe,
  "RE 2-way"       = m_re,
  "Reduced form"   = m_rf,
  "First stage"    = m_fs,
  "2SLS"           = m_iv_fx,
  "2SLS (trends)"  = m_iv_fx_tr
)
vlist <- list(v_pols_cl, v_fe_cl, v_re_cl, v_rf_cl, v_fs_cl,
              vcov(m_iv_fx), vcov(m_iv_fx_tr))

# Map both gdp_g and fit_gdp_g to a single row so the table is readable.
coef_map <- c(
  "gdp_g"     = "GDP growth",
  "fit_gdp_g" = "GDP growth",                    # 2SLS fitted endogenous
  "gdp_g_l"   = "GDP growth (lag)",
  "GPCP_g"    = "Rainfall growth",
  "GPCP_g_l"  = "Rainfall growth (lag)",
  "polity2l"  = "Polity2 (lag)"
)

modelsummary(models_main,
             vcov     = vlist,
             coef_map = coef_map,
             coef_omit = "factor\\(ccode\\)|y_0|ethfrac|relfrac|Oil|lmtnest|lpopl1|Intercept",
             stars    = c('*' = .1, '**' = .05, '***' = .01),
             output   = "tables/main_results.tex",
             gof_omit = "AIC|BIC|Log|RMSE|Adj|Pseudo|Std",
             title    = NULL,                    # caption set in essay.tex
             notes    = "Cluster-robust SE at country level in parentheses. Models 1, 5 and 6 use 2-way FE; model 7 uses country FE plus country-specific linear time trends. Time-invariant controls included in POLS only.",
             add_rows = data.frame(
                term = "Controls (POLS only)",
                "POLS (LPM)" = "Yes", "FE 2-way" = "—", "RE 2-way" = "—",
                "Reduced form" = "—", "First stage" = "—",
                "2SLS" = "—", "2SLS (trends)" = "—",
                check.names = FALSE))

# Dynamic-panel GMM table -----------------------------------------------------
gmm_models <- list()
gmm_models[["Arellano-Bond"]]   <- m_abond
gmm_models[["Blundell-Bond"]]   <- m_bbond
if (length(gmm_models) > 0) {
  gmm_coef_map <- c(
    "lag(any_prio, 1)" = "Conflict (lag)",
    "gdp_g"            = "GDP growth",
    "gdp_g_l"          = "GDP growth (lag)",
    "polity2l"         = "Polity2 (lag)"
  )
  modelsummary(gmm_models,
               coef_map = gmm_coef_map,
               stars    = c('*' = .1, '**' = .05, '***' = .01),
               output   = "tables/gmm_results.tex",
               gof_omit = "AIC|BIC|Log|RMSE|Adj|Pseudo|Std",
               title    = NULL,
               notes    = "Two-step robust SE. Collapsed instrument matrix with lags 2--3.")
}

# Diagnostics summary
diag_tab <- tibble(
  Statistic = c(
    "F-test country FE (H0: c_i=0)",
    "Hausman FE vs RE",
    "Breusch-Pagan heteroscedasticity",
    "First-stage F on excluded IVs (cluster-robust)",
    "First-stage F (MSS-style country trends)",
    "Stock-Yogo 10% relative-bias critical value",
    "Wu-Hausman endogeneity (AER)",
    "Sargan overid (AER)",
    "Control-function v_hat z (probit)",
    "Bootstrap 95% CI for beta_gdp_g"
  ),
  Value = c(
    round(F_FE$statistic, 3),
    round(H_FE_RE$statistic, 3),
    round(bp_test$statistic, 3),
    round(F_excl$F[2], 3),
    ifelse(is.null(F_excl_tr), "—", round(F_excl_tr$F[2], 3)),
    "19.93",
    round(diag_iv$diagnostics["Wu-Hausman", "statistic"], 3),
    round(diag_iv$diagnostics["Sargan",     "statistic"], 3),
    round(v_hat_z, 3),
    paste0("[", round(quantile(boot_b, 0.025), 3), ", ",
                  round(quantile(boot_b, 0.975), 3), "]")
  ),
  p.value = c(
    format.pval(F_FE$p.value,    digits = 3),
    format.pval(H_FE_RE$p.value, digits = 3),
    format.pval(bp_test$p.value, digits = 3),
    format.pval(F_excl$`Pr(>F)`[2], digits = 3),
    ifelse(is.null(F_excl_tr), "—",
           format.pval(F_excl_tr$`Pr(>F)`[2], digits = 3)),
    "—",
    format.pval(diag_iv$diagnostics["Wu-Hausman", "p-value"], digits = 3),
    format.pval(diag_iv$diagnostics["Sargan",     "p-value"], digits = 3),
    format.pval(v_hat_p, digits = 3),
    "—"
  )
)
write.csv(diag_tab, "tables/diagnostics.csv", row.names = FALSE)
print(diag_tab)

# Generate Figures
# (i) Conflict share over time
g1 <- df_est %>%
  group_by(year) %>%
  summarise(share = mean(any_prio, na.rm = TRUE)) %>%
  ggplot(aes(year, share)) +
    geom_line(linewidth = .7) + geom_point() +
    labs(x = "Year", y = "Share of countries in conflict (>= 25 deaths)") +
    theme_minimal()
ggsave("figures/fig1_conflict_share.pdf", g1, width = 6, height = 3.5)

# (ii) First-stage scatter
g2 <- ggplot(df_est, aes(GPCP_g, gdp_g)) +
  geom_point(alpha = .3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Rainfall growth", y = "GDP per capita growth") +
  theme_minimal()
ggsave("figures/fig2_first_stage.pdf", g2, width = 6, height = 3.5)

# (iii) Bootstrap distribution
g3 <- ggplot(tibble(b = boot_b), aes(b)) +
  geom_histogram(bins = 30, alpha = .7) +
  geom_vline(xintercept = coef(m_iv_fx)["fit_gdp_g"], linetype = 2) +
  labs(x = expression(hat(beta)[gdp_g]^IV),
       y = "Frequency",
       title = "Cluster-bootstrap distribution (B=500)") +
  theme_minimal()
ggsave("figures/fig3_bootstrap.pdf", g3, width = 6, height = 3.5)