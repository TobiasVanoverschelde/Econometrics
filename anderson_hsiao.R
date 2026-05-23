# Anderson-Hsiao FD-IV for MSS data, no AER dependency.
# First-difference equation (1) to remove country FE.
# Instrument Delta gdp_g with the lagged level gdp_g_{t-1}.
# Matches lecture 9 panel2.do.

set.seed(42)
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr)
  library(plm); library(sandwich); library(lmtest)
})

df <- haven::read_dta("data/mss_repdata.dta") |>
  haven::zap_labels() |> haven::zap_formats() |> as_tibble()
df <- df |> mutate(across(where(~ !is.character(.)), as.numeric))

df <- df |>
  arrange(ccode, year) |>
  group_by(ccode) |>
  mutate(gdp_g_l = dplyr::lag(gdp_g)) |>
  ungroup()

modvars <- c("any_prio","gdp_g","gdp_g_l","GPCP_g","GPCP_g_l",
             "y_0","polity2l","ethfrac","relfrac","Oil","lmtnest","lpopl1")
df_est <- df |> tidyr::drop_na(dplyr::all_of(modvars))

# Manual first-differences within country
df_ah <- df_est |>
  arrange(ccode, year) |>
  group_by(ccode) |>
  mutate(
    d_any_prio = any_prio - dplyr::lag(any_prio),
    d_gdp_g    = gdp_g    - dplyr::lag(gdp_g),
    d_polity2l = polity2l - dplyr::lag(polity2l),
    l_gdp_g    = dplyr::lag(gdp_g)
  ) |>
  ungroup() |>
  filter(!is.na(d_any_prio) & !is.na(d_gdp_g) &
         !is.na(l_gdp_g)    & !is.na(d_polity2l))

cat("Anderson-Hsiao sample:", nrow(df_ah), "rows,",
    length(unique(df_ah$ccode)), "countries\n\n")

# Manual 2SLS via two regressions
# Stage 1: d_gdp_g on l_gdp_g + d_polity2l
m_fs <- lm(d_gdp_g ~ l_gdp_g + d_polity2l, data = df_ah)
df_ah$d_gdp_g_hat <- fitted(m_fs)
v_fs <- sandwich::vcovCL(m_fs, cluster = df_ah$ccode, type = "HC1")
cat("--- First stage (cluster-robust SE) ---\n")
print(lmtest::coeftest(m_fs, vcov = v_fs))
F_l <- ( (coef(m_fs)["l_gdp_g"]) / sqrt(v_fs["l_gdp_g","l_gdp_g"]) )^2
cat("\nFirst-stage F on l_gdp_g (cluster-robust):", round(F_l, 3), "\n")

# Stage 2 (gives correct point estimate; SE needs correction)
m_ss <- lm(d_any_prio ~ d_gdp_g_hat + d_polity2l, data = df_ah)
cat("\n--- Second-stage point estimates (SE need correction) ---\n")
print(coef(m_ss))

# Proper 2SLS via plm
pdata_ah <- pdata.frame(df_est, index = c("ccode","year"))
m_iv <- plm(any_prio ~ gdp_g + polity2l |
                       polity2l + lag(gdp_g),
            data = pdata_ah, model = "fd")
v_iv <- vcovHC(m_iv, type = "HC1", cluster = "group", method = "arellano")
cat("\n--- plm FD-IV (Anderson-Hsiao), cluster-robust SE ---\n")
print(lmtest::coeftest(m_iv, vcov = v_iv))
cat("\nN observations:", nobs(m_iv), "\n")
