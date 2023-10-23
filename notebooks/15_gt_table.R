library(readr)
library(writexl)
library(gt)
library(tidyr)

path_50 = "/Users/niccolo/Desktop/r_projects/efficient-pairwise-likelihood/data/simres_50it_2xdata.rds"
path_100 =  "~/Desktop/r_projects/efficient-pairwise-likelihood/data/simres_100it_2xdata.rds"

simres =read_rds(path_100)


# for GPT commenting
# write_xlsx(simres, "/Users/niccolo/Desktop/simres.xlsx")

simres_tab = simres %>%
  mutate(across(is.numeric, ~round(., 5))) %>%
  mutate(phi = if_else(phi == 1, "φ	 = 1",  "φ = 0.8")) %>%
  gt(
    rowname_col = "sample_size",
    groupname_col = "phi"
  )  %>%
  cols_hide(columns = c(true_sigma, true_psi, true_beta, time_elapsed)) %>%
  tab_header(
    title = md("**MC** simulation results"),
    subtitle = md("100 iterations, **PL** (CKD-TPL) **FL** (Full Likelihood)")
  )   %>%
  cols_label(
    # true_sigma = "{{:sigma: = 1}}",
    # true_psi = "ψ",
    # true_beta = "β = 1",
    # sigma_hat_fullLik = "σ_FL",
    relative_bias_beta_hat_pairwise = md("β&#x0302;"),
    mse_beta_hat_pairwise = md("β&#x0302;"),
    relative_bias_beta_hat_fullLik = md("RB β&#x0302;<sub>FL</sub>"),
    relative_bias_sigma_hat_sq_pairwise = md("σ&#x0302;"),
    mse_sigma_hat_sq_pairwise = md("σ&#x0302;"),
    relative_bias_sigma_hat_sq_fullLik = md("RB σ&#x0302;<sub>PL</sub>"),
    relative_bias_psi_hat_pairwise = md("ψ&#x0302;"),
    mse_psi_hat_pairwise = md("ψ&#x0302;"),
    relative_bias_psi_hat_fullLik = md("RB ψ&#x0302;<sub>FL</sub>"),
    sigma_hat_pairwise = md("σ&#x0302;<sub>PL</sub>"),
    sigma_hat_fullLik = md("σ&#x0302;<sub>FL</sub>"),
    psi_hat_pairwise = md("ψ&#x0302;<sub>PL</sub>"),
    psi_hat_fullLik = md("ψ&#x0302;<sub>FL</sub>"),
    beta_hat_pairwise = md("β&#x0302;<sub>PL</sub>"),
    beta_hat_fullLik = md("β&#x0302;<sub>FL</sub>"),
  ) %>%
  tab_spanner(
    label = "σ = 1",
    columns = c("sigma_hat_pairwise", "sigma_hat_fullLik"),
  ) %>%
  tab_spanner(
    label = "β = 1",
    columns = c("beta_hat_pairwise", "beta_hat_fullLik"),
  ) %>%
  tab_spanner(
    label = "ψ",
    columns = c("psi_hat_pairwise", "psi_hat_fullLik"),
  ) %>%
  # tab_spanner(
  #   label = md("Parameter"),
  #   columns = c("relative_bias_beta_hat_pairwise", "relative_bias_beta_hat_fullLik")
  # ) %>%
  tab_spanner(
    label = md("Relative Bias (**PL**, **FL**)"),
    columns = c("relative_bias_beta_hat_pairwise", "relative_bias_beta_hat_fullLik",
                "relative_bias_sigma_hat_sq_pairwise", "relative_bias_sigma_hat_sq_fullLik",
                "relative_bias_psi_hat_pairwise", "relative_bias_psi_hat_fullLik"
    )
  ) %>%
  cols_merge(
    columns = c(relative_bias_beta_hat_pairwise, relative_bias_beta_hat_fullLik), pattern = "({1}, {2})"
  ) %>%
  cols_merge(
    columns = c(relative_bias_sigma_hat_sq_pairwise, relative_bias_sigma_hat_sq_fullLik), pattern = "({1}, {2})"
  ) %>%
  cols_merge(
    columns = c(relative_bias_psi_hat_pairwise, relative_bias_psi_hat_fullLik), pattern = "({1}, {2})"
  ) %>%
  tab_spanner(
    label = md("MSE (**PL**, **FL**)"),
    columns = c("mse_beta_hat_pairwise", "mse_beta_hat_fullLik",
                "mse_sigma_hat_sq_pairwise", "mse_sigma_hat_sq_fuillLik",
                "mse_psi_hat_pairwise", "mse_psi_hat_fullLik"
    )
  ) %>%
  cols_merge(
    columns = c(mse_beta_hat_pairwise, mse_beta_hat_fullLik), pattern = "({1}, {2})"
  ) %>%
  cols_merge(
    columns = c(mse_sigma_hat_sq_pairwise, mse_sigma_hat_sq_fuillLik), pattern = "({1}, {2})"
  ) %>%
  cols_merge(
    columns = c(mse_psi_hat_pairwise, mse_psi_hat_fullLik), pattern = "({1}, {2})"
  ) %>%
  cols_align(
    align ="center",
    columns = everything()
  )

gtsave(data = simres_tab, filename = here("img", "tables", "simres_tab_100it_2xdata.png"))




## seconda tabella con simulazione BUFFER ----
path = "/Users/niccolo/Desktop/r_projects/efficient-pairwise-likelihood/data/simres_buff2023-10-16.rds"

simres_buff=read_rds(path)


simres_buff_tab = simres_buff %>%
  mutate(across(is.numeric, ~round(., 5))) %>%
  mutate(phi = if_else(phi == 1, "φ	 = 1",  "φ = 0.8")) %>%
  gt(
    rowname_col = "sample_size",
    groupname_col = "phi"
  )  %>%
  cols_hide(columns = c(true_sigma, true_psi, true_beta, time_elapsed)) %>%
  tab_header(
    title = md("**MC** simulation results"),
    subtitle = md("50 iterations, **PL** (CKD-TPL) **FL** (Full Likelihood)")
  )   %>%
  cols_label(
    # true_sigma = "{{:sigma: = 1}}",
    # true_psi = "ψ",
    # true_beta = "β = 1",
    # sigma_hat_fullLik = "σ_FL",
    relative_bias_beta_hat_pairwise = md("β&#x0302;"),
    mse_beta_hat_pairwise = md("β&#x0302;"),
    relative_bias_beta_hat_fullLik = md("RB β&#x0302;<sub>FL</sub>"),
    relative_bias_sigma_hat_sq_pairwise = md("σ&#x0302;"),
    mse_sigma_hat_sq_pairwise = md("σ&#x0302;"),
    relative_bias_sigma_hat_sq_fullLik = md("RB σ&#x0302;<sub>PL</sub>"),
    relative_bias_psi_hat_pairwise = md("ψ&#x0302;"),
    mse_psi_hat_pairwise = md("ψ&#x0302;"),
    relative_bias_psi_hat_fullLik = md("RB ψ&#x0302;<sub>FL</sub>"),
    sigma_hat_pairwise = md("σ&#x0302;<sub>PL</sub>"),
    sigma_hat_fullLik = md("σ&#x0302;<sub>FL</sub>"),
    psi_hat_pairwise = md("ψ&#x0302;<sub>PL</sub>"),
    psi_hat_fullLik = md("ψ&#x0302;<sub>FL</sub>"),
    beta_hat_pairwise = md("β&#x0302;<sub>PL</sub>"),
    beta_hat_fullLik = md("β&#x0302;<sub>FL</sub>"),
  ) %>%
  tab_spanner(
    label = "σ = 1",
    columns = c("sigma_hat_pairwise", "sigma_hat_fullLik"),
  ) %>%
  tab_spanner(
    label = "β = 1",
    columns = c("beta_hat_pairwise", "beta_hat_fullLik"),
  ) %>%
  tab_spanner(
    label = "ψ",
    columns = c("psi_hat_pairwise", "psi_hat_fullLik"),
  ) %>%
  tab_spanner(
    label = md("Relative Bias (**PL**, **FL**)"),
    columns = c("relative_bias_beta_hat_pairwise", "relative_bias_beta_hat_fullLik",
                "relative_bias_sigma_hat_sq_pairwise", "relative_bias_sigma_hat_sq_fullLik",
                "relative_bias_psi_hat_pairwise", "relative_bias_psi_hat_fullLik"
    )
  ) %>%
  cols_merge(
    columns = c(relative_bias_beta_hat_pairwise, relative_bias_beta_hat_fullLik), pattern = "({1}, {2})"
  ) %>%
  cols_merge(
    columns = c(relative_bias_sigma_hat_sq_pairwise, relative_bias_sigma_hat_sq_fullLik), pattern = "({1}, {2})"
  ) %>%
  cols_merge(
    columns = c(relative_bias_psi_hat_pairwise, relative_bias_psi_hat_fullLik), pattern = "({1}, {2})"
  ) %>%
  tab_spanner(
    label = md("MSE (**PL**, **FL**)"),
    columns = c("mse_beta_hat_pairwise", "mse_beta_hat_fullLik",
                "mse_sigma_hat_sq_pairwise", "mse_sigma_hat_sq_fuillLik",
                "mse_psi_hat_pairwise", "mse_psi_hat_fullLik"
    )
  ) %>%
  cols_merge(
    columns = c(mse_beta_hat_pairwise, mse_beta_hat_fullLik), pattern = "({1}, {2})"
  ) %>%
  cols_merge(
    columns = c(mse_sigma_hat_sq_pairwise, mse_sigma_hat_sq_fuillLik), pattern = "({1}, {2})"
  ) %>%
  cols_merge(
    columns = c(mse_psi_hat_pairwise, mse_psi_hat_fullLik), pattern = "({1}, {2})"
  ) %>%
  cols_align(
    align ="center",
    columns = everything()
  )

gtsave(data = simres_buff_tab, filename = here("img", "tables", "simres_buff_tab.png"))



