# Diagnostic experiments for T-zone mechanisms

## Experiment 1A: fixed 1000 deaths

This experiment sweeps lambda_d across three ERK regimes while keeping the activation regime fixed at
lambda_a_T=0.5, lambda_a_c=0.005, beta_a_R=2.5, beta_a_T=3.0.
Each run uses the original stopping rule: burn-in 500 deaths, then analyze the next 1000 deaths.
Use experiment1A_fixed_1000_deaths_group_summary.csv and experiment1_spearman_trend_tests.csv to check whether lambda_d has a stable main effect across ERK regimes.

## Experiment 1B: fixed physical time

The fixed physical observation window was estimated from Experiment 1A at the reference ERK regime and lambda_d=1.0:
T_obs = 61.349906 time units after burn-in.
All lambda_d values are then compared over the same physical time window.
If lambda_d mainly changes the speed of reaching 1000 deaths, its effect on T-zone fraction should be weaker here than in Experiment 1A,
while the number of observed deaths should increase with lambda_d.

## Experiment 2: activation localization x death/ERK robustness

This experiment keeps lambda_a_T=0.5 and lambda_a_c=0.005 fixed, then compares three activation-localization regimes:
non-local (beta_a_R=1.0, beta_a_T=1.2), medium (beta_a_R=2.5, beta_a_T=3.0), and local (beta_a_R=5.0, beta_a_T=3.0).
Inside each activation regime, it sweeps the same lambda_d values and ERK regimes.
Use experiment2_activation_regime_robustness_summary.csv and experiment2_activation_localization_deathERK_pooled_bonferroni.csv to decide whether
activation localization controls the existence of T-zone enrichment, while death/ERK parameters modulate strength and speed.

## Key output files

- experiment1A_fixed_1000_deaths_replicates.csv
- experiment1A_fixed_1000_deaths_group_summary.csv
- experiment1A_fixed_1000_deaths_pooled_bonferroni.csv
- experiment1B_fixed_physical_time_replicates.csv
- experiment1B_fixed_physical_time_group_summary.csv
- experiment1B_fixed_physical_time_pooled_bonferroni.csv
- experiment1_spearman_trend_tests.csv
- experiment2_activation_localization_deathERK_replicates.csv
- experiment2_activation_localization_deathERK_group_summary.csv
- experiment2_activation_localization_deathERK_pooled_bonferroni.csv
- experiment2_activation_regime_robustness_summary.csv
- experiment2_lambda_d_trend_tests_by_activation_and_ERK.csv

