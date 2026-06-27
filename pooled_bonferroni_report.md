# Pooled T-zone binomial test with Bonferroni correction for death/ERK-only sweep

## Method

For each parameter setting, the five independent replicates were pooled. If a parameter setting has five replicates of 1000 analyzed observed deaths, then the pooled sample size is N = 5000. The pooled number of T-zone deaths was tested against the null probability p0 = |T|/|W| using a one-sided exact binomial test.

The hypotheses are:

$$H_0:\pi=p_0,\qquad H_A:\pi>p_0.$$

Since 27 parameter settings were tested, Bonferroni correction was applied. With global alpha = 0.05, the Bonferroni threshold is:

$$\alpha_\mathrm{Bonferroni} = \frac{0.05}{27} = 0.00185185.$$


A parameter setting is considered significant only if:


$$p_\mathrm{pooled} < 0.00185185.$$

## Results

| parameter setting | replicates | pooled k/n | pooled fraction | pooled p-value | Bonferroni adjusted p-value | significant after Bonferroni |
|---|---:|---:|---:|---:|---:|---:|
| deathERK_ld0p5_bdR4p0_bdT1p6 | 5 | 3321/5000 | 0.6642 | 9.06e-215 | 2.45e-213 | True |
| deathERK_ld0p5_bdR4p0_bdT0p8 | 5 | 3262/5000 | 0.6524 | 3.09e-192 | 8.35e-191 | True |
| deathERK_ld0p5_bdR2p0_bdT1p6 | 5 | 3232/5000 | 0.6464 | 2.73e-181 | 7.37e-180 | True |
| deathERK_ld1p0_bdR4p0_bdT1p6 | 5 | 3231/5000 | 0.6462 | 6.23e-181 | 1.68e-179 | True |
| deathERK_ld0p5_bdR1p0_bdT1p6 | 5 | 3224/5000 | 0.6448 | 1.98e-178 | 5.35e-177 | True |
| deathERK_ld0p5_bdR4p0_bdT0p4 | 5 | 3204/5000 | 0.6408 | 2.2e-171 | 5.95e-170 | True |
| deathERK_ld0p5_bdR2p0_bdT0p4 | 5 | 3170/5000 | 0.6340 | 9.48e-160 | 2.56e-158 | True |
| deathERK_ld1p0_bdR2p0_bdT0p8 | 5 | 3165/5000 | 0.6330 | 4.48e-158 | 1.21e-156 | True |
| deathERK_ld1p0_bdR4p0_bdT0p8 | 5 | 3158/5000 | 0.6316 | 9.53e-156 | 2.57e-154 | True |
| deathERK_ld2p0_bdR4p0_bdT0p8 | 5 | 3157/5000 | 0.6314 | 2.04e-155 | 5.52e-154 | True |
| deathERK_ld0p5_bdR1p0_bdT0p8 | 5 | 3146/5000 | 0.6292 | 8.48e-152 | 2.29e-150 | True |
| deathERK_ld1p0_bdR1p0_bdT0p8 | 5 | 3128/5000 | 0.6256 | 5.64e-146 | 1.52e-144 | True |
| deathERK_ld2p0_bdR2p0_bdT1p6 | 5 | 3128/5000 | 0.6256 | 5.64e-146 | 1.52e-144 | True |
| deathERK_ld0p5_bdR2p0_bdT0p8 | 5 | 3122/5000 | 0.6244 | 4.63e-144 | 1.25e-142 | True |
| deathERK_ld1p0_bdR4p0_bdT0p4 | 5 | 3119/5000 | 0.6238 | 4.15e-143 | 1.12e-141 | True |
| deathERK_ld1p0_bdR2p0_bdT0p4 | 5 | 3101/5000 | 0.6202 | 1.82e-137 | 4.92e-136 | True |
| deathERK_ld2p0_bdR4p0_bdT0p4 | 5 | 3099/5000 | 0.6198 | 7.6e-137 | 2.05e-135 | True |
| deathERK_ld1p0_bdR1p0_bdT1p6 | 5 | 3089/5000 | 0.6178 | 9.05e-134 | 2.44e-132 | True |
| deathERK_ld1p0_bdR2p0_bdT1p6 | 5 | 3070/5000 | 0.6140 | 5.02e-128 | 1.35e-126 | True |
| deathERK_ld2p0_bdR4p0_bdT1p6 | 5 | 3066/5000 | 0.6132 | 7.81e-127 | 2.11e-125 | True |
| deathERK_ld2p0_bdR1p0_bdT1p6 | 5 | 3031/5000 | 0.6062 | 1.19e-116 | 3.22e-115 | True |
| deathERK_ld2p0_bdR2p0_bdT0p8 | 5 | 3017/5000 | 0.6034 | 1.06e-112 | 2.86e-111 | True |
| deathERK_ld2p0_bdR1p0_bdT0p8 | 5 | 3016/5000 | 0.6032 | 2.02e-112 | 5.45e-111 | True |
| deathERK_ld0p5_bdR1p0_bdT0p4 | 5 | 2984/5000 | 0.5968 | 1.13e-103 | 3.05e-102 | True |
| deathERK_ld1p0_bdR1p0_bdT0p4 | 5 | 2983/5000 | 0.5966 | 2.09e-103 | 5.65e-102 | True |
| deathERK_ld2p0_bdR2p0_bdT0p4 | 5 | 2979/5000 | 0.5958 | 2.44e-102 | 6.58e-101 | True |
| deathERK_ld2p0_bdR1p0_bdT0p4 | 5 | 2871/5000 | 0.5742 | 1.02e-75 | 2.74e-74 | True |
