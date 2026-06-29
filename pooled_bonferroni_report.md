# Pooled T-zone binomial test with Bonferroni correction: rejection-sampling version

## Method

For each parameter setting, independent replicates are pooled. If a parameter setting has five replicates of 1000 analyzed observed deaths, then the pooled sample size is N = 5000. The pooled T-zone count is tested against p0 = |T|/|W| by a one-sided exact binomial test.

The hypotheses are:

```latex
H_0:\pi=p_0,\qquad H_A:\pi>p_0.
```

Since 27 parameter settings are tested, the Bonferroni threshold at global alpha=0.05 is alpha_Bonferroni = 0.00185185.

## Results

| parameter setting | replicates | pooled k/n | pooled fraction | pooled density ratio | pooled p-value | Bonferroni adjusted p-value | significant |
|---|---:|---:|---:|---:|---:|---:|---:|
| stage2_deathERK_ld1p0_bdR4p0_bdT1p6 | 5 | 3222/5000 | 0.6444 | 2.2652 | 1.02e-177 | 2.75e-176 | True |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p8 | 5 | 3217/5000 | 0.6434 | 2.2553 | 6.03e-176 | 1.63e-174 | True |
| stage2_deathERK_ld2p0_bdR2p0_bdT1p6 | 5 | 3210/5000 | 0.6420 | 2.2416 | 1.76e-173 | 4.75e-172 | True |
| stage2_deathERK_ld0p5_bdR1p0_bdT1p6 | 5 | 3201/5000 | 0.6402 | 2.2242 | 2.44e-170 | 6.59e-169 | True |
| stage2_deathERK_ld0p5_bdR4p0_bdT1p6 | 5 | 3181/5000 | 0.6362 | 2.1860 | 1.82e-163 | 4.92e-162 | True |
| stage2_deathERK_ld1p0_bdR2p0_bdT1p6 | 5 | 3168/5000 | 0.6336 | 2.1616 | 4.44e-159 | 1.2e-157 | True |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p8 | 5 | 3165/5000 | 0.6330 | 2.1560 | 4.48e-158 | 1.21e-156 | True |
| stage2_deathERK_ld0p5_bdR2p0_bdT1p6 | 5 | 3160/5000 | 0.6320 | 2.1467 | 2.07e-156 | 5.59e-155 | True |
| stage2_deathERK_ld2p0_bdR4p0_bdT1p6 | 5 | 3156/5000 | 0.6312 | 2.1394 | 4.38e-155 | 1.18e-153 | True |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p4 | 5 | 3148/5000 | 0.6296 | 2.1247 | 1.88e-152 | 5.07e-151 | True |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p8 | 5 | 3147/5000 | 0.6294 | 2.1229 | 3.99e-152 | 1.08e-150 | True |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p4 | 5 | 3145/5000 | 0.6290 | 2.1193 | 1.8e-151 | 4.85e-150 | True |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p8 | 5 | 3143/5000 | 0.6286 | 2.1156 | 8.07e-151 | 2.18e-149 | True |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p4 | 5 | 3132/5000 | 0.6264 | 2.0958 | 2.93e-147 | 7.92e-146 | True |
| stage2_deathERK_ld1p0_bdR1p0_bdT1p6 | 5 | 3127/5000 | 0.6254 | 2.0869 | 1.18e-145 | 3.18e-144 | True |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p4 | 5 | 3118/5000 | 0.6236 | 2.0709 | 8.6e-143 | 2.32e-141 | True |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p8 | 5 | 3115/5000 | 0.6230 | 2.0656 | 7.62e-142 | 2.06e-140 | True |
| stage2_deathERK_ld2p0_bdR1p0_bdT1p6 | 5 | 3108/5000 | 0.6216 | 2.0534 | 1.2e-139 | 3.25e-138 | True |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p8 | 5 | 3093/5000 | 0.6186 | 2.0274 | 5.38e-135 | 1.45e-133 | True |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p8 | 5 | 3072/5000 | 0.6144 | 1.9917 | 1.26e-128 | 3.42e-127 | True |
| stage2_deathERK_ld0p5_bdR1p0_bdT0p8 | 5 | 3069/5000 | 0.6138 | 1.9867 | 9.98e-128 | 2.69e-126 | True |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p8 | 5 | 3022/5000 | 0.6044 | 1.9098 | 4.2e-114 | 1.13e-112 | True |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p4 | 5 | 3020/5000 | 0.6040 | 1.9066 | 1.53e-113 | 4.14e-112 | True |
| stage2_deathERK_ld0p5_bdR1p0_bdT0p4 | 5 | 3005/5000 | 0.6010 | 1.8828 | 2.26e-109 | 6.1e-108 | True |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p4 | 5 | 2958/5000 | 0.5916 | 1.8107 | 7.75e-97 | 2.09e-95 | True |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p4 | 5 | 2929/5000 | 0.5858 | 1.7679 | 1.69e-89 | 4.57e-88 | True |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p4 | 5 | 2902/5000 | 0.5804 | 1.7290 | 6.18e-83 | 1.67e-81 | True |
