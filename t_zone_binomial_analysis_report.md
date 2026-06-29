# T-zone binomial enrichment analysis after burn-in: rejection-sampling version

All 2D histograms use the same colorbar scale: 0 to 13 observed deaths per bin.

Burn-in removes the initial transient period: the first 500 observed deaths are simulated but not used for statistics. The reported statistics use the next 1000 observed deaths.

## Interpretation rules

- `t_zone_density_ratio`: values > 1 mean deaths are denser in the T-zone than outside.
- `t_zone_binomial_p_value_greater`: one-sided exact binomial p-value for `H0: pi = |T|/|W|` against `HA: pi > |T|/|W|`.
- `pooled_binomial_p_value_greater`: same test after pooling independent replicates of the same parameter setting.
- Old `spatial CV` and `chi-square uniformity` are intentionally not used as the central evidence, because they test global non-uniformity rather than enrichment in the fixed T-zone.

## Replicate-level results

| scenario | time for analyzed 1000 deaths | T-zone deaths k/n | area fraction p0 | observed fraction | binomial p-value greater | T-zone density ratio | death rejection | by ERK |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| stage2_deathERK_ld0p5_bdR1p0_bdT0p4__rep01 | 209.148 | 597/1000 | 0.4444 | 0.5970 | 2.64e-22 | 1.8517 | 0.9051 | 0.5037 |
| stage2_deathERK_ld0p5_bdR1p0_bdT0p4__rep02 | 218.612 | 599/1000 | 0.4444 | 0.5990 | 7.57e-23 | 1.8672 | 0.9084 | 0.4899 |
| stage2_deathERK_ld0p5_bdR1p0_bdT0p4__rep03 | 202.693 | 611/1000 | 0.4444 | 0.6110 | 2.97e-26 | 1.9634 | 0.8986 | 0.4803 |
| stage2_deathERK_ld0p5_bdR1p0_bdT0p4__rep04 | 193.016 | 582/1000 | 0.4444 | 0.5820 | 1.84e-18 | 1.7404 | 0.8959 | 0.5230 |
| stage2_deathERK_ld0p5_bdR1p0_bdT0p4__rep05 | 204.874 | 616/1000 | 0.4444 | 0.6160 | 9.49e-28 | 2.0052 | 0.9042 | 0.4915 |
| stage2_deathERK_ld0p5_bdR1p0_bdT0p8__rep01 | 144.353 | 615/1000 | 0.4444 | 0.6150 | 1.91e-27 | 1.9968 | 0.8707 | 0.3811 |
| stage2_deathERK_ld0p5_bdR1p0_bdT0p8__rep02 | 151.780 | 581/1000 | 0.4444 | 0.5810 | 3.2e-18 | 1.7333 | 0.8657 | 0.3840 |
| stage2_deathERK_ld0p5_bdR1p0_bdT0p8__rep03 | 146.315 | 643/1000 | 0.4444 | 0.6430 | 1.28e-36 | 2.2514 | 0.8545 | 0.3563 |
| stage2_deathERK_ld0p5_bdR1p0_bdT0p8__rep04 | 153.596 | 603/1000 | 0.4444 | 0.6030 | 5.93e-24 | 1.8986 | 0.8610 | 0.4014 |
| stage2_deathERK_ld0p5_bdR1p0_bdT0p8__rep05 | 157.060 | 627/1000 | 0.4444 | 0.6270 | 3.35e-31 | 2.1012 | 0.8706 | 0.3490 |
| stage2_deathERK_ld0p5_bdR1p0_bdT1p6__rep01 | 125.704 | 662/1000 | 0.4444 | 0.6620 | 1.11e-43 | 2.4482 | 0.8361 | 0.2438 |
| stage2_deathERK_ld0p5_bdR1p0_bdT1p6__rep02 | 112.456 | 640/1000 | 0.4444 | 0.6400 | 1.44e-35 | 2.2222 | 0.8245 | 0.2554 |
| stage2_deathERK_ld0p5_bdR1p0_bdT1p6__rep03 | 115.518 | 649/1000 | 0.4444 | 0.6490 | 8.91e-39 | 2.3113 | 0.8224 | 0.2514 |
| stage2_deathERK_ld0p5_bdR1p0_bdT1p6__rep04 | 119.180 | 624/1000 | 0.4444 | 0.6240 | 3.08e-30 | 2.0745 | 0.8242 | 0.2941 |
| stage2_deathERK_ld0p5_bdR1p0_bdT1p6__rep05 | 112.223 | 626/1000 | 0.4444 | 0.6260 | 7.06e-31 | 2.0922 | 0.8309 | 0.2443 |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p4__rep01 | 125.828 | 648/1000 | 0.4444 | 0.6480 | 2.06e-38 | 2.3011 | 0.8336 | 0.2601 |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p4__rep02 | 119.140 | 629/1000 | 0.4444 | 0.6290 | 7.48e-32 | 2.1193 | 0.8308 | 0.2941 |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p4__rep03 | 111.436 | 619/1000 | 0.4444 | 0.6190 | 1.14e-28 | 2.0308 | 0.8153 | 0.2694 |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p4__rep04 | 116.607 | 607/1000 | 0.4444 | 0.6070 | 4.34e-25 | 1.9307 | 0.8174 | 0.3024 |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p4__rep05 | 109.262 | 645/1000 | 0.4444 | 0.6450 | 2.48e-37 | 2.2711 | 0.8229 | 0.2323 |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p8__rep01 | 105.437 | 654/1000 | 0.4444 | 0.6540 | 1.26e-40 | 2.3627 | 0.8043 | 0.1646 |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p8__rep02 | 97.747 | 644/1000 | 0.4444 | 0.6440 | 5.64e-37 | 2.2612 | 0.7928 | 0.1872 |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p8__rep03 | 91.870 | 603/1000 | 0.4444 | 0.6030 | 5.93e-24 | 1.8986 | 0.7826 | 0.1703 |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p8__rep04 | 89.618 | 609/1000 | 0.4444 | 0.6090 | 1.15e-25 | 1.9469 | 0.7750 | 0.1883 |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p8__rep05 | 90.351 | 637/1000 | 0.4444 | 0.6370 | 1.56e-34 | 2.1935 | 0.7937 | 0.1361 |
| stage2_deathERK_ld0p5_bdR2p0_bdT1p6__rep01 | 95.013 | 663/1000 | 0.4444 | 0.6630 | 4.5e-44 | 2.4592 | 0.7829 | 0.0969 |
| stage2_deathERK_ld0p5_bdR2p0_bdT1p6__rep02 | 88.802 | 664/1000 | 0.4444 | 0.6640 | 1.82e-44 | 2.4702 | 0.7672 | 0.1068 |
| stage2_deathERK_ld0p5_bdR2p0_bdT1p6__rep03 | 83.381 | 586/1000 | 0.4444 | 0.5860 | 1.9e-19 | 1.7693 | 0.7624 | 0.1090 |
| stage2_deathERK_ld0p5_bdR2p0_bdT1p6__rep04 | 79.716 | 625/1000 | 0.4444 | 0.6250 | 1.48e-30 | 2.0833 | 0.7488 | 0.1193 |
| stage2_deathERK_ld0p5_bdR2p0_bdT1p6__rep05 | 83.396 | 622/1000 | 0.4444 | 0.6220 | 1.32e-29 | 2.0569 | 0.7818 | 0.0906 |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p4__rep01 | 90.421 | 648/1000 | 0.4444 | 0.6480 | 2.06e-38 | 2.3011 | 0.7710 | 0.1079 |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p4__rep02 | 87.739 | 660/1000 | 0.4444 | 0.6600 | 6.6e-43 | 2.4265 | 0.7649 | 0.1169 |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p4__rep03 | 81.286 | 603/1000 | 0.4444 | 0.6030 | 5.93e-24 | 1.8986 | 0.7548 | 0.1052 |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p4__rep04 | 77.280 | 613/1000 | 0.4444 | 0.6130 | 7.59e-27 | 1.9800 | 0.7390 | 0.1381 |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p4__rep05 | 80.392 | 594/1000 | 0.4444 | 0.5940 | 1.67e-21 | 1.8288 | 0.7722 | 0.0905 |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p8__rep01 | 86.367 | 660/1000 | 0.4444 | 0.6600 | 6.6e-43 | 2.4265 | 0.7590 | 0.0594 |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p8__rep02 | 79.232 | 682/1000 | 0.4444 | 0.6820 | 7.19e-52 | 2.6808 | 0.7489 | 0.0663 |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p8__rep03 | 76.139 | 598/1000 | 0.4444 | 0.5980 | 1.42e-22 | 1.8595 | 0.7368 | 0.0638 |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p8__rep04 | 69.473 | 621/1000 | 0.4444 | 0.6210 | 2.73e-29 | 2.0482 | 0.7139 | 0.0787 |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p8__rep05 | 78.739 | 604/1000 | 0.4444 | 0.6040 | 3.1e-24 | 1.9066 | 0.7702 | 0.0486 |
| stage2_deathERK_ld0p5_bdR4p0_bdT1p6__rep01 | 80.591 | 666/1000 | 0.4444 | 0.6660 | 2.94e-45 | 2.4925 | 0.7474 | 0.0335 |
| stage2_deathERK_ld0p5_bdR4p0_bdT1p6__rep02 | 74.871 | 692/1000 | 0.4444 | 0.6920 | 2.92e-56 | 2.8084 | 0.7303 | 0.0350 |
| stage2_deathERK_ld0p5_bdR4p0_bdT1p6__rep03 | 73.570 | 606/1000 | 0.4444 | 0.6060 | 8.39e-25 | 1.9226 | 0.7259 | 0.0301 |
| stage2_deathERK_ld0p5_bdR4p0_bdT1p6__rep04 | 64.422 | 616/1000 | 0.4444 | 0.6160 | 9.49e-28 | 2.0052 | 0.6956 | 0.0420 |
| stage2_deathERK_ld0p5_bdR4p0_bdT1p6__rep05 | 73.617 | 601/1000 | 0.4444 | 0.6010 | 2.14e-23 | 1.8828 | 0.7539 | 0.0258 |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p4__rep01 | 164.329 | 574/1000 | 0.4444 | 0.5740 | 1.41e-16 | 1.6843 | 0.9385 | 0.6043 |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p4__rep02 | 177.809 | 633/1000 | 0.4444 | 0.6330 | 3.54e-33 | 2.1560 | 0.9434 | 0.5648 |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p4__rep03 | 141.778 | 589/1000 | 0.4444 | 0.5890 | 3.32e-20 | 1.7914 | 0.9323 | 0.5876 |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p4__rep04 | 151.158 | 562/1000 | 0.4444 | 0.5620 | 5.89e-14 | 1.6039 | 0.9305 | 0.6005 |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p4__rep05 | 163.815 | 571/1000 | 0.4444 | 0.5710 | 6.74e-16 | 1.6638 | 0.9410 | 0.6460 |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p8__rep01 | 114.698 | 588/1000 | 0.4444 | 0.5880 | 5.96e-20 | 1.7840 | 0.9095 | 0.4784 |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p8__rep02 | 119.846 | 654/1000 | 0.4444 | 0.6540 | 1.26e-40 | 2.3627 | 0.9169 | 0.4210 |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p8__rep03 | 111.350 | 607/1000 | 0.4444 | 0.6070 | 4.34e-25 | 1.9307 | 0.9032 | 0.4515 |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p8__rep04 | 109.224 | 584/1000 | 0.4444 | 0.5840 | 5.95e-19 | 1.7548 | 0.9075 | 0.5113 |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p8__rep05 | 104.909 | 589/1000 | 0.4444 | 0.5890 | 3.32e-20 | 1.7914 | 0.9052 | 0.4514 |
| stage2_deathERK_ld1p0_bdR1p0_bdT1p6__rep01 | 78.831 | 611/1000 | 0.4444 | 0.6110 | 2.97e-26 | 1.9634 | 0.8680 | 0.3135 |
| stage2_deathERK_ld1p0_bdR1p0_bdT1p6__rep02 | 89.506 | 639/1000 | 0.4444 | 0.6390 | 3.2e-35 | 2.2126 | 0.8854 | 0.3300 |
| stage2_deathERK_ld1p0_bdR1p0_bdT1p6__rep03 | 80.718 | 627/1000 | 0.4444 | 0.6270 | 3.35e-31 | 2.1012 | 0.8756 | 0.3280 |
| stage2_deathERK_ld1p0_bdR1p0_bdT1p6__rep04 | 71.020 | 599/1000 | 0.4444 | 0.5990 | 7.57e-23 | 1.8672 | 0.8613 | 0.3688 |
| stage2_deathERK_ld1p0_bdR1p0_bdT1p6__rep05 | 80.801 | 651/1000 | 0.4444 | 0.6510 | 1.64e-39 | 2.3317 | 0.8739 | 0.3138 |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p4__rep01 | 87.543 | 595/1000 | 0.4444 | 0.5950 | 9.06e-22 | 1.8364 | 0.8738 | 0.3742 |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p4__rep02 | 85.691 | 604/1000 | 0.4444 | 0.6040 | 3.1e-24 | 1.9066 | 0.8823 | 0.3304 |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p4__rep03 | 80.079 | 609/1000 | 0.4444 | 0.6090 | 1.15e-25 | 1.9469 | 0.8726 | 0.3582 |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p4__rep04 | 72.140 | 593/1000 | 0.4444 | 0.5930 | 3.06e-21 | 1.8213 | 0.8626 | 0.4120 |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p4__rep05 | 88.171 | 619/1000 | 0.4444 | 0.6190 | 1.14e-28 | 2.0308 | 0.8841 | 0.3649 |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p8__rep01 | 71.474 | 639/1000 | 0.4444 | 0.6390 | 3.2e-35 | 2.2126 | 0.8496 | 0.2470 |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p8__rep02 | 65.215 | 623/1000 | 0.4444 | 0.6230 | 6.4e-30 | 2.0656 | 0.8415 | 0.2310 |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p8__rep03 | 63.001 | 614/1000 | 0.4444 | 0.6140 | 3.81e-27 | 1.9883 | 0.8348 | 0.2193 |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p8__rep04 | 53.737 | 601/1000 | 0.4444 | 0.6010 | 2.14e-23 | 1.8828 | 0.8120 | 0.3027 |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p8__rep05 | 67.293 | 666/1000 | 0.4444 | 0.6660 | 2.94e-45 | 2.4925 | 0.8453 | 0.2457 |
| stage2_deathERK_ld1p0_bdR2p0_bdT1p6__rep01 | 59.596 | 662/1000 | 0.4444 | 0.6620 | 1.11e-43 | 2.4482 | 0.8214 | 0.1491 |
| stage2_deathERK_ld1p0_bdR2p0_bdT1p6__rep02 | 56.246 | 617/1000 | 0.4444 | 0.6170 | 4.71e-28 | 2.0137 | 0.8111 | 0.1478 |
| stage2_deathERK_ld1p0_bdR2p0_bdT1p6__rep03 | 54.139 | 619/1000 | 0.4444 | 0.6190 | 1.14e-28 | 2.0308 | 0.8117 | 0.1359 |
| stage2_deathERK_ld1p0_bdR2p0_bdT1p6__rep04 | 44.791 | 605/1000 | 0.4444 | 0.6050 | 1.62e-24 | 1.9146 | 0.7798 | 0.2028 |
| stage2_deathERK_ld1p0_bdR2p0_bdT1p6__rep05 | 55.514 | 665/1000 | 0.4444 | 0.6650 | 7.34e-45 | 2.4813 | 0.8118 | 0.1530 |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p4__rep01 | 57.726 | 655/1000 | 0.4444 | 0.6550 | 5.32e-41 | 2.3732 | 0.8142 | 0.1736 |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p4__rep02 | 52.925 | 616/1000 | 0.4444 | 0.6160 | 9.49e-28 | 2.0052 | 0.8069 | 0.1598 |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p4__rep03 | 49.820 | 614/1000 | 0.4444 | 0.6140 | 3.81e-27 | 1.9883 | 0.8014 | 0.1605 |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p4__rep04 | 42.963 | 597/1000 | 0.4444 | 0.5970 | 2.64e-22 | 1.8517 | 0.7721 | 0.2158 |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p4__rep05 | 52.085 | 663/1000 | 0.4444 | 0.6630 | 4.5e-44 | 2.4592 | 0.8089 | 0.1580 |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p8__rep01 | 46.168 | 657/1000 | 0.4444 | 0.6570 | 9.31e-42 | 2.3943 | 0.7788 | 0.0961 |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p8__rep02 | 48.529 | 628/1000 | 0.4444 | 0.6280 | 1.59e-31 | 2.1102 | 0.7807 | 0.0943 |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p8__rep03 | 42.263 | 631/1000 | 0.4444 | 0.6310 | 1.64e-32 | 2.1375 | 0.7656 | 0.0939 |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p8__rep04 | 35.048 | 623/1000 | 0.4444 | 0.6230 | 6.4e-30 | 2.0656 | 0.7336 | 0.1310 |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p8__rep05 | 46.652 | 678/1000 | 0.4444 | 0.6780 | 3.6e-50 | 2.6320 | 0.7897 | 0.0939 |
| stage2_deathERK_ld1p0_bdR4p0_bdT1p6__rep01 | 42.080 | 663/1000 | 0.4444 | 0.6630 | 4.5e-44 | 2.4592 | 0.7617 | 0.0534 |
| stage2_deathERK_ld1p0_bdR4p0_bdT1p6__rep02 | 46.988 | 635/1000 | 0.4444 | 0.6350 | 7.5e-34 | 2.1747 | 0.7676 | 0.0544 |
| stage2_deathERK_ld1p0_bdR4p0_bdT1p6__rep03 | 38.861 | 632/1000 | 0.4444 | 0.6320 | 7.64e-33 | 2.1467 | 0.7444 | 0.0504 |
| stage2_deathERK_ld1p0_bdR4p0_bdT1p6__rep04 | 31.914 | 629/1000 | 0.4444 | 0.6290 | 7.48e-32 | 2.1193 | 0.7200 | 0.0719 |
| stage2_deathERK_ld1p0_bdR4p0_bdT1p6__rep05 | 41.614 | 663/1000 | 0.4444 | 0.6630 | 4.5e-44 | 2.4592 | 0.7639 | 0.0500 |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p4__rep01 | 130.330 | 553/1000 | 0.4444 | 0.5530 | 3.72e-12 | 1.5464 | 0.9612 | 0.6619 |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p4__rep02 | 131.815 | 575/1000 | 0.4444 | 0.5750 | 8.33e-17 | 1.6912 | 0.9611 | 0.6561 |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p4__rep03 | 118.009 | 598/1000 | 0.4444 | 0.5980 | 1.42e-22 | 1.8595 | 0.9545 | 0.6380 |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p4__rep04 | 142.874 | 612/1000 | 0.4444 | 0.6120 | 1.51e-26 | 1.9716 | 0.9626 | 0.6553 |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p4__rep05 | 143.704 | 564/1000 | 0.4444 | 0.5640 | 2.24e-14 | 1.6170 | 0.9613 | 0.6704 |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p8__rep01 | 94.548 | 598/1000 | 0.4444 | 0.5980 | 1.42e-22 | 1.8595 | 0.9391 | 0.5423 |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p8__rep02 | 95.030 | 647/1000 | 0.4444 | 0.6470 | 4.74e-38 | 2.2911 | 0.9398 | 0.5384 |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p8__rep03 | 83.717 | 618/1000 | 0.4444 | 0.6180 | 2.32e-28 | 2.0223 | 0.9372 | 0.5615 |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p8__rep04 | 89.088 | 641/1000 | 0.4444 | 0.6410 | 6.45e-36 | 2.2319 | 0.9407 | 0.5577 |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p8__rep05 | 83.938 | 589/1000 | 0.4444 | 0.5890 | 3.32e-20 | 1.7914 | 0.9330 | 0.5588 |
| stage2_deathERK_ld2p0_bdR1p0_bdT1p6__rep01 | 59.719 | 608/1000 | 0.4444 | 0.6080 | 2.23e-25 | 1.9388 | 0.8995 | 0.4442 |
| stage2_deathERK_ld2p0_bdR1p0_bdT1p6__rep02 | 66.603 | 679/1000 | 0.4444 | 0.6790 | 1.36e-50 | 2.6441 | 0.9175 | 0.4217 |
| stage2_deathERK_ld2p0_bdR1p0_bdT1p6__rep03 | 57.150 | 643/1000 | 0.4444 | 0.6430 | 1.28e-36 | 2.2514 | 0.9061 | 0.4380 |
| stage2_deathERK_ld2p0_bdR1p0_bdT1p6__rep04 | 58.825 | 615/1000 | 0.4444 | 0.6150 | 1.91e-27 | 1.9968 | 0.9173 | 0.4653 |
| stage2_deathERK_ld2p0_bdR1p0_bdT1p6__rep05 | 57.553 | 563/1000 | 0.4444 | 0.5630 | 3.64e-14 | 1.6104 | 0.9135 | 0.4642 |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p4__rep01 | 54.974 | 578/1000 | 0.4444 | 0.5780 | 1.66e-17 | 1.7121 | 0.9028 | 0.4908 |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p4__rep02 | 72.006 | 638/1000 | 0.4444 | 0.6380 | 7.09e-35 | 2.2030 | 0.9188 | 0.4547 |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p4__rep03 | 60.947 | 611/1000 | 0.4444 | 0.6110 | 2.97e-26 | 1.9634 | 0.9146 | 0.4740 |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p4__rep04 | 57.949 | 570/1000 | 0.4444 | 0.5700 | 1.13e-15 | 1.6570 | 0.9078 | 0.4932 |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p4__rep05 | 55.280 | 561/1000 | 0.4444 | 0.5610 | 9.49e-14 | 1.5974 | 0.9035 | 0.4696 |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p8__rep01 | 40.488 | 620/1000 | 0.4444 | 0.6200 | 5.59e-29 | 2.0395 | 0.8710 | 0.3515 |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p8__rep02 | 46.744 | 655/1000 | 0.4444 | 0.6550 | 5.32e-41 | 2.3732 | 0.8788 | 0.3326 |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p8__rep03 | 43.453 | 657/1000 | 0.4444 | 0.6570 | 9.31e-42 | 2.3943 | 0.8636 | 0.3421 |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p8__rep04 | 40.066 | 564/1000 | 0.4444 | 0.5640 | 2.24e-14 | 1.6170 | 0.8733 | 0.3570 |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p8__rep05 | 41.997 | 576/1000 | 0.4444 | 0.5760 | 4.89e-17 | 1.6981 | 0.8757 | 0.3367 |
| stage2_deathERK_ld2p0_bdR2p0_bdT1p6__rep01 | 33.880 | 663/1000 | 0.4444 | 0.6630 | 4.5e-44 | 2.4592 | 0.8317 | 0.2341 |
| stage2_deathERK_ld2p0_bdR2p0_bdT1p6__rep02 | 34.089 | 688/1000 | 0.4444 | 0.6880 | 1.76e-54 | 2.7564 | 0.8537 | 0.2409 |
| stage2_deathERK_ld2p0_bdR2p0_bdT1p6__rep03 | 35.944 | 677/1000 | 0.4444 | 0.6770 | 9.48e-50 | 2.6200 | 0.8416 | 0.2381 |
| stage2_deathERK_ld2p0_bdR2p0_bdT1p6__rep04 | 28.594 | 551/1000 | 0.4444 | 0.5510 | 8.96e-12 | 1.5340 | 0.8411 | 0.2496 |
| stage2_deathERK_ld2p0_bdR2p0_bdT1p6__rep05 | 34.652 | 631/1000 | 0.4444 | 0.6310 | 1.64e-32 | 2.1375 | 0.8462 | 0.2297 |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p4__rep01 | 33.533 | 652/1000 | 0.4444 | 0.6520 | 7.02e-40 | 2.3420 | 0.8315 | 0.2574 |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p4__rep02 | 30.846 | 645/1000 | 0.4444 | 0.6450 | 2.48e-37 | 2.2711 | 0.8353 | 0.2655 |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p4__rep03 | 34.996 | 669/1000 | 0.4444 | 0.6690 | 1.85e-46 | 2.5264 | 0.8319 | 0.2561 |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p4__rep04 | 28.767 | 558/1000 | 0.4444 | 0.5580 | 3.87e-13 | 1.5781 | 0.8194 | 0.2633 |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p4__rep05 | 28.516 | 608/1000 | 0.4444 | 0.6080 | 2.23e-25 | 1.9388 | 0.8163 | 0.2359 |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p8__rep01 | 26.628 | 668/1000 | 0.4444 | 0.6680 | 4.67e-46 | 2.5151 | 0.7908 | 0.1633 |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p8__rep02 | 20.897 | 623/1000 | 0.4444 | 0.6230 | 6.4e-30 | 2.0656 | 0.7912 | 0.1819 |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p8__rep03 | 27.568 | 626/1000 | 0.4444 | 0.6260 | 7.06e-31 | 2.0922 | 0.7929 | 0.1691 |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p8__rep04 | 25.244 | 595/1000 | 0.4444 | 0.5950 | 9.06e-22 | 1.8364 | 0.8075 | 0.1633 |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p8__rep05 | 24.278 | 603/1000 | 0.4444 | 0.6030 | 5.93e-24 | 1.8986 | 0.8116 | 0.1513 |
| stage2_deathERK_ld2p0_bdR4p0_bdT1p6__rep01 | 21.938 | 654/1000 | 0.4444 | 0.6540 | 1.26e-40 | 2.3627 | 0.7626 | 0.1049 |
| stage2_deathERK_ld2p0_bdR4p0_bdT1p6__rep02 | 17.472 | 647/1000 | 0.4444 | 0.6470 | 4.74e-38 | 2.2911 | 0.7629 | 0.1039 |
| stage2_deathERK_ld2p0_bdR4p0_bdT1p6__rep03 | 23.432 | 644/1000 | 0.4444 | 0.6440 | 5.64e-37 | 2.2612 | 0.7598 | 0.1109 |
| stage2_deathERK_ld2p0_bdR4p0_bdT1p6__rep04 | 22.088 | 624/1000 | 0.4444 | 0.6240 | 3.08e-30 | 2.0745 | 0.7752 | 0.0944 |
| stage2_deathERK_ld2p0_bdR4p0_bdT1p6__rep05 | 20.830 | 587/1000 | 0.4444 | 0.5870 | 1.07e-19 | 1.7766 | 0.7837 | 0.0845 |

## Aggregated replicate summary

| parameter setting | replicates | mean T-zone fraction | mean T-zone ratio | mean binomial p-value | fraction p<0.05 | mean time | mean death rejection |
|---|---:|---:|---:|---:|---:|---:|---:|
| stage2_deathERK_ld0p5_bdR1p0_bdT0p4 | 5 | 0.6010 | 1.8856 | 3.67e-19 | 1.000 | 205.669 | 0.9024 |
| stage2_deathERK_ld0p5_bdR1p0_bdT0p8 | 5 | 0.6138 | 1.9963 | 6.41e-19 | 1.000 | 150.621 | 0.8645 |
| stage2_deathERK_ld0p5_bdR1p0_bdT1p6 | 5 | 0.6402 | 2.2297 | 7.58e-31 | 1.000 | 117.016 | 0.8276 |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p4 | 5 | 0.6296 | 2.1306 | 8.68e-26 | 1.000 | 116.455 | 0.8240 |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p8 | 5 | 0.6294 | 2.1326 | 1.21e-24 | 1.000 | 95.005 | 0.7897 |
| stage2_deathERK_ld0p5_bdR2p0_bdT1p6 | 5 | 0.6320 | 2.1678 | 3.8e-20 | 1.000 | 86.061 | 0.7686 |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p4 | 5 | 0.6236 | 2.0870 | 3.35e-22 | 1.000 | 83.424 | 0.7604 |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p8 | 5 | 0.6330 | 2.1843 | 2.9e-23 | 1.000 | 77.990 | 0.7457 |
| stage2_deathERK_ld0p5_bdR4p0_bdT1p6 | 5 | 0.6362 | 2.2223 | 4.44e-24 | 1.000 | 73.414 | 0.7306 |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p4 | 5 | 0.5858 | 1.7799 | 1.19e-14 | 1.000 | 159.778 | 0.9371 |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p8 | 5 | 0.6044 | 1.9247 | 1.38e-19 | 1.000 | 112.005 | 0.9085 |
| stage2_deathERK_ld1p0_bdR1p0_bdT1p6 | 5 | 0.6254 | 2.0952 | 1.52e-23 | 1.000 | 80.175 | 0.8728 |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p4 | 5 | 0.6040 | 1.9084 | 7.93e-22 | 1.000 | 82.725 | 0.8751 |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p8 | 5 | 0.6286 | 2.1284 | 4.27e-24 | 1.000 | 64.144 | 0.8366 |
| stage2_deathERK_ld1p0_bdR2p0_bdT1p6 | 5 | 0.6336 | 2.1777 | 3.23e-25 | 1.000 | 54.057 | 0.8072 |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p4 | 5 | 0.6290 | 2.1355 | 5.28e-23 | 1.000 | 51.104 | 0.8007 |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p8 | 5 | 0.6434 | 2.2679 | 1.32e-30 | 1.000 | 43.732 | 0.7697 |
| stage2_deathERK_ld1p0_bdR4p0_bdT1p6 | 5 | 0.6444 | 2.2718 | 1.66e-32 | 1.000 | 40.291 | 0.7515 |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p4 | 5 | 0.5804 | 1.7371 | 7.49e-13 | 1.000 | 133.347 | 0.9601 |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p8 | 5 | 0.6186 | 2.0392 | 6.66e-21 | 1.000 | 89.264 | 0.9380 |
| stage2_deathERK_ld2p0_bdR1p0_bdT1p6 | 5 | 0.6216 | 2.0883 | 7.29e-15 | 1.000 | 59.970 | 0.9108 |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p4 | 5 | 0.5916 | 1.8266 | 1.92e-14 | 1.000 | 60.231 | 0.9095 |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p8 | 5 | 0.6144 | 2.0244 | 4.5e-15 | 1.000 | 42.550 | 0.8725 |
| stage2_deathERK_ld2p0_bdR2p0_bdT1p6 | 5 | 0.6420 | 2.3014 | 1.79e-12 | 1.000 | 33.432 | 0.8429 |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p4 | 5 | 0.6264 | 2.1313 | 7.74e-14 | 1.000 | 31.332 | 0.8269 |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p8 | 5 | 0.6230 | 2.0816 | 1.82e-22 | 1.000 | 24.923 | 0.7988 |
| stage2_deathERK_ld2p0_bdR4p0_bdT1p6 | 5 | 0.6312 | 2.1532 | 2.13e-20 | 1.000 | 21.152 | 0.7688 |

## Pooled Bonferroni summary

| parameter setting | pooled k/n | pooled fraction | pooled ratio | pooled p-value | Bonferroni alpha | significant |
|---|---:|---:|---:|---:|---:|---:|
| stage2_deathERK_ld1p0_bdR4p0_bdT1p6 | 3222/5000 | 0.6444 | 2.2652 | 1.02e-177 | 0.00185185 | True |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p8 | 3217/5000 | 0.6434 | 2.2553 | 6.03e-176 | 0.00185185 | True |
| stage2_deathERK_ld2p0_bdR2p0_bdT1p6 | 3210/5000 | 0.6420 | 2.2416 | 1.76e-173 | 0.00185185 | True |
| stage2_deathERK_ld0p5_bdR1p0_bdT1p6 | 3201/5000 | 0.6402 | 2.2242 | 2.44e-170 | 0.00185185 | True |
| stage2_deathERK_ld0p5_bdR4p0_bdT1p6 | 3181/5000 | 0.6362 | 2.1860 | 1.82e-163 | 0.00185185 | True |
| stage2_deathERK_ld1p0_bdR2p0_bdT1p6 | 3168/5000 | 0.6336 | 2.1616 | 4.44e-159 | 0.00185185 | True |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p8 | 3165/5000 | 0.6330 | 2.1560 | 4.48e-158 | 0.00185185 | True |
| stage2_deathERK_ld0p5_bdR2p0_bdT1p6 | 3160/5000 | 0.6320 | 2.1467 | 2.07e-156 | 0.00185185 | True |
| stage2_deathERK_ld2p0_bdR4p0_bdT1p6 | 3156/5000 | 0.6312 | 2.1394 | 4.38e-155 | 0.00185185 | True |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p4 | 3148/5000 | 0.6296 | 2.1247 | 1.88e-152 | 0.00185185 | True |
| stage2_deathERK_ld0p5_bdR2p0_bdT0p8 | 3147/5000 | 0.6294 | 2.1229 | 3.99e-152 | 0.00185185 | True |
| stage2_deathERK_ld1p0_bdR4p0_bdT0p4 | 3145/5000 | 0.6290 | 2.1193 | 1.8e-151 | 0.00185185 | True |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p8 | 3143/5000 | 0.6286 | 2.1156 | 8.07e-151 | 0.00185185 | True |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p4 | 3132/5000 | 0.6264 | 2.0958 | 2.93e-147 | 0.00185185 | True |
| stage2_deathERK_ld1p0_bdR1p0_bdT1p6 | 3127/5000 | 0.6254 | 2.0869 | 1.18e-145 | 0.00185185 | True |
| stage2_deathERK_ld0p5_bdR4p0_bdT0p4 | 3118/5000 | 0.6236 | 2.0709 | 8.6e-143 | 0.00185185 | True |
| stage2_deathERK_ld2p0_bdR4p0_bdT0p8 | 3115/5000 | 0.6230 | 2.0656 | 7.62e-142 | 0.00185185 | True |
| stage2_deathERK_ld2p0_bdR1p0_bdT1p6 | 3108/5000 | 0.6216 | 2.0534 | 1.2e-139 | 0.00185185 | True |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p8 | 3093/5000 | 0.6186 | 2.0274 | 5.38e-135 | 0.00185185 | True |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p8 | 3072/5000 | 0.6144 | 1.9917 | 1.26e-128 | 0.00185185 | True |
| stage2_deathERK_ld0p5_bdR1p0_bdT0p8 | 3069/5000 | 0.6138 | 1.9867 | 9.98e-128 | 0.00185185 | True |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p8 | 3022/5000 | 0.6044 | 1.9098 | 4.2e-114 | 0.00185185 | True |
| stage2_deathERK_ld1p0_bdR2p0_bdT0p4 | 3020/5000 | 0.6040 | 1.9066 | 1.53e-113 | 0.00185185 | True |
| stage2_deathERK_ld0p5_bdR1p0_bdT0p4 | 3005/5000 | 0.6010 | 1.8828 | 2.26e-109 | 0.00185185 | True |
| stage2_deathERK_ld2p0_bdR2p0_bdT0p4 | 2958/5000 | 0.5916 | 1.8107 | 7.75e-97 | 0.00185185 | True |
| stage2_deathERK_ld1p0_bdR1p0_bdT0p4 | 2929/5000 | 0.5858 | 1.7679 | 1.69e-89 | 0.00185185 | True |
| stage2_deathERK_ld2p0_bdR1p0_bdT0p4 | 2902/5000 | 0.5804 | 1.7290 | 6.18e-83 | 0.00185185 | True |

## Short conclusion

The first setting `stage2_deathERK_ld0p5_bdR1p0_bdT0p4__rep01` has observed T-zone fraction 0.5970, area fraction 0.4444, one-sided binomial p-value 2.64e-22, and density ratio 1.8517.
The largest replicate-level T-zone density ratio occurs in `stage2_deathERK_ld0p5_bdR4p0_bdT1p6__rep02` with ratio 2.8084.
The smallest replicate-level one-sided binomial p-value occurs in `stage2_deathERK_ld0p5_bdR4p0_bdT1p6__rep02` with p-value 2.92e-56.
