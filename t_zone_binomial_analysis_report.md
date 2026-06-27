# T-zone binomial enrichment analysis after burn-in
All 2D histograms use the same colorbar scale: 0 to 14 observed deaths per bin.
Burn-in removes the initial transient period: the first 500 observed deaths are simulated but not used for statistics. The reported statistics use the next 1000 observed deaths.
## Interpretation rules
- `t_zone_density_ratio`: values > 1 mean deaths are denser in the T-zone than outside; values < 1 mean the opposite.
- `t_zone_binomial_p_value_greater` is the one-sided exact binomial p-value for `H0: pi = |T|/|W|` against `HA: pi > |T|/|W|`.
- A small p-value, for example < 0.05, means that the observed number of deaths in the T-zone is unusually large under a spatially uniform death-location model.

## Results
| scenario | time for analyzed 1000 deaths | T-zone deaths k/n | area fraction p0 | observed fraction | binomial p-value greater | T-zone density ratio |
|---|---:|---:|---:|---:|---:|---:|
| deathERK_ld0p5_bdR1p0_bdT0p4__rep01 | 226.286 | 590/1000 | 0.444 | 0.590 | 1.84e-20 | 1.799 |
| deathERK_ld0p5_bdR1p0_bdT0p4__rep02 | 194.848 | 598/1000 | 0.444 | 0.598 | 1.42e-22 | 1.859 |
| deathERK_ld0p5_bdR1p0_bdT0p4__rep03 | 216.604 | 598/1000 | 0.444 | 0.598 | 1.42e-22 | 1.859 |
| deathERK_ld0p5_bdR1p0_bdT0p4__rep04 | 225.870 | 608/1000 | 0.444 | 0.608 | 2.23e-25 | 1.939 |
| deathERK_ld0p5_bdR1p0_bdT0p4__rep05 | 210.579 | 590/1000 | 0.444 | 0.590 | 1.84e-20 | 1.799 |
| deathERK_ld0p5_bdR1p0_bdT0p8__rep01 | 157.026 | 604/1000 | 0.444 | 0.604 | 3.1e-24 | 1.907 |
| deathERK_ld0p5_bdR1p0_bdT0p8__rep02 | 180.519 | 684/1000 | 0.444 | 0.684 | 9.87e-53 | 2.706 |
| deathERK_ld0p5_bdR1p0_bdT0p8__rep03 | 163.462 | 641/1000 | 0.444 | 0.641 | 6.45e-36 | 2.232 |
| deathERK_ld0p5_bdR1p0_bdT0p8__rep04 | 147.334 | 600/1000 | 0.444 | 0.600 | 4.03e-23 | 1.875 |
| deathERK_ld0p5_bdR1p0_bdT0p8__rep05 | 146.887 | 617/1000 | 0.444 | 0.617 | 4.71e-28 | 2.014 |
| deathERK_ld0p5_bdR1p0_bdT1p6__rep01 | 120.139 | 656/1000 | 0.444 | 0.656 | 2.23e-41 | 2.384 |
| deathERK_ld0p5_bdR1p0_bdT1p6__rep02 | 121.962 | 686/1000 | 0.444 | 0.686 | 1.33e-53 | 2.731 |
| deathERK_ld0p5_bdR1p0_bdT1p6__rep03 | 112.919 | 639/1000 | 0.444 | 0.639 | 3.2e-35 | 2.213 |
| deathERK_ld0p5_bdR1p0_bdT1p6__rep04 | 119.250 | 649/1000 | 0.444 | 0.649 | 8.91e-39 | 2.311 |
| deathERK_ld0p5_bdR1p0_bdT1p6__rep05 | 105.909 | 594/1000 | 0.444 | 0.594 | 1.67e-21 | 1.829 |
| deathERK_ld0p5_bdR2p0_bdT0p4__rep01 | 118.708 | 649/1000 | 0.444 | 0.649 | 8.91e-39 | 2.311 |
| deathERK_ld0p5_bdR2p0_bdT0p4__rep02 | 126.737 | 659/1000 | 0.444 | 0.659 | 1.6e-42 | 2.416 |
| deathERK_ld0p5_bdR2p0_bdT0p4__rep03 | 125.692 | 647/1000 | 0.444 | 0.647 | 4.74e-38 | 2.291 |
| deathERK_ld0p5_bdR2p0_bdT0p4__rep04 | 108.421 | 607/1000 | 0.444 | 0.607 | 4.34e-25 | 1.931 |
| deathERK_ld0p5_bdR2p0_bdT0p4__rep05 | 119.108 | 608/1000 | 0.444 | 0.608 | 2.23e-25 | 1.939 |
| deathERK_ld0p5_bdR2p0_bdT0p8__rep01 | 93.133 | 617/1000 | 0.444 | 0.617 | 4.71e-28 | 2.014 |
| deathERK_ld0p5_bdR2p0_bdT0p8__rep02 | 108.793 | 661/1000 | 0.444 | 0.661 | 2.71e-43 | 2.437 |
| deathERK_ld0p5_bdR2p0_bdT0p8__rep03 | 102.612 | 628/1000 | 0.444 | 0.628 | 1.59e-31 | 2.110 |
| deathERK_ld0p5_bdR2p0_bdT0p8__rep04 | 105.239 | 637/1000 | 0.444 | 0.637 | 1.56e-34 | 2.194 |
| deathERK_ld0p5_bdR2p0_bdT0p8__rep05 | 88.440 | 579/1000 | 0.444 | 0.579 | 9.65e-18 | 1.719 |
| deathERK_ld0p5_bdR2p0_bdT1p6__rep01 | 84.105 | 655/1000 | 0.444 | 0.655 | 5.32e-41 | 2.373 |
| deathERK_ld0p5_bdR2p0_bdT1p6__rep02 | 80.957 | 626/1000 | 0.444 | 0.626 | 7.06e-31 | 2.092 |
| deathERK_ld0p5_bdR2p0_bdT1p6__rep03 | 96.119 | 680/1000 | 0.444 | 0.680 | 5.14e-51 | 2.656 |
| deathERK_ld0p5_bdR2p0_bdT1p6__rep04 | 84.501 | 656/1000 | 0.444 | 0.656 | 2.23e-41 | 2.384 |
| deathERK_ld0p5_bdR2p0_bdT1p6__rep05 | 81.592 | 615/1000 | 0.444 | 0.615 | 1.91e-27 | 1.997 |
| deathERK_ld0p5_bdR4p0_bdT0p4__rep01 | 89.123 | 696/1000 | 0.444 | 0.696 | 4.49e-58 | 2.862 |
| deathERK_ld0p5_bdR4p0_bdT0p4__rep02 | 86.982 | 627/1000 | 0.444 | 0.627 | 3.35e-31 | 2.101 |
| deathERK_ld0p5_bdR4p0_bdT0p4__rep03 | 88.231 | 617/1000 | 0.444 | 0.617 | 4.71e-28 | 2.014 |
| deathERK_ld0p5_bdR4p0_bdT0p4__rep04 | 82.642 | 640/1000 | 0.444 | 0.640 | 1.44e-35 | 2.222 |
| deathERK_ld0p5_bdR4p0_bdT0p4__rep05 | 74.768 | 624/1000 | 0.444 | 0.624 | 3.08e-30 | 2.074 |
| deathERK_ld0p5_bdR4p0_bdT0p8__rep01 | 79.037 | 667/1000 | 0.444 | 0.667 | 1.18e-45 | 2.504 |
| deathERK_ld0p5_bdR4p0_bdT0p8__rep02 | 71.282 | 623/1000 | 0.444 | 0.623 | 6.4e-30 | 2.066 |
| deathERK_ld0p5_bdR4p0_bdT0p8__rep03 | 74.194 | 666/1000 | 0.444 | 0.666 | 2.94e-45 | 2.493 |
| deathERK_ld0p5_bdR4p0_bdT0p8__rep04 | 79.789 | 656/1000 | 0.444 | 0.656 | 2.23e-41 | 2.384 |
| deathERK_ld0p5_bdR4p0_bdT0p8__rep05 | 72.424 | 650/1000 | 0.444 | 0.650 | 3.84e-39 | 2.321 |
| deathERK_ld0p5_bdR4p0_bdT1p6__rep01 | 71.179 | 642/1000 | 0.444 | 0.642 | 2.87e-36 | 2.242 |
| deathERK_ld0p5_bdR4p0_bdT1p6__rep02 | 72.428 | 683/1000 | 0.444 | 0.683 | 2.67e-52 | 2.693 |
| deathERK_ld0p5_bdR4p0_bdT1p6__rep03 | 82.179 | 662/1000 | 0.444 | 0.662 | 1.11e-43 | 2.448 |
| deathERK_ld0p5_bdR4p0_bdT1p6__rep04 | 84.780 | 691/1000 | 0.444 | 0.691 | 8.2e-56 | 2.795 |
| deathERK_ld0p5_bdR4p0_bdT1p6__rep05 | 70.030 | 643/1000 | 0.444 | 0.643 | 1.28e-36 | 2.251 |
| deathERK_ld1p0_bdR1p0_bdT0p4__rep01 | 156.461 | 610/1000 | 0.444 | 0.610 | 5.85e-26 | 1.955 |
| deathERK_ld1p0_bdR1p0_bdT0p4__rep02 | 192.081 | 585/1000 | 0.444 | 0.585 | 3.37e-19 | 1.762 |
| deathERK_ld1p0_bdR1p0_bdT0p4__rep03 | 169.649 | 598/1000 | 0.444 | 0.598 | 1.42e-22 | 1.859 |
| deathERK_ld1p0_bdR1p0_bdT0p4__rep04 | 163.079 | 578/1000 | 0.444 | 0.578 | 1.66e-17 | 1.712 |
| deathERK_ld1p0_bdR1p0_bdT0p4__rep05 | 162.898 | 612/1000 | 0.444 | 0.612 | 1.51e-26 | 1.972 |
| deathERK_ld1p0_bdR1p0_bdT0p8__rep01 | 109.183 | 608/1000 | 0.444 | 0.608 | 2.23e-25 | 1.939 |
| deathERK_ld1p0_bdR1p0_bdT0p8__rep02 | 124.545 | 631/1000 | 0.444 | 0.631 | 1.64e-32 | 2.138 |
| deathERK_ld1p0_bdR1p0_bdT0p8__rep03 | 116.225 | 644/1000 | 0.444 | 0.644 | 5.64e-37 | 2.261 |
| deathERK_ld1p0_bdR1p0_bdT0p8__rep04 | 107.027 | 610/1000 | 0.444 | 0.610 | 5.85e-26 | 1.955 |
| deathERK_ld1p0_bdR1p0_bdT0p8__rep05 | 106.145 | 635/1000 | 0.444 | 0.635 | 7.5e-34 | 2.175 |
| deathERK_ld1p0_bdR1p0_bdT1p6__rep01 | 67.946 | 600/1000 | 0.444 | 0.600 | 4.03e-23 | 1.875 |
| deathERK_ld1p0_bdR1p0_bdT1p6__rep02 | 71.462 | 608/1000 | 0.444 | 0.608 | 2.23e-25 | 1.939 |
| deathERK_ld1p0_bdR1p0_bdT1p6__rep03 | 88.663 | 638/1000 | 0.444 | 0.638 | 7.09e-35 | 2.203 |
| deathERK_ld1p0_bdR1p0_bdT1p6__rep04 | 84.545 | 612/1000 | 0.444 | 0.612 | 1.51e-26 | 1.972 |
| deathERK_ld1p0_bdR1p0_bdT1p6__rep05 | 71.792 | 631/1000 | 0.444 | 0.631 | 1.64e-32 | 2.138 |
| deathERK_ld1p0_bdR2p0_bdT0p4__rep01 | 79.864 | 624/1000 | 0.444 | 0.624 | 3.08e-30 | 2.074 |
| deathERK_ld1p0_bdR2p0_bdT0p4__rep02 | 74.841 | 645/1000 | 0.444 | 0.645 | 2.48e-37 | 2.271 |
| deathERK_ld1p0_bdR2p0_bdT0p4__rep03 | 68.223 | 585/1000 | 0.444 | 0.585 | 3.37e-19 | 1.762 |
| deathERK_ld1p0_bdR2p0_bdT0p4__rep04 | 84.058 | 625/1000 | 0.444 | 0.625 | 1.48e-30 | 2.083 |
| deathERK_ld1p0_bdR2p0_bdT0p4__rep05 | 83.221 | 622/1000 | 0.444 | 0.622 | 1.32e-29 | 2.057 |
| deathERK_ld1p0_bdR2p0_bdT0p8__rep01 | 59.099 | 636/1000 | 0.444 | 0.636 | 3.43e-34 | 2.184 |
| deathERK_ld1p0_bdR2p0_bdT0p8__rep02 | 61.560 | 644/1000 | 0.444 | 0.644 | 5.64e-37 | 2.261 |
| deathERK_ld1p0_bdR2p0_bdT0p8__rep03 | 63.616 | 651/1000 | 0.444 | 0.651 | 1.64e-39 | 2.332 |
| deathERK_ld1p0_bdR2p0_bdT0p8__rep04 | 64.880 | 636/1000 | 0.444 | 0.636 | 3.43e-34 | 2.184 |
| deathERK_ld1p0_bdR2p0_bdT0p8__rep05 | 53.389 | 598/1000 | 0.444 | 0.598 | 1.42e-22 | 1.859 |
| deathERK_ld1p0_bdR2p0_bdT1p6__rep01 | 47.849 | 589/1000 | 0.444 | 0.589 | 3.32e-20 | 1.791 |
| deathERK_ld1p0_bdR2p0_bdT1p6__rep02 | 48.525 | 625/1000 | 0.444 | 0.625 | 1.48e-30 | 2.083 |
| deathERK_ld1p0_bdR2p0_bdT1p6__rep03 | 45.821 | 610/1000 | 0.444 | 0.610 | 5.85e-26 | 1.955 |
| deathERK_ld1p0_bdR2p0_bdT1p6__rep04 | 50.764 | 597/1000 | 0.444 | 0.597 | 2.64e-22 | 1.852 |
| deathERK_ld1p0_bdR2p0_bdT1p6__rep05 | 49.383 | 649/1000 | 0.444 | 0.649 | 8.91e-39 | 2.311 |
| deathERK_ld1p0_bdR4p0_bdT0p4__rep01 | 50.093 | 632/1000 | 0.444 | 0.632 | 7.64e-33 | 2.147 |
| deathERK_ld1p0_bdR4p0_bdT0p4__rep02 | 47.953 | 592/1000 | 0.444 | 0.592 | 5.58e-21 | 1.814 |
| deathERK_ld1p0_bdR4p0_bdT0p4__rep03 | 48.300 | 646/1000 | 0.444 | 0.646 | 1.09e-37 | 2.281 |
| deathERK_ld1p0_bdR4p0_bdT0p4__rep04 | 50.525 | 625/1000 | 0.444 | 0.625 | 1.48e-30 | 2.083 |
| deathERK_ld1p0_bdR4p0_bdT0p4__rep05 | 44.479 | 624/1000 | 0.444 | 0.624 | 3.08e-30 | 2.074 |
| deathERK_ld1p0_bdR4p0_bdT0p8__rep01 | 44.700 | 677/1000 | 0.444 | 0.677 | 9.48e-50 | 2.620 |
| deathERK_ld1p0_bdR4p0_bdT0p8__rep02 | 51.112 | 650/1000 | 0.444 | 0.650 | 3.84e-39 | 2.321 |
| deathERK_ld1p0_bdR4p0_bdT0p8__rep03 | 46.492 | 651/1000 | 0.444 | 0.651 | 1.64e-39 | 2.332 |
| deathERK_ld1p0_bdR4p0_bdT0p8__rep04 | 42.940 | 584/1000 | 0.444 | 0.584 | 5.95e-19 | 1.755 |
| deathERK_ld1p0_bdR4p0_bdT0p8__rep05 | 38.265 | 596/1000 | 0.444 | 0.596 | 4.9e-22 | 1.844 |
| deathERK_ld1p0_bdR4p0_bdT1p6__rep01 | 41.912 | 649/1000 | 0.444 | 0.649 | 8.91e-39 | 2.311 |
| deathERK_ld1p0_bdR4p0_bdT1p6__rep02 | 35.149 | 630/1000 | 0.444 | 0.630 | 3.51e-32 | 2.128 |
| deathERK_ld1p0_bdR4p0_bdT1p6__rep03 | 42.812 | 678/1000 | 0.444 | 0.678 | 3.6e-50 | 2.632 |
| deathERK_ld1p0_bdR4p0_bdT1p6__rep04 | 34.770 | 632/1000 | 0.444 | 0.632 | 7.64e-33 | 2.147 |
| deathERK_ld1p0_bdR4p0_bdT1p6__rep05 | 37.809 | 642/1000 | 0.444 | 0.642 | 2.87e-36 | 2.242 |
| deathERK_ld2p0_bdR1p0_bdT0p4__rep01 | 146.410 | 567/1000 | 0.444 | 0.567 | 5.12e-15 | 1.637 |
| deathERK_ld2p0_bdR1p0_bdT0p4__rep02 | 132.161 | 540/1000 | 0.444 | 0.540 | 8.43e-10 | 1.467 |
| deathERK_ld2p0_bdR1p0_bdT0p4__rep03 | 146.268 | 593/1000 | 0.444 | 0.593 | 3.06e-21 | 1.821 |
| deathERK_ld2p0_bdR1p0_bdT0p4__rep04 | 140.491 | 606/1000 | 0.444 | 0.606 | 8.39e-25 | 1.923 |
| deathERK_ld2p0_bdR1p0_bdT0p4__rep05 | 127.909 | 565/1000 | 0.444 | 0.565 | 1.38e-14 | 1.624 |
| deathERK_ld2p0_bdR1p0_bdT0p8__rep01 | 94.021 | 600/1000 | 0.444 | 0.600 | 4.03e-23 | 1.875 |
| deathERK_ld2p0_bdR1p0_bdT0p8__rep02 | 81.370 | 610/1000 | 0.444 | 0.610 | 5.85e-26 | 1.955 |
| deathERK_ld2p0_bdR1p0_bdT0p8__rep03 | 85.379 | 615/1000 | 0.444 | 0.615 | 1.91e-27 | 1.997 |
| deathERK_ld2p0_bdR1p0_bdT0p8__rep04 | 75.833 | 590/1000 | 0.444 | 0.590 | 1.84e-20 | 1.799 |
| deathERK_ld2p0_bdR1p0_bdT0p8__rep05 | 78.721 | 601/1000 | 0.444 | 0.601 | 2.14e-23 | 1.883 |
| deathERK_ld2p0_bdR1p0_bdT1p6__rep01 | 55.217 | 599/1000 | 0.444 | 0.599 | 7.57e-23 | 1.867 |
| deathERK_ld2p0_bdR1p0_bdT1p6__rep02 | 53.544 | 621/1000 | 0.444 | 0.621 | 2.73e-29 | 2.048 |
| deathERK_ld2p0_bdR1p0_bdT1p6__rep03 | 63.428 | 630/1000 | 0.444 | 0.630 | 3.51e-32 | 2.128 |
| deathERK_ld2p0_bdR1p0_bdT1p6__rep04 | 61.250 | 632/1000 | 0.444 | 0.632 | 7.64e-33 | 2.147 |
| deathERK_ld2p0_bdR1p0_bdT1p6__rep05 | 52.836 | 549/1000 | 0.444 | 0.549 | 2.12e-11 | 1.522 |
| deathERK_ld2p0_bdR2p0_bdT0p4__rep01 | 51.906 | 555/1000 | 0.444 | 0.555 | 1.52e-12 | 1.559 |
| deathERK_ld2p0_bdR2p0_bdT0p4__rep02 | 58.349 | 603/1000 | 0.444 | 0.603 | 5.93e-24 | 1.899 |
| deathERK_ld2p0_bdR2p0_bdT0p4__rep03 | 65.568 | 626/1000 | 0.444 | 0.626 | 7.06e-31 | 2.092 |
| deathERK_ld2p0_bdR2p0_bdT0p4__rep04 | 55.436 | 604/1000 | 0.444 | 0.604 | 3.1e-24 | 1.907 |
| deathERK_ld2p0_bdR2p0_bdT0p4__rep05 | 55.120 | 591/1000 | 0.444 | 0.591 | 1.02e-20 | 1.806 |
| deathERK_ld2p0_bdR2p0_bdT0p8__rep01 | 33.293 | 558/1000 | 0.444 | 0.558 | 3.87e-13 | 1.578 |
| deathERK_ld2p0_bdR2p0_bdT0p8__rep02 | 35.428 | 555/1000 | 0.444 | 0.555 | 1.52e-12 | 1.559 |
| deathERK_ld2p0_bdR2p0_bdT0p8__rep03 | 42.368 | 631/1000 | 0.444 | 0.631 | 1.64e-32 | 2.138 |
| deathERK_ld2p0_bdR2p0_bdT0p8__rep04 | 46.688 | 656/1000 | 0.444 | 0.656 | 2.23e-41 | 2.384 |
| deathERK_ld2p0_bdR2p0_bdT0p8__rep05 | 39.989 | 617/1000 | 0.444 | 0.617 | 4.71e-28 | 2.014 |
| deathERK_ld2p0_bdR2p0_bdT1p6__rep01 | 29.039 | 573/1000 | 0.444 | 0.573 | 2.39e-16 | 1.677 |
| deathERK_ld2p0_bdR2p0_bdT1p6__rep02 | 29.787 | 640/1000 | 0.444 | 0.640 | 1.44e-35 | 2.222 |
| deathERK_ld2p0_bdR2p0_bdT1p6__rep03 | 29.969 | 649/1000 | 0.444 | 0.649 | 8.91e-39 | 2.311 |
| deathERK_ld2p0_bdR2p0_bdT1p6__rep04 | 37.334 | 641/1000 | 0.444 | 0.641 | 6.45e-36 | 2.232 |
| deathERK_ld2p0_bdR2p0_bdT1p6__rep05 | 32.177 | 625/1000 | 0.444 | 0.625 | 1.48e-30 | 2.083 |
| deathERK_ld2p0_bdR4p0_bdT0p4__rep01 | 31.626 | 626/1000 | 0.444 | 0.626 | 7.06e-31 | 2.092 |
| deathERK_ld2p0_bdR4p0_bdT0p4__rep02 | 26.075 | 583/1000 | 0.444 | 0.583 | 1.05e-18 | 1.748 |
| deathERK_ld2p0_bdR4p0_bdT0p4__rep03 | 25.681 | 586/1000 | 0.444 | 0.586 | 1.9e-19 | 1.769 |
| deathERK_ld2p0_bdR4p0_bdT0p4__rep04 | 29.636 | 631/1000 | 0.444 | 0.631 | 1.64e-32 | 2.138 |
| deathERK_ld2p0_bdR4p0_bdT0p4__rep05 | 32.788 | 673/1000 | 0.444 | 0.673 | 4.34e-48 | 2.573 |
| deathERK_ld2p0_bdR4p0_bdT0p8__rep01 | 18.612 | 575/1000 | 0.444 | 0.575 | 8.33e-17 | 1.691 |
| deathERK_ld2p0_bdR4p0_bdT0p8__rep02 | 20.733 | 645/1000 | 0.444 | 0.645 | 2.48e-37 | 2.271 |
| deathERK_ld2p0_bdR4p0_bdT0p8__rep03 | 23.008 | 622/1000 | 0.444 | 0.622 | 1.32e-29 | 2.057 |
| deathERK_ld2p0_bdR4p0_bdT0p8__rep04 | 23.236 | 662/1000 | 0.444 | 0.662 | 1.11e-43 | 2.448 |
| deathERK_ld2p0_bdR4p0_bdT0p8__rep05 | 27.969 | 653/1000 | 0.444 | 0.653 | 2.98e-40 | 2.352 |
| deathERK_ld2p0_bdR4p0_bdT1p6__rep01 | 19.623 | 653/1000 | 0.444 | 0.653 | 2.98e-40 | 2.352 |
| deathERK_ld2p0_bdR4p0_bdT1p6__rep02 | 19.493 | 601/1000 | 0.444 | 0.601 | 2.14e-23 | 1.883 |
| deathERK_ld2p0_bdR4p0_bdT1p6__rep03 | 23.702 | 635/1000 | 0.444 | 0.635 | 7.5e-34 | 2.175 |
| deathERK_ld2p0_bdR4p0_bdT1p6__rep04 | 17.240 | 549/1000 | 0.444 | 0.549 | 2.12e-11 | 1.522 |
| deathERK_ld2p0_bdR4p0_bdT1p6__rep05 | 22.603 | 628/1000 | 0.444 | 0.628 | 1.59e-31 | 2.110 |

## Short conclusion
The baseline case has T-zone observed fraction = 0.590, area fraction = 0.444, one-sided binomial p-value = 1.84e-20, and T-zone density ratio = 1.799.
The largest T-zone density ratio occurs in `deathERK_ld0p5_bdR4p0_bdT0p4__rep01` with ratio = 2.862 and p-value = 4.49e-58.
The smallest one-sided binomial p-value occurs in `deathERK_ld0p5_bdR4p0_bdT0p4__rep01` with p-value = 4.49e-58 and ratio = 2.862.

## Aggregated replicate summary
| parameter setting | replicates | mean T-zone ratio | mean binomial p-value | fraction p<0.05 |
|---|---:|---:|---:|---:|
| deathERK_ld0p5_bdR1p0_bdT0p4 | 5 | 1.851 | 7.41e-21 | 1.000 |
| deathERK_ld0p5_bdR1p0_bdT0p8 | 5 | 2.147 | 8.68e-24 | 1.000 |
| deathERK_ld0p5_bdR1p0_bdT1p6 | 5 | 2.293 | 3.34e-22 | 1.000 |
| deathERK_ld0p5_bdR2p0_bdT0p4 | 5 | 2.177 | 1.31e-25 | 1.000 |
| deathERK_ld0p5_bdR2p0_bdT0p8 | 5 | 2.095 | 1.93e-18 | 1.000 |
| deathERK_ld0p5_bdR2p0_bdT1p6 | 5 | 2.300 | 3.81e-28 | 1.000 |
| deathERK_ld0p5_bdR4p0_bdT0p4 | 5 | 2.255 | 9.48e-29 | 1.000 |
| deathERK_ld0p5_bdR4p0_bdT0p8 | 5 | 2.353 | 1.28e-30 | 1.000 |
| deathERK_ld0p5_bdR4p0_bdT1p6 | 5 | 2.486 | 8.3e-37 | 1.000 |
| deathERK_ld1p0_bdR1p0_bdT0p4 | 5 | 1.852 | 3.4e-18 | 1.000 |
| deathERK_ld1p0_bdR1p0_bdT0p8 | 5 | 2.093 | 5.64e-26 | 1.000 |
| deathERK_ld1p0_bdR1p0_bdT1p6 | 5 | 2.025 | 8.11e-24 | 1.000 |
| deathERK_ld1p0_bdR2p0_bdT0p4 | 5 | 2.050 | 6.73e-20 | 1.000 |
| deathERK_ld1p0_bdR2p0_bdT0p8 | 5 | 2.164 | 2.83e-23 | 1.000 |
| deathERK_ld1p0_bdR2p0_bdT1p6 | 5 | 1.999 | 6.68e-21 | 1.000 |
| deathERK_ld1p0_bdR4p0_bdT0p4 | 5 | 2.080 | 1.12e-21 | 1.000 |
| deathERK_ld1p0_bdR4p0_bdT0p8 | 5 | 2.174 | 1.19e-19 | 1.000 |
| deathERK_ld1p0_bdR4p0_bdT1p6 | 5 | 2.292 | 8.55e-33 | 1.000 |
| deathERK_ld2p0_bdR1p0_bdT0p4 | 5 | 1.694 | 1.69e-10 | 1.000 |
| deathERK_ld2p0_bdR1p0_bdT0p8 | 5 | 1.902 | 3.69e-21 | 1.000 |
| deathERK_ld2p0_bdR1p0_bdT1p6 | 5 | 1.942 | 4.24e-12 | 1.000 |
| deathERK_ld2p0_bdR2p0_bdT0p4 | 5 | 1.853 | 3.05e-13 | 1.000 |
| deathERK_ld2p0_bdR2p0_bdT0p8 | 5 | 1.934 | 3.82e-13 | 1.000 |
| deathERK_ld2p0_bdR2p0_bdT1p6 | 5 | 2.105 | 4.78e-17 | 1.000 |
| deathERK_ld2p0_bdR4p0_bdT0p4 | 5 | 2.064 | 2.47e-19 | 1.000 |
| deathERK_ld2p0_bdR4p0_bdT0p8 | 5 | 2.164 | 1.67e-17 | 1.000 |
| deathERK_ld2p0_bdR4p0_bdT1p6 | 5 | 2.008 | 4.24e-12 | 1.000 |
