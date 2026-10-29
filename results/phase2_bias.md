# Phase 2b: energy bias of exact dyadic rotations (Lanczos)

## Other error sources on a reflection patch (kagome 2x3 cells, 18 sites)

| test | eps / eta | ||A'-A|| or max dw | A'^2 = 1 err | E0' - E0 | first-order bound |
|---|---|---|---|---|---|
| kagome 2x3, gridsynth W (exact inverse pair) | 0.01 | 1.17e-02 | 1.8e-14 | 9.68e-05 | 4.22e-01 |
| kagome 2x3, gridsynth W (exact inverse pair) | 0.001 | 1.13e-03 | 2.7e-14 | 2.46e-06 | 4.05e-02 |
| kagome 2x3, gridsynth W (exact inverse pair) | 0.0001 | 4.40e-05 | 3.8e-14 | 2.32e-09 | 1.58e-03 |
| kagome 2x3, PREPARE weights 1 + dw, sum dw = 0 | 0.01 | 1.00e-02 | - | -1.28e-04 | 1.31e-01 |
| kagome 2x3, PREPARE weights 1 + dw, sum dw = 0 | 0.001 | 1.00e-03 | - | -1.53e-06 | 1.94e-02 |
| kagome 2x3, PREPARE weights 1 + dw, sum dw = 0 | 0.0001 | 1.00e-04 | - | -1.24e-08 | 1.19e-03 |

## square J1-J2 J2=0.5 4x4 (tri_0.5)

N = 16, triangles = 64, decomposition check ||hw A + shift - h|| = 3.2e-15, E0/N = -2.114481

| bits | max angle error | ||A'-A|| | A' Hermitian err | A'^2 = 1 err | E0' - E0 | bound sum hw ||A'-A|| | bound / actual |
|---|---|---|---|---|---|---|---|
| 6 | 4.9e-02 | 5.29e-02 | 1.6e-16 | 8.4e-15 | -2.89e-04 | 2.54e+00 | 8767.5 |
| 8 | 1.2e-02 | 3.77e-03 | 9.5e-17 | 7.9e-15 | -3.48e-06 | 1.81e-01 | 52107.9 |
| 10 | 3.1e-03 | 3.77e-03 | 9.5e-17 | 7.9e-15 | -3.48e-06 | 1.81e-01 | 52107.9 |
| 12 | 7.7e-04 | 7.07e-04 | 9.9e-17 | 7.4e-15 | -1.27e-07 | 3.39e-02 | 267080.8 |
| 14 | 1.9e-04 | 6.02e-05 | 1.5e-16 | 8.2e-15 | -9.30e-10 | 2.89e-03 | 3105972.1 |
| 16 | 4.8e-05 | 6.02e-05 | 1.5e-16 | 8.2e-15 | -9.30e-10 | 2.89e-03 | 3105924.6 |

## square J1-J2 J2=0.25 4x4 (tri_0.25)

N = 16, triangles = 64, decomposition check ||hw A + shift - h|| = 3.3e-15, E0/N = -2.403818

| bits | max angle error | ||A'-A|| | A' Hermitian err | A'^2 = 1 err | E0' - E0 | bound sum hw ||A'-A|| | bound / actual |
|---|---|---|---|---|---|---|---|
| 6 | 4.9e-02 | 5.29e-02 | 1.8e-16 | - | -1.32e-01 | 2.54e+00 | 19.2 |
| 8 | 1.2e-02 | 3.77e-03 | 1.8e-16 | - | -1.17e-02 | 1.81e-01 | 15.5 |
| 10 | 3.1e-03 | 3.77e-03 | 1.8e-16 | - | -2.64e-03 | 1.81e-01 | 68.7 |
| 12 | 7.7e-04 | 7.07e-04 | 1.8e-16 | - | 1.58e-04 | 3.39e-02 | 214.4 |
| 14 | 1.9e-04 | 6.02e-05 | 1.4e-16 | - | 1.77e-04 | 2.89e-03 | 16.3 |
| 16 | 4.8e-05 | 6.02e-05 | 1.4e-16 | - | 3.59e-05 | 2.89e-03 | 80.4 |

## kagome 2x2 cells (tri_0.5 x 4)

N = 12, triangles = 8, decomposition check ||hw A + shift - h|| = 1.3e-14, E0/N = -1.814958

| bits | max angle error | ||A'-A|| | A' Hermitian err | A'^2 = 1 err | E0' - E0 | bound sum hw ||A'-A|| | bound / actual |
|---|---|---|---|---|---|---|---|
| 6 | 4.9e-02 | 5.29e-02 | 1.6e-16 | 8.4e-15 | 6.53e-03 | 1.27e+00 | 194.4 |
| 8 | 1.2e-02 | 3.77e-03 | 9.5e-17 | 7.9e-15 | 3.27e-05 | 9.06e-02 | 2771.9 |
| 10 | 3.1e-03 | 3.77e-03 | 9.5e-17 | 7.9e-15 | 3.27e-05 | 9.06e-02 | 2771.9 |
| 12 | 7.7e-04 | 7.07e-04 | 9.9e-17 | 7.4e-15 | 1.14e-06 | 1.70e-02 | 14822.1 |
| 14 | 1.9e-04 | 6.02e-05 | 1.5e-16 | 8.2e-15 | 8.29e-09 | 1.44e-03 | 174187.6 |
| 16 | 4.8e-05 | 6.02e-05 | 1.5e-16 | 8.2e-15 | 8.29e-09 | 1.44e-03 | 174188.9 |

## kagome 2x3 cells (tri_0.5 x 4)

N = 18, triangles = 12, decomposition check ||hw A + shift - h|| = 1.3e-14, E0/N = -1.788505

| bits | max angle error | ||A'-A|| | A' Hermitian err | A'^2 = 1 err | E0' - E0 | bound sum hw ||A'-A|| | bound / actual |
|---|---|---|---|---|---|---|---|
| 6 | 4.9e-02 | 5.29e-02 | 1.6e-16 | 8.4e-15 | 5.36e-03 | 1.90e+00 | 355.1 |
| 8 | 1.2e-02 | 3.77e-03 | 9.5e-17 | 7.9e-15 | 2.64e-05 | 1.36e-01 | 5154.8 |
| 10 | 3.1e-03 | 3.77e-03 | 9.5e-17 | 7.9e-15 | 2.64e-05 | 1.36e-01 | 5154.8 |
| 12 | 7.7e-04 | 7.07e-04 | 9.9e-17 | 7.4e-15 | 9.22e-07 | 2.54e-02 | 27598.1 |
| 14 | 1.9e-04 | 6.02e-05 | 1.5e-16 | 8.2e-15 | 6.68e-09 | 2.17e-03 | 324419.3 |
| 16 | 4.8e-05 | 6.02e-05 | 1.5e-16 | 8.2e-15 | 6.68e-09 | 2.17e-03 | 324422.1 |

