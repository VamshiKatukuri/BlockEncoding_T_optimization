# Phase 3: patch designs and magic-state supply

## 1. Exactness obstruction

Optimal block of the isotropic triangle (lambda = 3): distinct entries [-0.333333, 0.0, 0.666667, 1.0]. Entries 1/3, 2/3 are outside Z[1/sqrt2, i], so no exact Clifford+T circuit encodes it (Giles-Selinger).

## 2. Smallest ring-compatible normalization (necessary condition) and its cost

| lattice | L | half-width | ring lambda | alpha penalty | lambda_eff penalty (re-centred) | (same center) |
|---|---|---|---|---|---|---|
| square J1-J2 0.5 | 8 | 0.750 | 1.000 | 1.333 | 1.184 | 1.545 |
| square J1-J2 0.5 | 10 | 0.750 | 1.000 | 1.333 | 1.184 | 1.545 |
| square J1-J2 0.5 | 16 | 0.750 | 1.000 | 1.333 | 1.184 | 1.545 |
| square J1-J2 0.5 | 20 | 0.750 | 1.000 | 1.333 | 1.184 | 1.545 |
| square J1-J2 0.5 | 30 | 0.750 | 1.000 | 1.333 | 1.184 | 1.545 |
| square J1-J2 0.5 | 32 | 0.750 | 1.000 | 1.333 | 1.184 | 1.545 |
| square J1-J2 0.5 | 64 | 0.750 | 1.000 | 1.333 | 1.184 | 1.545 |
| kagome | 4 | 3.000 | 4.000 | 1.333 | 1.164 | 2.091 |
| kagome | 6 | 3.000 | 4.000 | 1.333 | 1.164 | 2.091 |
| kagome | 8 | 3.000 | 4.000 | 1.333 | 1.164 | 2.091 |
| kagome | 10 | 3.000 | 4.000 | 1.333 | 1.164 | 2.091 |
| kagome | 16 | 3.000 | 4.000 | 1.333 | 1.164 | 2.091 |
| kagome | 18 | 3.000 | 4.000 | 1.333 | 1.164 | 2.091 |

## 3. Synthillation upper bound and rotation share (BEST pricing, p = 1e-3, rotations by rule C/B)

Removing every rotation pays for a lambda_eff penalty P only if rotation share > 1 - 1/P.

| lattice | method | L | eps | local-block share | synthillation bound | rotation share |
|---|---|---|---|---|---|---|
| square J1-J2 0.5 | swap | 8 | 0.01 | 0.005 | 1.005 | 0.319 |
| square J1-J2 0.5 | swap | 8 | 0.001 | 0.004 | 1.004 | 0.396 |
| square J1-J2 0.5 | swap | 10 | 0.01 | 0.003 | 1.003 | 0.339 |
| square J1-J2 0.5 | swap | 10 | 0.001 | 0.002 | 1.002 | 0.371 |
| square J1-J2 0.5 | swap | 16 | 0.01 | 0.002 | 1.002 | 0.135 |
| square J1-J2 0.5 | swap | 16 | 0.001 | 0.002 | 1.002 | 0.155 |
| square J1-J2 0.5 | swap | 20 | 0.01 | 0.001 | 1.001 | 0.148 |
| square J1-J2 0.5 | swap | 20 | 0.001 | 0.001 | 1.001 | 0.087 |
| square J1-J2 0.5 | swap | 30 | 0.01 | 0.000 | 1.000 | 0.076 |
| square J1-J2 0.5 | swap | 30 | 0.001 | 0.001 | 1.001 | 0.042 |
| square J1-J2 0.5 | swap | 32 | 0.01 | 0.000 | 1.000 | 0.022 |
| square J1-J2 0.5 | swap | 32 | 0.001 | 0.000 | 1.000 | 0.024 |
| square J1-J2 0.5 | swap | 64 | 0.01 | 0.000 | 1.000 | 0.006 |
| square J1-J2 0.5 | swap | 64 | 0.001 | 0.000 | 1.000 | 0.007 |
| square J1-J2 0.5 | triangle | 8 | 0.01 | 0.337 | 1.507 | 0.298 |
| square J1-J2 0.5 | triangle | 8 | 0.001 | 0.384 | 1.623 | 0.348 |
| square J1-J2 0.5 | triangle | 10 | 0.01 | 0.202 | 1.252 | 0.362 |
| square J1-J2 0.5 | triangle | 10 | 0.001 | 0.217 | 1.278 | 0.396 |
| square J1-J2 0.5 | triangle | 16 | 0.01 | 0.126 | 1.144 | 0.113 |
| square J1-J2 0.5 | triangle | 16 | 0.001 | 0.152 | 1.179 | 0.139 |
| square J1-J2 0.5 | triangle | 20 | 0.01 | 0.088 | 1.097 | 0.161 |
| square J1-J2 0.5 | triangle | 20 | 0.001 | 0.092 | 1.102 | 0.169 |
| square J1-J2 0.5 | triangle | 30 | 0.01 | 0.043 | 1.045 | 0.079 |
| square J1-J2 0.5 | triangle | 30 | 0.001 | 0.024 | 1.025 | 0.042 |
| square J1-J2 0.5 | triangle | 32 | 0.01 | 0.040 | 1.042 | 0.036 |
| square J1-J2 0.5 | triangle | 32 | 0.001 | 0.023 | 1.024 | 0.020 |
| square J1-J2 0.5 | triangle | 64 | 0.01 | 0.006 | 1.006 | 0.005 |
| square J1-J2 0.5 | triangle | 64 | 0.001 | 0.006 | 1.006 | 0.005 |
| kagome | swap | 4 | 0.01 | 0.004 | 1.004 | 0.535 |
| kagome | swap | 4 | 0.001 | 0.003 | 1.003 | 0.617 |
| kagome | swap | 6 | 0.01 | 0.003 | 1.003 | 0.344 |
| kagome | swap | 6 | 0.001 | 0.003 | 1.003 | 0.358 |
| kagome | swap | 8 | 0.01 | 0.002 | 1.002 | 0.288 |
| kagome | swap | 8 | 0.001 | 0.002 | 1.002 | 0.324 |
| kagome | swap | 10 | 0.01 | 0.001 | 1.001 | 0.178 |
| kagome | swap | 10 | 0.001 | 0.001 | 1.001 | 0.198 |
| kagome | swap | 16 | 0.01 | 0.001 | 1.001 | 0.106 |
| kagome | swap | 16 | 0.001 | 0.001 | 1.001 | 0.062 |
| kagome | swap | 18 | 0.01 | 0.000 | 1.000 | 0.066 |
| kagome | swap | 18 | 0.001 | 0.000 | 1.000 | 0.042 |
| kagome | triangle | 4 | 0.01 | 0.243 | 1.321 | 0.200 |
| kagome | triangle | 4 | 0.001 | 0.418 | 1.717 | 0.374 |
| kagome | triangle | 6 | 0.01 | 0.192 | 1.237 | 0.343 |
| kagome | triangle | 6 | 0.001 | 0.213 | 1.271 | 0.389 |
| kagome | triangle | 8 | 0.01 | 0.163 | 1.194 | 0.146 |
| kagome | triangle | 8 | 0.001 | 0.176 | 1.214 | 0.160 |
| kagome | triangle | 10 | 0.01 | 0.105 | 1.117 | 0.191 |
| kagome | triangle | 10 | 0.001 | 0.112 | 1.126 | 0.204 |
| kagome | triangle | 16 | 0.01 | 0.046 | 1.049 | 0.042 |
| kagome | triangle | 16 | 0.001 | 0.030 | 1.031 | 0.027 |
| kagome | triangle | 18 | 0.01 | 0.043 | 1.045 | 0.078 |
| kagome | triangle | 18 | 0.001 | 0.025 | 1.025 | 0.043 |

## 4. Nested-reflection decomposition

| patch | class eigenvalues | b_k | sum b = half-width | eigen rotations (spectral) | PREP rotations (nested) |
|---|---|---|---|---|---|
| J1-J2 triangle r=0.25 (lcvcbe tri_0.25) | [-0.875, -0.375, 0.625] | [0.25, 0.5] | 0.7500 = 0.7500 | 2 | 2 |
| isotropic tetrahedron | [-6.0, -2.0, 6.0] | [2.0, 4.0] | 6.0000 = 6.0000 | 2 | 2 |
