# Phase 4: routing layers as shared-control CCZ batches

## 11. Bell-basis rewrite (dense check)

|| (prod V_i)^dag prod_i CCZ(c, a_i, b_i) (prod V_i) - prod_i CSWAP(c; a_i, b_i) || = 6.3e-16 (m = 1), 2.5e-15 (m = 2), 7.5e-15 (m = 3)

## 12. T count of m CCZ gates sharing one control

| m | 7m | merged phase polynomial 6m + (m mod 2) | PyZX input | PyZX full_reduce | measurement-based 4m |
|---|---|---|---|---|---|
| 1 | 7 | 7 | 7 | 7 | 4 |
| 2 | 14 | 12 | 14 | 12 | 8 |
| 3 | 21 | 19 | 21 | 19 | 12 |
| 4 | 28 | 24 | 28 | 24 | 16 |
| 5 | 35 | 31 | 35 | 31 | 20 |
| 6 | 42 | 36 | 42 | 36 | 24 |
| 7 | 49 | 43 | 49 | 43 | 28 |
| 8 | 56 | 48 | 56 | 48 | 32 |

## 13. CCZ supply needed to run one routing layer per logical step

| lattice | method | L | p | widest layer (CCZ) | CSWAP layers / query | factories for one layer per 15 / 25 / 35 cycles | factory qubits (25 cycles) |
|---|---|---|---|---|---|---|---|
| square J1-J2 0.5 | swap | 8 | 0.001 | 32 | 8 | 128 / 77 / 55 | 3.62e+06 |
| square J1-J2 0.5 | swap | 16 | 0.001 | 128 | 12 | 512 / 308 / 220 | 1.45e+07 |
| square J1-J2 0.5 | swap | 32 | 0.001 | 512 | 16 | 2048 / 1229 / 878 | 5.78e+07 |
| square J1-J2 0.5 | swap | 64 | 0.001 | 2048 | 20 | 8192 / 4916 / 3511 | 2.31e+08 |
| square J1-J2 0.5 | triangle | 8 | 0.001 | 32 | 8 | 128 / 77 / 55 | 3.62e+06 |
| square J1-J2 0.5 | triangle | 16 | 0.001 | 128 | 12 | 512 / 308 / 220 | 1.45e+07 |
| square J1-J2 0.5 | triangle | 32 | 0.001 | 512 | 16 | 2048 / 1229 / 878 | 5.78e+07 |
| square J1-J2 0.5 | triangle | 64 | 0.001 | 2048 | 20 | 8192 / 4916 / 3511 | 2.31e+08 |
| kagome | swap | 4 | 0.001 | 24 | 8 | 96 / 58 / 42 | 2.73e+06 |
| kagome | swap | 8 | 0.001 | 96 | 12 | 384 / 231 / 165 | 1.09e+07 |
| kagome | swap | 16 | 0.001 | 384 | 16 | 1536 / 922 / 659 | 4.33e+07 |
| kagome | triangle | 4 | 0.001 | 24 | 8 | 96 / 58 / 42 | 2.73e+06 |
| kagome | triangle | 8 | 0.001 | 96 | 12 | 384 / 231 / 165 | 1.09e+07 |
| kagome | triangle | 16 | 0.001 | 384 | 16 | 1536 / 922 / 659 | 4.33e+07 |
| square J1-J2 0.5 | swap | 8 | 0.0001 | 32 | 8 | 78 / 47 / 34 | 5.83e+05 |
| square J1-J2 0.5 | swap | 16 | 0.0001 | 128 | 12 | 309 / 185 / 133 | 2.29e+06 |
| square J1-J2 0.5 | swap | 32 | 0.0001 | 512 | 16 | 1233 / 740 / 529 | 9.18e+06 |
| square J1-J2 0.5 | swap | 64 | 0.0001 | 2048 | 20 | 4929 / 2958 / 2113 | 3.67e+07 |
| square J1-J2 0.5 | triangle | 8 | 0.0001 | 32 | 8 | 78 / 47 / 34 | 5.83e+05 |
| square J1-J2 0.5 | triangle | 16 | 0.0001 | 128 | 12 | 309 / 185 / 133 | 2.29e+06 |
| square J1-J2 0.5 | triangle | 32 | 0.0001 | 512 | 16 | 1233 / 740 / 529 | 9.18e+06 |
| square J1-J2 0.5 | triangle | 64 | 0.0001 | 2048 | 20 | 4929 / 2958 / 2113 | 3.67e+07 |
| kagome | swap | 4 | 0.0001 | 24 | 8 | 58 / 35 / 25 | 4.34e+05 |
| kagome | swap | 8 | 0.0001 | 96 | 12 | 232 / 139 / 100 | 1.72e+06 |
| kagome | swap | 16 | 0.0001 | 384 | 16 | 925 / 555 / 397 | 6.88e+06 |
| kagome | triangle | 4 | 0.0001 | 24 | 8 | 58 / 35 / 25 | 4.34e+05 |
| kagome | triangle | 8 | 0.0001 | 96 | 12 | 232 / 139 / 100 | 1.72e+06 |
| kagome | triangle | 16 | 0.0001 | 384 | 16 | 925 / 555 / 397 | 6.88e+06 |
