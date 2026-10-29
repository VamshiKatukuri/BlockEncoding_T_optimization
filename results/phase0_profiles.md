# Phase 0: per-query magic demand profiles

Counts per query (both routing directions, PREPARE and its inverse). CSWAP = 1 Toffoli-class gate.

| lattice | method | L | N | CSWAP routing | AND addr | PREP (AND, T, rot) | local (AND+Tof, T, rot) | routing share of non-Clifford gates | alpha/N | lambda_eff/N | T/query (lcvcbe rule) | lcvcbe ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| square J1-J2 0.5 | swap | 8 | 64 | 126 | 12 | 0, 0, 2 | 1, 0, 0 | 0.894 | 6.000 | 5.914 | 738 | 738 |
| square J1-J2 0.5 | swap | 10 | 100 | 198 | 16 | 46, 0, 4 | 1, 0, 0 | 0.747 | 6.000 | 5.914 | 1352 | 1352 |
| square J1-J2 0.5 | swap | 16 | 256 | 510 | 16 | 0, 0, 2 | 1, 0, 0 | 0.964 | 6.000 | 5.914 | 2302 | 2302 |
| square J1-J2 0.5 | swap | 20 | 400 | 798 | 20 | 46, 0, 4 | 1, 0, 0 | 0.918 | 6.000 | 5.914 | 3788 | - |
| square J1-J2 0.5 | swap | 30 | 900 | 1798 | 20 | 64, 0, 4 | 1, 0, 0 | 0.953 | 6.000 | 5.914 | 7860 | - |
| square J1-J2 0.5 | swap | 32 | 1024 | 2046 | 20 | 0, 0, 2 | 1, 0, 0 | 0.989 | 6.000 | 5.914 | 8474 | 8474 |
| square J1-J2 0.5 | swap | 64 | 4096 | 8190 | 24 | 0, 0, 2 | 1, 0, 0 | 0.997 | 6.000 | 5.914 | 33094 | 33094 |
| square J1-J2 0.5 | triangle | 8 | 64 | 126 | 12 | 0, 0, 0 | 4, 8, 4 | 0.818 | 3.000 | 2.247 | 868 | 868 |
| square J1-J2 0.5 | triangle | 10 | 100 | 198 | 16 | 28, 0, 4 | 4, 8, 4 | 0.756 | 3.000 | 2.247 | 1632 | 1632 |
| square J1-J2 0.5 | triangle | 16 | 256 | 510 | 16 | 0, 0, 0 | 4, 8, 4 | 0.941 | 3.000 | 2.247 | 2448 | 2448 |
| square J1-J2 0.5 | triangle | 20 | 400 | 798 | 20 | 28, 0, 4 | 4, 8, 4 | 0.921 | 3.000 | 2.247 | 4088 | - |
| square J1-J2 0.5 | triangle | 30 | 900 | 1798 | 20 | 46, 0, 4 | 4, 8, 4 | 0.954 | 3.000 | 2.247 | 8192 | - |
| square J1-J2 0.5 | triangle | 32 | 1024 | 2046 | 20 | 0, 0, 0 | 4, 8, 4 | 0.983 | 3.000 | 2.247 | 8628 | 8628 |
| square J1-J2 0.5 | triangle | 64 | 4096 | 8190 | 24 | 0, 0, 0 | 4, 8, 4 | 0.995 | 3.000 | 2.247 | 33252 | 33252 |
| kagome | swap | 4 | 48 | 94 | 8 | 10, 0, 4 | 1, 0, 0 | 0.803 | 4.000 | 3.992 | 816 |  |
| kagome | swap | 6 | 108 | 214 | 12 | 28, 0, 4 | 1, 0, 0 | 0.826 | 4.000 | 3.992 | 1304 |  |
| kagome | swap | 8 | 192 | 382 | 12 | 10, 0, 4 | 1, 0, 0 | 0.934 | 4.000 | 3.992 | 2008 |  |
| kagome | swap | 10 | 300 | 598 | 16 | 40, 0, 4 | 1, 0, 0 | 0.907 | 4.000 | 3.992 | 2944 |  |
| kagome | swap | 16 | 768 | 1534 | 16 | 10, 0, 4 | 1, 0, 0 | 0.980 | 4.000 | 3.992 | 6656 |  |
| kagome | swap | 18 | 972 | 1942 | 20 | 46, 0, 4 | 1, 0, 0 | 0.965 | 4.000 | 3.992 | 8392 |  |
| kagome | triangle | 4 | 48 | 90 | 16 | 0, 0, 0 | 4, 8, 4 | 0.738 | 2.000 | 0.960 | 740 |  |
| kagome | triangle | 6 | 108 | 210 | 24 | 22, 0, 4 | 4, 8, 4 | 0.761 | 2.000 | 0.960 | 1672 |  |
| kagome | triangle | 8 | 192 | 378 | 24 | 0, 0, 0 | 4, 8, 4 | 0.904 | 2.000 | 0.960 | 1940 |  |
| kagome | triangle | 10 | 300 | 594 | 32 | 28, 0, 4 | 4, 8, 4 | 0.881 | 2.000 | 0.960 | 3292 |  |
| kagome | triangle | 16 | 768 | 1530 | 32 | 0, 0, 0 | 4, 8, 4 | 0.970 | 2.000 | 0.960 | 6612 |  |
| kagome | triangle | 18 | 972 | 1938 | 40 | 40, 0, 4 | 4, 8, 4 | 0.951 | 2.000 | 0.960 | 8856 |  |

Routing layers (power-of-two trees): each layer j of a tree holds 2^j CSWAPs sharing one address bit;
all trees run in parallel. Example rows:
- square J1-J2 0.5 swap L=8: widths per direction [32, 16, 8, 4], 8 CSWAP layers per query
- square J1-J2 0.5 swap L=16: widths per direction [128, 64, 32, 16, 8, 4], 12 CSWAP layers per query
- square J1-J2 0.5 triangle L=8: widths per direction [32, 16, 8, 4], 8 CSWAP layers per query
- square J1-J2 0.5 triangle L=16: widths per direction [128, 64, 32, 16, 8, 4], 12 CSWAP layers per query
- kagome swap L=8: widths per direction [96, 48, 24, 12, 6, 3], 12 CSWAP layers per query
- kagome swap L=16: widths per direction [384, 192, 96, 48, 24, 12, 6, 3], 16 CSWAP layers per query
- kagome triangle L=8: widths per direction [96, 48, 24, 12, 6, 3], 12 CSWAP layers per query
- kagome triangle L=16: widths per direction [384, 192, 96, 48, 24, 12, 6, 3], 16 CSWAP layers per query
