---
---

# Pleural Dilution Factor Analysis

## Background

Initial pleural protein results, `GB390725-RB_pleural fluid_NPX.csv` showed poor correlation with gold standard ELISA measured IL6 values. 

Several re-runs were performed in order to arrive at the correct dilution:
1. A dilution series analysis with the same samples run repeatedly across several dilution series (`GB390725-RB_dilutionseries_Extended_2026-05-07.csv`)
2. 1 plate run at 16x dilution when this was first selected at the appropriate dilution, but which didn't show good correlation with ELISA IL6 (`GB390725-RB_1_16_Extended_2026-07-06.csv`)
3. a re-run of 1 plate run at 16x dilution (above) because Kate Orr reported that the previous run used a bad batch of Olink reagents (`GB390725-RB_1_1in16_rpt_Standard_2026-08-27.csv`). This run still didn't show good correlation with ELISA IL6 so we decided to re-run the full experiment at 64x dilution

## `run-dilution.r`

This has been generated as a copy from the original full directory `run-all.r` on 20 October 2026.

As I didn't expect to re-generate these reports again, I didn't update the directory structure to be able re-run these directly from `run-dilution.r`.

However, all relevant reports were copied into `docs/dilution` so that a record of the analysis woueld be preserved

