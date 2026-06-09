# Random Source Audit (2026-05-07)

- Fixed random seed: `20260507`
- Sampling unit: one random supported source row from each of 12 randomly selected galaxies.
- Supported in this pass: local VizieR TSV snapshots and local Simon & Geha jsimon data files.
- Tolerances: coordinates <= 5e-4 deg, velocities/errors <= 5e-3 km/s.

## Summary

audit_status
PASS    12

## Sampled Records

### Antlia 2 / `5432925619884476928`
- CSV: `data/galaxies/01_Antlia_2.csv`
- Source: `J/ApJ/921/32/table4`; link: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/921/32
- RA: source=145.14466 csv=145.14466 match=True
- Dec: source=-38.22482 csv=-38.22482 match=True
- v_los: source=30.45 csv=30.45 match=True
- v_los_err: source=3.29 csv_down=3.29 csv_up=3.29 down_match=True up_match=True
- member: source_raw= mapped=2 csv=2 match=True
- Status: **PASS** 

### Ursa Major I / `96_513`
- CSV: `data/galaxies/25_Ursa_Major_I.csv`
- Source: `jsimon_data/UMaI_feh.dat`; link: https://users.obs.carnegiescience.edu/jsimon/data.html
- RA: source=158.48374 csv=158.48374 match=True
- Dec: source=51.893743 csv=51.893743 match=True
- v_los: source=16.76 csv=16.76 match=True
- v_los_err: source=2.69 csv_down=2.69 csv_up=2.69 down_match=True up_match=True
- member: source_raw=0 mapped=0 csv=0 match=True
- Status: **PASS** 

### Leo V / `Leo5_1046`
- CSV: `data/galaxies/15_Leo_V.csv`
- Source: `J/ApJ/920/92/table6`; link: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/920/92
- RA: source=172.7569167 csv=172.7569167 match=True
- Dec: source=2.1903056 csv=2.1903056 match=True
- v_los: source=173.6 csv=173.6 match=True
- v_los_err: source=0.9 csv_down=0.9 csv_up=0.9 down_match=True up_match=True
- member: source_raw=1 mapped=1 csv=1 match=True
- Status: **PASS** 

### Tucana 3 / `J235709.04-593400.3`
- CSV: `data/galaxies/23_Tucana_3.csv`
- Source: `J/ApJ/838/11/table2`; link: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/838/11
- RA: source=359.28767 csv=359.28767 match=True
- Dec: source=-59.56675 csv=-59.56675 match=True
- v_los: source=-193.88 csv=-193.88 match=True
- v_los_err: source=7.22 csv_down=7.22 csv_up=7.22 down_match=True up_match=True
- member: source_raw=0 mapped=0 csv=0 match=True
- Status: **PASS** 

### Ursa Major II / `200_329`
- CSV: `data/galaxies/26_Ursa_Major_II.csv`
- Source: `jsimon_data/UMa2_feh.dat`; link: https://users.obs.carnegiescience.edu/jsimon/data.html
- RA: source=132.45797 csv=132.45797 match=True
- Dec: source=63.153612 csv=63.153612 match=True
- v_los: source=16.6 csv=16.6 match=True
- v_los_err: source=2.28 csv_down=2.28 csv_up=2.28 down_match=True up_match=True
- member: source_raw=0 mapped=0 csv=0 match=True
- Status: **PASS** 

### Bootes I / `Boo1_44`
- CSV: `data/galaxies/02_Bootes_I.csv`
- Source: `J/ApJ/920/92/table6`; link: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/920/92
- RA: source=209.98325 csv=209.98325 match=True
- Dec: source=14.5737778 csv=14.5737778 match=True
- v_los: source=114.8 csv=114.8 match=True
- v_los_err: source=0.8 csv_down=0.8 csv_up=0.8 down_match=True up_match=True
- member: source_raw=1 mapped=1 csv=1 match=True
- Status: **PASS** 

### Triangulum II / `166.0`
- CSV: `data/galaxies/21_Triangulum_II.csv`
- Source: `J/ApJ/838/83/table2`; link: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/838/83
- RA: source=33.297583333333336 csv=33.297583333333336 match=True
- Dec: source=36.209 csv=36.209 match=True
- v_los: source=11.0 csv=11.0 match=True
- v_los_err: source=5.9 csv_down=5.9 csv_up=5.9 down_match=True up_match=True
- member: source_raw=N mapped=0 csv=0 match=True
- Status: **PASS** 

### Tucana 2 / `TucII-022`
- CSV: `data/galaxies/22_Tucana_2.csv`
- Source: `J/AJ/165/55/table6`; link: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/AJ/165/55
- RA: source=343.08908 csv=343.08908 match=True
- Dec: source=-58.51869 csv=-58.51869 match=True
- v_los: source=-120.8 csv=-120.8 match=True
- v_los_err: source=1.1 csv_down=1.1 csv_up=1.1 down_match=True up_match=True
- member: source_raw= mapped=2 csv=2 match=True
- Status: **PASS** 

### Tucana 4 / `J000501.31-604504.2`
- CSV: `data/galaxies/24_Tucana_4.csv`
- Source: `J/ApJ/892/137/table3`; link: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/892/137
- RA: source=1.25544 csv=1.25544 match=True
- Dec: source=-60.75118 csv=-60.75118 match=True
- v_los: source=155.2 csv=155.2 match=True
- v_los_err: source=1.6 csv_down=1.6 csv_up=1.6 down_match=True up_match=True
- member: source_raw=0 mapped=0 csv=0 match=True
- Status: **PASS** 

### Crater 2 / `3.5440043245630147e+18`
- CSV: `data/galaxies/06_Crater_2.csv`
- Source: `J/ApJ/921/32/table5`; link: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/921/32
- RA: source=177.27017 csv=177.27017 match=True
- Dec: source=-18.21341 csv=-18.21341 match=True
- v_los: source=87.65 csv=87.65 match=True
- v_los_err: source=3.67 csv_down=3.67 csv_up=3.67 down_match=True up_match=True
- member: source_raw= mapped=2 csv=2 match=True
- Status: **PASS** 

### Canes Venatici I / `542_907`
- CSV: `data/galaxies/03_Canes_Venatici_I.csv`
- Source: `jsimon_data/CVn1_with_feh.dat`; link: https://users.obs.carnegiescience.edu/jsimon/data.html
- RA: source=201.97332 csv=201.97332 match=True
- Dec: source=33.531877 csv=33.531877 match=True
- v_los: source=27.4 csv=27.4 match=True
- v_los_err: source=4.71 csv_down=4.71 csv_up=4.71 down_match=True up_match=True
- member: source_raw=1 mapped=1 csv=1 match=True
- Status: **PASS** 

### Canes Venatici II / `499_94`
- CSV: `data/galaxies/04_Canes_Venatici_II.csv`
- Source: `jsimon_data/CVnII_feh.dat`; link: https://users.obs.carnegiescience.edu/jsimon/data.html
- RA: source=194.28579 csv=194.28579 match=True
- Dec: source=34.344509 csv=34.344509 match=True
- v_los: source=-128.05 csv=-128.05 match=True
- v_los_err: source=2.87 csv_down=2.87 csv_up=2.87 down_match=True up_match=True
- member: source_raw=1 mapped=1 csv=1 match=True
- Status: **PASS** 
