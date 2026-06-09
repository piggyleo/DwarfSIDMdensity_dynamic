# Member Count Discrepancy Audit Log (2026-05-06)

本日志按用户定义规则检查 `n_member_flag_1_minus_hayashi2023` 与 `n_member_flag_12_minus_hayashi2023`。原始星系 CSV 未在本轮修改。

## Summary

- `case1`: 0 galaxies
- `case2`: 2 galaxies
- `case3`: 10 galaxies
- `case4`: 3 galaxies
- `skip`: 12 galaxies
- Additional follow-up among skipped galaxies: 1 galaxy (`Reticulum II`)

## Checked Galaxies

### Antlia 2 (case4)

- Counts: Hayashi=283; flag1=0 (diff=-283); flag2=729; flag12=729 (diff=446).
- Duplicate check: unique flag1 coords=0; duplicate flag1 rows=0; unique flag12 coords=726; duplicate flag12 rows=6.
- MemProb thresholds among all star rows: >=0.1:292; >=0.5:290; >=0.75:287; >=0.8:287; >=0.85:285; >=0.9:284; >=0.95:284; >=0.99:279.
- Source counts: J/ApJ/921/32/table4: flag0=0, flag1=0, flag2=508; J/MNRAS/488/2743/table2: flag0=0, flag1=0, flag2=221.
- Result: 触发规则4。member_flag=2 不是“全部成员星”：Ji et al. 2021 的 Mm 是成员概率，Torrealba et al. 2019 table2 没有硬成员列。当前 729 行不应直接作为 Hayashi 成员样本。
- Likely reason: 当前表保留了概率样本/建模样本的全部逐星记录；Hayashi 的 n=283 应来自额外筛选。Ji 表中 Mm>=0.90 有 284 行，接近但不完全等于 283；Torrealba 表没有硬成员标记。
- Recommended next action: 后续建模前需回到 Ji/Torrealba 正文确认 Antlia 2 的最终 283 星筛选口径；在未确认前保持 member_flag=2，不把全部 flag2 计入似然。
- Evidence notes: VizieR J/ApJ/921/32 ReadMe: Mm=[0/1] membership probability；J/MNRAS/488/2743 table2: spectroscopic modelling sample, no membership code.
- Links: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/921/32; https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/MNRAS/488/2743

### Bootes I (case3)

- Counts: Hayashi=66; flag1=73 (diff=7); flag2=0; flag12=73 (diff=7).
- Duplicate check: unique flag1 coords=73; duplicate flag1 rows=0; unique flag12 coords=73; duplicate flag12 rows=0.
- MemProb thresholds among member_flag=1 rows: >=0.1:71; >=0.5:70; >=0.75:70; >=0.8:70; >=0.85:69; >=0.9:69; >=0.95:66; >=0.99:57.
- MemProb thresholds among all star rows: >=0.1:73; >=0.5:72; >=0.75:72; >=0.8:72; >=0.85:71; >=0.9:70; >=0.95:67; >=0.99:57.
- Source counts: J/ApJ/920/92/table6: flag0=45, flag1=73, flag2=0.
- Result: 触发规则3。member_flag=1 有 73 行、无重复；但 Jenkins 表的 MemProb>=0.95 正好为 66，与 Hayashi n=66 一致。
- Likely reason: 当前 member_flag=1 采用 Jenkins 的 subjective Member=1，范围比 Hayashi 采用的高概率成员样本更宽。
- Recommended next action: 保留原始 Member=1，同时新增/派生 Hayashi 等价筛选：Jenkins Bootes I 使用 MemProb>=0.95。
- Evidence notes: VizieR J/ApJ/920/92 ReadMe: Member is subjective membership classification; MemProb is membership probability.
- Links: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/920/92

### Crater 2 (case4)

- Counts: Hayashi=141; flag1=0 (diff=-141); flag2=207; flag12=207 (diff=66).
- Duplicate check: unique flag1 coords=0; duplicate flag1 rows=0; unique flag12 coords=207; duplicate flag12 rows=0.
- MemProb thresholds among all star rows: >=0.1:141; >=0.5:141; >=0.75:141; >=0.8:141; >=0.85:141; >=0.9:141; >=0.95:141; >=0.99:141.
- Source counts: J/ApJ/921/32/table5: flag0=0, flag1=0, flag2=207.
- Result: 触发规则4。member_flag=2 不是全部成员星；Ji 表提供的是成员概率。Mm>=0.1 到 Mm>=0.99 的计数均为 141，精确等于 Hayashi n=141。
- Likely reason: 当前表保留了 Ji 概率表的全部 207 行；Hayashi 很可能只取非零/高概率成员星。
- Recommended next action: 后续建模可把 Crater 2 的 Hayashi 等价样本定义为 Ji Mm>=0.1（等价于非零概率样本），但正式改动前建议人工确认源文献正文。
- Evidence notes: VizieR J/ApJ/921/32 ReadMe: Crater 2 table has Mm membership probability, not a hard flag.
- Links: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/921/32

### Draco 2 (case3)

- Counts: Hayashi=9; flag1=14 (diff=5); flag2=0; flag12=14 (diff=5).
- Duplicate check: unique flag1 coords=14; duplicate flag1 rows=0; unique flag12 coords=14; duplicate flag12 rows=0.
- MemProb thresholds among member_flag=1 rows: >=0.1:13; >=0.5:10; >=0.75:10; >=0.8:10; >=0.85:9; >=0.9:6; >=0.95:5; >=0.99:2.
- MemProb thresholds among all star rows: >=0.1:18; >=0.5:11; >=0.75:10; >=0.8:10; >=0.85:9; >=0.9:6; >=0.95:5; >=0.99:2.
- Source counts: arXiv:1807.10655/table_spectro: flag0=37, flag1=14, flag2=0.
- Result: 触发规则3。当前源表中 Member=Y 的 14 行无重复；源文献正文也说明用 CaHK/CMD 清理后得到 14 likely/dynamical members。Hayashi Table 1 为 9，差异未由重复行解释。
- Likely reason: Hayashi 可能采用了更严格的子样本或旧/不同成员口径；当前已抽取源文献硬成员，但与 Hayashi n 不一致。
- Recommended next action: 暂不改 flag；后续需要专门追踪 Hayashi 为何取 n=9，尤其检查是否有只用于 Jeans 似然的速度质量/空间/成员概率截断。
- Evidence notes: Longeard et al. 2018 source text: final spectroscopic cleaning isolates 14 likely member stars; current CSV 14 Member=Y.
- Links: https://arxiv.org/abs/1807.10655

### Grus 1 (case3)

- Counts: Hayashi=8; flag1=16 (diff=8); flag2=0; flag12=16 (diff=8).
- Duplicate check: unique flag1 coords=8; duplicate flag1 rows=14; unique flag12 coords=8; duplicate flag12 rows=14.
- Source counts: arXiv:2206.04580/table_gru1spec: flag0=63, flag1=16, flag2=0.
- Result: 触发规则3。member_flag=1 原始行数 16，但按坐标去重后为 8，正好等于 Hayashi n=8。
- Likely reason: 当前 Grus 1 表包含多历元速度测量行（MJD 不同），不是唯一恒星表。
- Recommended next action: 后续建模应按唯一恒星合并多历元速度，或按源文献方式构造联合/加权速度；不要把 16 行直接当作 16 颗成员星。
- Evidence notes: Chiti et al. 2022 source text: summary table gives N_spectroscopic_members=8; raw table contains MJD-tagged repeated observations.
- Links: https://arxiv.org/abs/2206.04580

### Grus 2 (case3)

- Counts: Hayashi=19; flag1=49 (diff=30); flag2=115; flag12=164 (diff=145).
- Duplicate check: unique flag1 coords=21; duplicate flag1 rows=43; unique flag12 coords=135; duplicate flag12 rows=44.
- Source counts: J/ApJ/892/137/table2: flag0=0, flag1=0, flag2=115; J/ApJ/892/137/table3: flag0=186, flag1=49, flag2=0.
- Result: 触发规则3。member_flag=1 原始行数 49；按坐标/ID 去重后为 21，重复观测解释了大部分差异，但仍比 Hayashi n=19 多 2 颗唯一成员。
- Likely reason: Simon et al. 2020 table3 是速度测量表，含多次观测；当前 Mm=1 的唯一星数仍略高于 Hayashi，可能有二元星/质量/最终动力学样本剔除。
- Recommended next action: 先在建模样本中合并重复观测；再逐星核查 21 个唯一 Mm=1 中哪 2 个未进入 Hayashi 的 19 星样本。
- Evidence notes: VizieR J/ApJ/892/137 ReadMe: table3 is IMACS velocity and metallicity measurements; Nrv gives number of HRV measurements; Mm is membership code.
- Links: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/892/137

### Hercules (case3)

- Counts: Hayashi=18; flag1=21 (diff=3); flag2=0; flag12=21 (diff=3).
- Duplicate check: unique flag1 coords=21; duplicate flag1 rows=0; unique flag12 coords=21; duplicate flag12 rows=0.
- MemProb thresholds among member_flag=1 rows: >=0.1:21; >=0.5:21; >=0.75:21; >=0.8:21; >=0.85:21; >=0.9:21; >=0.95:21; >=0.99:21.
- Source counts: jsimon_data/Herc_feh.dat: flag0=0, flag1=21, flag2=0.
- Result: 触发规则3。当前 21 个 member_flag=1 无重复；Hayashi n=18，差异不是重复观测造成。
- Likely reason: 当前速度由 jsimon_data/Herc_feh.dat 回填，且 CSV 星表源自 Kirby metallicity overlap；Hayashi Table 1 对 Hercules 列的是 refs 2,5,19，成员数可能来自更窄的 18 星动力学/成员样本。
- Recommended next action: 暂不删星；后续需要定位 Hercules 的 18 星原始成员表或 Hayashi 采用的筛选口径。
- Evidence notes: 本地 CSV notes 显示 velocity from Simon & Geha via jsimon_data；无 member_flag=2 和无重复成员 ID。
- Links: https://users.obs.carnegiescience.edu/jsimon/data.html

### Leo IV (case3)

- Counts: Hayashi=18; flag1=20 (diff=2); flag2=0; flag12=20 (diff=2).
- Duplicate check: unique flag1 coords=20; duplicate flag1 rows=0; unique flag12 coords=20; duplicate flag12 rows=0.
- MemProb thresholds among member_flag=1 rows: >=0.1:20; >=0.5:19; >=0.75:18; >=0.8:17; >=0.85:17; >=0.9:17; >=0.95:15; >=0.99:10.
- MemProb thresholds among all star rows: >=0.1:24; >=0.5:23; >=0.75:19; >=0.8:18; >=0.85:18; >=0.9:18; >=0.95:16; >=0.99:11.
- Source counts: J/ApJ/920/92/table6: flag0=84, flag1=20, flag2=0.
- Result: 触发规则3。Jenkins subjective Member=1 有 20 行、无重复；MemProb>=0.75 为 18，与 Hayashi n=18 一致。
- Likely reason: 当前 member_flag=1 采用 subjective Member=1，Hayashi 更可能采用概率阈值后的高概率成员样本。
- Recommended next action: 后续建模前派生 Hayashi 等价筛选：Leo IV 可从 MemProb>=0.75 开始人工确认。
- Evidence notes: VizieR J/ApJ/920/92 ReadMe: Member is subjective membership classification; MemProb is membership probability.
- Links: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/920/92

### Leo V (case3)

- Counts: Hayashi=7; flag1=11 (diff=4); flag2=0; flag12=11 (diff=4).
- Duplicate check: unique flag1 coords=11; duplicate flag1 rows=0; unique flag12 coords=11; duplicate flag12 rows=0.
- MemProb thresholds among member_flag=1 rows: >=0.1:9; >=0.5:9; >=0.75:7; >=0.8:7; >=0.85:7; >=0.9:7; >=0.95:7; >=0.99:5.
- MemProb thresholds among all star rows: >=0.1:11; >=0.5:9; >=0.75:7; >=0.8:7; >=0.85:7; >=0.9:7; >=0.95:7; >=0.99:5.
- Source counts: J/ApJ/920/92/table6: flag0=94, flag1=11, flag2=0.
- Result: 触发规则3。Jenkins subjective Member=1 有 11 行、无重复；有 MemProb 的行中 MemProb>=0.95 为 7，与 Hayashi n=7 一致。
- Likely reason: 当前 member_flag=1 比 Hayashi 的高概率成员样本更宽；另有两个 subjective members 缺 MemProb。
- Recommended next action: 后续建模前派生 Hayashi 等价筛选：Leo V 可从 MemProb>=0.95 开始人工确认。
- Evidence notes: VizieR J/ApJ/920/92 ReadMe: Member is subjective membership classification; MemProb is membership probability.
- Links: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/920/92

### Segue 1 (case2)

- Counts: Hayashi=71; flag1=71 (diff=0); flag2=129; flag12=200 (diff=129).
- Duplicate check: unique flag1 coords=71; duplicate flag1 rows=0; unique flag12 coords=137; duplicate flag12 rows=104.
- MemProb thresholds among member_flag=1 rows: >=0.1:71; >=0.5:71; >=0.75:70; >=0.8:70; >=0.85:69; >=0.9:67; >=0.95:62; >=0.99:38.
- MemProb thresholds among all star rows: >=0.1:76; >=0.5:74; >=0.75:71; >=0.8:71; >=0.85:70; >=0.9:68; >=0.95:63; >=0.99:38.
- Source counts: J/ApJ/733/46/table3: flag0=322, flag1=71, flag2=129.
- Result: 触发规则2。member_flag=1 正好为 71，与 Hayashi n=71 一致；member_flag=2 来自源表 Mm 空白/非硬成员状态，文献没有说明这些都是成员星。
- Likely reason: flag2 是未给硬成员标记或概率但未被 Mm=1 接纳的观测行；不应计入 Hayashi 等价成员样本。
- Recommended next action: 保持 flag2 作为不确定/未分类行；后续似然默认只用 member_flag=1。
- Evidence notes: VizieR J/ApJ/733/46 ReadMe: Mm=[0/1]? member status, with blank values possible; EM/Bayesian probabilities are separate columns.
- Links: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/733/46

### Segue 2 (case3)

- Counts: Hayashi=26; flag1=27 (diff=1); flag2=126; flag12=153 (diff=127).
- Duplicate check: unique flag1 coords=26; duplicate flag1 rows=2; unique flag12 coords=135; duplicate flag12 rows=35.
- Source counts: J/ApJ/770/16/table2: flag0=590, flag1=22, flag2=35; J/MNRAS/397/1748/table2: flag0=326, flag1=5, flag2=91.
- Result: 触发规则3。member_flag=1 原始行数 27；按坐标去重后为 26，正好等于 Hayashi n=26。存在 1 个跨源重复成员（Kirby ref29 与 Belokurov ref28 同一坐标）。
- Likely reason: 当前合并了两个来源的成员表，重复计入了同一颗成员星；flag2 多数是 ?/B/未知口径，不是硬成员样本。
- Recommended next action: 后续建模时对 ref28/ref29 跨表成员按坐标去重；flag2 不计入默认 Hayashi 等价似然样本。
- Evidence notes: VizieR J/ApJ/770/16 ReadMe: Mm Y/N/?/B; J/MNRAS/397/1748 source has overlapping member at RA=34.77054, Dec=20.12094.
- Links: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/770/16; https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/MNRAS/397/1748

### Triangulum II (case2)

- Counts: Hayashi=13; flag1=13 (diff=0); flag2=1; flag12=14 (diff=1).
- Duplicate check: unique flag1 coords=13; duplicate flag1 rows=0; unique flag12 coords=14; duplicate flag12 rows=0.
- Source counts: J/ApJ/838/83/table2: flag0=20, flag1=13, flag2=1.
- Result: 触发规则2。member_flag=1 为 13，与 Hayashi n=13 一致；唯一 member_flag=2 是 [MIC2016] 25，源表标注 f_Mm=b 且 HRV 为空。
- Likely reason: flag2 是源表中的特殊/不确定行，不是成员星；它也没有速度，不能进入动力学似然。
- Recommended next action: 保持该行 member_flag=2；后续似然只用 13 个 member_flag=1 且有速度的成员。
- Evidence notes: VizieR J/ApJ/838/83 ReadMe: Mm member status has Y:13 and N:20; f_Mm=[b] flags M16 25.
- Links: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/838/83

### Tucana 2 (case4)

- Counts: Hayashi=19; flag1=0 (diff=-19); flag2=60; flag12=60 (diff=41).
- Duplicate check: unique flag1 coords=0; duplicate flag1 rows=0; unique flag12 coords=19; duplicate flag12 rows=55.
- Source counts: J/AJ/165/55/table6: flag0=0, flag1=0, flag2=60.
- Result: 触发规则4。member_flag=2 原始行数 60；按 star_id/坐标去重后为 19，正好等于 Hayashi n=19。源表 table6 无硬成员列。
- Likely reason: Chiti et al. 2023 table6 是速度测量表，包含 19 颗星的多次速度测量；行数差异来自重复观测，不是额外成员。
- Recommended next action: 保持原始多历元行；后续建模需按 19 个唯一 star_id 合并速度，或实现多历元观测合并接口。
- Evidence notes: VizieR J/AJ/165/55 ReadMe: table6 title is Velocity measurements; columns include Name/MJD/Inst/RVel but no membership flag.
- Links: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/AJ/165/55

### Tucana 3 (case3)

- Counts: Hayashi=26; flag1=44 (diff=18); flag2=0; flag12=44 (diff=18).
- Duplicate check: unique flag1 coords=26; duplicate flag1 rows=28; unique flag12 coords=26; duplicate flag12 rows=28.
- Source counts: J/ApJ/838/11/table2: flag0=107, flag1=44, flag2=0.
- Result: 触发规则3。member_flag=1 原始行数 44；按 star_id/坐标去重后为 26，正好等于 Hayashi n=26。
- Likely reason: Simon et al. 2017 table2 是逐次速度测量表，同一成员星有多次观测。
- Recommended next action: 后续建模前按唯一恒星合并多次速度测量；不要把 44 行当作 44 颗成员。
- Evidence notes: 当前 CSV 中 10 个成员 ID 有重复；唯一 member_flag=1 坐标数为 26。
- Links: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/838/11

### Tucana 4 (case3)

- Counts: Hayashi=11; flag1=39 (diff=28); flag2=63; flag12=102 (diff=91).
- Duplicate check: unique flag1 coords=11; duplicate flag1 rows=36; unique flag12 coords=74; duplicate flag12 rows=36.
- Source counts: J/ApJ/892/137/table2: flag0=0, flag1=0, flag2=63; J/ApJ/892/137/table3: flag0=184, flag1=39, flag2=0.
- Result: 触发规则3。member_flag=1 原始行数 39；按 star_id/坐标去重后为 11，正好等于 Hayashi n=11。member_flag=2 来自 AAT velocity table，非默认分析样本。
- Likely reason: Simon et al. 2020 table3 包含多次 IMACS 观测；table2 AAT 数据被保留为不确定行但不应计入 Hayashi 等价成员样本。
- Recommended next action: 后续建模使用 table3 的 11 个唯一 Mm=1 成员并合并重复观测；table2 flag2 暂不纳入默认似然。
- Evidence notes: VizieR J/ApJ/892/137 ReadMe: table3 has Nrv and Mm; table2 is AAT velocity measurements.
- Links: https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/892/137

## Additional Follow-Up

### Reticulum II (member-definition mismatch)

- Counts: Hayashi=25; flag1=18 (diff=-7); flag2=0; flag12=18 (diff=-7).
- Original rule status: This galaxy fell under rule 5 in the first pass, but was rechecked later because both `diff1` and `diff12` are negative.
- Source counts: Koposov et al. 2015b Table 2 contains 25 Reticulum II candidate observation rows; the source text identifies 18 confirmed members.
- Result: No star-row omission was found. The current CSV already contains all 25 candidate rows from the source table, but only the 18 source-confirmed members have `member_flag=1`.
- Likely reason: Hayashi Table 1 reports `n=25`, which appears to correspond to the candidate observation count or a different sample definition, not the source paper's confirmed-member count.
- Recommended next action: Do not add new rows. Before reproducing Hayashi for Reticulum II, confirm whether the modeling sample should use the 18 confirmed members or all 25 candidates.
- Evidence notes: Koposov et al. 2015b states that Reticulum II has 18 identified members; its Table 2 lists 25 candidates with membership annotations.
- Links: https://arxiv.org/abs/1504.07916

## Skipped By Rule 5

- Canes Venatici I: diff1=0, diff12=0; no further check this round.
- Canes Venatici II: diff1=0, diff12=0; no further check this round.
- Coma Berenices: diff1=0, diff12=0; no further check this round.
- Eridanus II: diff1=0, diff12=0; no further check this round.
- Horologium I: diff1=0, diff12=0; no further check this round.
- Hydra II: diff1=0, diff12=0; no further check this round.
- Leo T: diff1=0, diff12=0; no further check this round.
- Pisces II: diff1=0, diff12=0; no further check this round.
- Reticulum II: diff1=-7, diff12=-7; follow-up result recorded above.
- Ursa Major I: diff1=0, diff12=0; no further check this round.
- Ursa Major II: diff1=0, diff12=0; no further check this round.
- Willman 1: diff1=0, diff12=0; no further check this round.
