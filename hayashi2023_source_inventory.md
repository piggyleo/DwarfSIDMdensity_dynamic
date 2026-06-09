# Hayashi et al. (2023) 数据来源定位清单

这份清单服务于项目第 2 步的第一阶段：先为全部 27 个目标星系定位论文使用的数据来源，再进入实际的数据抓取。

配套文件：

- [hayashi2023_reference_catalog.csv](</Users/wangkaihao/Documents/New project 4/hayashi2023_reference_catalog.csv>)
- [hayashi2023_galaxy_data_sources.csv](</Users/wangkaihao/Documents/New project 4/hayashi2023_galaxy_data_sources.csv>)
- [hayashi2023_table1_observables.csv](</Users/wangkaihao/Documents/New project 4/hayashi2023_table1_observables.csv>)
- [hayashi2023_reference_vizier_availability.csv](</Users/wangkaihao/Documents/New project 4/hayashi2023_reference_vizier_availability.csv>)
- [data/raw/vizier/manifest.csv](</Users/wangkaihao/Documents/New project 4/data/raw/vizier/manifest.csv>)
- [data/raw/vizier/galaxy_local_source_map.csv](</Users/wangkaihao/Documents/New project 4/data/raw/vizier/galaxy_local_source_map.csv>)
- [data/processed/kinematics_ref3_jenkins2021_bootes1_leo4_leo5.csv](</Users/wangkaihao/Documents/New project 4/data/processed/kinematics_ref3_jenkins2021_bootes1_leo4_leo5.csv>)
- [data/processed/kinematics_ref9_ji2021_antlia2_crater2.csv](</Users/wangkaihao/Documents/New project 4/data/processed/kinematics_ref9_ji2021_antlia2_crater2.csv>)
- [data/processed/star_observations_master.csv](</Users/wangkaihao/Documents/New project 4/data/processed/star_observations_master.csv>)
- [data/galaxies/_coverage_summary.csv](</Users/wangkaihao/Documents/New project 4/data/galaxies/_coverage_summary.csv>)
- [data/galaxies/_field_completeness_summary.csv](</Users/wangkaihao/Documents/New project 4/data/galaxies/_field_completeness_summary.csv>)
- [data/galaxies/_kinematic_velocity_gaps.csv](</Users/wangkaihao/Documents/New project 4/data/galaxies/_kinematic_velocity_gaps.csv>)
- [data/raw/arxiv/1506.01021_src/dwarfs2015.tex](</Users/wangkaihao/Documents/New project 4/data/raw/arxiv/1506.01021_src/dwarfs2015.tex>)
- [data/raw/arxiv/1007.3499_src/ms.tex](</Users/wangkaihao/Documents/New project 4/data/raw/arxiv/1007.3499_src/ms.tex>)

## 说明

- `hayashi2023_reference_catalog.csv` 是 Table 1 脚注中 `ref_id -> 文献` 的总表，并附上我对文献用途的初步分类。
- `hayashi2023_galaxy_data_sources.csv` 现在是严格的 `Table 1` 版本：每个星系只保留该星系在论文 Table 1 明确给出的 references。
- `status = table1-strict` 表示这行数据来源与论文 Table 1 引用完全一致，不引入额外替代文献。

## 当前结论

- 27 个目标星系都已经有了第一轮来源定位。
- 每个星系都可以直接通过 `table1_ref_ids` 回溯到论文脚注的 `References (1) ... (35)`。
- 我已撤回先前对部分星系的“样本量匹配推断来源”，避免与论文文本不一致。
- 后续抓取时将按这份严格清单执行：先从该星系 Table 1 所列文献入手，不额外换源。
- 我已把 Table 1 的观测量转存为结构化 CSV，误差统一拆成 `err_down`/`err_up` 两列（即使对称误差也拆分）。
- 我已完成 35 篇参考文献的 VizieR 可用性探测：13 篇可直接下载机读表，21 篇当前在 VizieR 未找到，1 篇连接异常待复核。
- 可下载机读表已落地到 `data/raw/vizier/`，并生成下载清单与星系映射清单。
- 已完成第一份逐星运动学规范化样例：`Jenkins et al. 2021` 的 Bootes I / Leo IV / Leo V 数据。
- 已完成第二份逐星运动学规范化样例：`Ji et al. 2021` 的 Antlia 2 / Crater 2 数据。
- 该样例中 `member_flag==1` 的数量与 Hayashi Table 1 的 `Nsample` 不完全一致，说明还需按原文 membership 规则进一步筛选，不能直接等同。
- 在 `Ji et al. 2021` 样例中，直接使用 `Mm>=0.95` 的成员定义时，样本量也与 Hayashi Table 1 的 `Nsample` 不完全一致，后续会按论文具体筛选口径复刻。
- 我已按论文目标顺序生成 `27` 个“每星系一个 CSV”文件，存放在 `data/galaxies/`。
- 每个星系文件都包含：Table 1 全局参数 + 当前已抓取的逐星数据（若暂无逐星数据则明确标记 `star_data_not_yet_collected`）。
- 当前逐星覆盖统计见 `_coverage_summary.csv`，可用于后续按星系逐个补齐。
- 已补齐此前为空的 `Hydra II`、`Pisces II`、`Willman 1`：
- `Hydra II` 与 `Pisces II` 使用 `Kirby et al. 2015 (ref 21)` 的目标星表（arXiv:1506.01021 Table `Target List`）。
- `Willman 1` 使用 `Willman et al. 2011 (ref 34)` 的候选成员星表（arXiv:1007.3499 Table 2），并保留文中 5 颗“likely non-member”标记。
- 覆盖汇总现已达到 `27/27` 星系都有逐星行（均为 `star_data_available_partial`）。
- 已完成字段完整性审计；当前主要缺口是若干星系 `v_los` 列为空（见 `_kinematic_velocity_gaps.csv`）。
- 我已验证 `Simon & Geha 2007` 的 arXiv 源码与 PDF 都不直接内嵌完整逐星速度长表；并尝试访问期刊补充材料页面，但遇到网站验证码拦截（自动化抓取受限）。

## 下一步

下一轮我会重点做两件事：

- 按 Hayashi 建模所需字段检查每个星系 CSV 的完整性（坐标、速度、速度误差、成员信息、金属丰度及误差）。
- 对“只有部分字段”的星系补抓对应来源，并把成员筛选口径固化为可复现脚本。
