# Negative Member Count Omission Audit (2026-05-06)

检查对象：`n_member_flag_1_minus_hayashi2023 < 0` 且 `n_member_flag_12_minus_hayashi2023 < 0` 的星系。

## Summary

- Checked galaxies: 1
- Galaxies with confirmed omitted source-member rows: 0
- Total omitted source-member rows listed: 0
- Reticulum II has no row omission; it is a member-definition mismatch.

## Results

### Reticulum II

- Counts: Hayashi=25; current flag1=18; current flag2=0; current flag12=18; diff1=-7; diff12=-7.
- Source checked: Koposov et al. 2015b arXiv:1504.07916 Table 2; source member rows=18; represented source members=18; missing source members=0; current duplicate coord groups=0.
- Conclusion: 当前 CSV 已包含 Koposov et al. 2015b Table 2 的 25 个 Reticulum II 候选观测行；源文献明确识别 18 个 confirmed members，当前 member_flag=1 也为 18。因此未发现逐星数据遗漏；Hayashi n=25 更像是采用候选观测总数或不同成员口径，而不是我们漏了 7 行。
- Recommended action: 不要补新行；如后续要复现 Hayashi n=25，需要和用户确认是否将 7 个源文献未确认/非成员候选也纳入 Hayashi 等价动力学样本，或将 Hayashi n 视为候选数而非 confirmed member 数。
- Link: https://arxiv.org/abs/1504.07916

