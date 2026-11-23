# SAAP Proteomics

从人类 DDA 蛋白质组的 MaxQuant dependent peptide 搜索中枚举单氨基酸替换肽（SAAP），用实测 MS2 和 Prosit 预测谱比较竞争解释，并生成逐替换证据总表。

- [技能入口](SKILL.md)
- [实验设计与 MaxQuant 参数](references/design.md)
- [完整工作流](references/workflow.md)
- [小样本替换证据总表](06_validation/saap_skill_r02/SAAP_evidence.tsv)
- [离线验收与复现](06_validation/saap_skill_r02/验收与复现.txt)

小样本来自 PXD034575 的既有搜索导出，覆盖 18 条扫描记录、69 条质量相容假设和 37 行替换。复现使用抽取的峰与预测缓存，可跳过 MaxQuant 搜索。运行前请按复现说明将 MaxQuant 配置目录改为本机路径。
