# 单独运行候选发现和竞争建库

以下命令从 skill 根目录执行。本机验收环境为 Python 3.11.16、pandas 3.0.6、NumPy 2.4.6。脚本要求 Python 3.11 或更新版本，冻结任务与校验返回文件只依赖标准库。研究软件与公共数据继续由项目清单管理。

## 从逐扫描 DP 生成候选

以实际输入的公共父目录替换 `DATA`，其余输入文件参数均相对于该目录。将作业清单保存在已有配置目录，将输出指向新的运行目录。命令在 PowerShell 和 Linux shell 中均可按单行执行。

```text
python scripts/saap.py freeze --stage discover --root DATA --scans r01/combined/txt/msmsScans.txt --mqpar r01/combined/txt/mqpar.xml --job discover.json
python scripts/saap.py run --job discover.json --root DATA --out DATA/derived/discovery-01
python scripts/saap.py verify --job discover.json --out DATA/derived/discovery-01
```

使用同轮完整搜索参数与导出表。只有聚合 `dependentPeptides.txt` 时先取得逐扫描导出。当前配置按同 RAW 的 DP 基底做适配，跨 RAW 基底记录进入 `dp_review.tsv`。脚本从实际 XML 读取固定 CAM(C) 或游离 C，及无标记或 TMT6。五例基准的 C 有效残基质量为 160.030649 Da。

继承五例的初筛条件为质量差误差绝对值 ≤0.005 Da、DP PEP ≤0.01、最大和所提位点定位概率均 ≥0.90、肽长差为零、簇内 forward 数大于 reverse。历史 PTM/制样警告规则参与严格筛选。前三类数值可分别用 `--mass-tolerance-da`、`--pep-max`、`--localization-min` 在 freeze 时配置。警告词表的覆盖有限，后续仍需化学竞争。

`candidates.tsv` 是进入下一阶段的表，保留物理谱键和 I/L 合并事件 ID。`saap_candidates_relaxed.tsv` 保留按质量枚举的假设及其严格筛选标志。`dp_review.tsv` 记录 DP 输入适配时排除的行，`completion.json` 汇总各层计数及文件哈希。短肽文件 `saap_candidates_strict.fasta` 供核对候选序列。首阶段只需两份导出文件，PIF、基底 PSM、峰图和普通肽竞争均另行复核。

## 从候选生成竞争库

```text
python scripts/saap.py freeze --stage competition --root DATA --mqpar r01/combined/txt/mqpar.xml --candidates derived/discovery-01/candidates.tsv --reference refs/search.fasta --isoforms refs/expanded.fasta.gz --contaminants refs/contaminants.fasta --job competition.json
python scripts/saap.py run --job competition.json --root DATA --out DATA/derived/competition-01
python scripts/saap.py verify --job competition.json --out DATA/derived/competition-01
```

提供实际版本的人类搜索 FASTA、包含同工型的扩展参考及 MaxQuant 污染库。提供 mqpar 时使用其中的酶、漏切、最小肽长和最大肽质量；历史调用省略 mqpar 时沿用五例的 Trypsin/P、两处漏切、最短 7 aa、5000 Da。起始 Met 切除须另行复核。内存需求随参考和候选规模增长。

- `competition.fasta` 含完整正常参考及所需单替换母蛋白，按完整蛋白序列去重。第二轮用独立参数、运行目录，关闭 DP，将污染物及 decoy 交由搜索引擎按已核实设置加入。
- `candidate_context_mapping.tsv` 将候选和原始谱图连到真实母蛋白位点，保留正常竞争 ID、候选蛋白 ID、切点和 unresolved_parent。上下文定位依赖 DP 源蛋白登录号，来源未解应复核参考版本。
- `normal_competitors.tsv` 保留参考和污染库中的精确与 I/L 等价肽命中。`reference_aliases.tsv` 保留完整参考去重前的原始名称；污染库条目随引擎加载。
- `completion.json` 证明该阶段文件已生成并记录校验值。第二轮检索结果、候选专属错误率以及生物学来源另有各自证据。
