# 从实验说明到替换证据总表

1. 读取用户口述、项目说明和实验清单。需要定位材料时运行 `python scripts/inspect_inputs.py 项目目录`。打开命中的说明或附件，照 [参数选择](design.md) 生成 design.json。可单独生成两轮 XML 检查。
2. 查已有搜索和中间产物，确定本次从哪步开始。初始化新的运行目录。学校主机使用其实际路径及已安装的 MaxQuant 2.8.1.0、Python 3.11+、NumPy、pandas、pyteomics。读 mzML 时安装其解析依赖，包括 psims。参考和污染库按版本及校验值定位。
3. 运行 pipeline.py。它执行已选择的搜索或读入已有导出，生成逐扫描候选、竞争 FASTA、同谱竞争及最终表。严格候选及其来源参与第二轮建库。宽松表的所有目标替换假设进入最终表，未通过严格筛选的行保留原条件和警告。输入适配剔除的 decoy 或无效定位另在 dp_review.tsv。
4. 提供峰与化学状态。新项目用同轮基底 PSM、第二轮 PSM 和每个 RAW 对应的质心 mzML。已复核的历史项目也可提供完整峰的 JSON 摘录、冻结的比较单元与链接、已有 Koina 响应。竞争双方在同一扫描、电荷和窗口比较。每个候选保留最佳竞争者及其身份相同、竞争更强、无可比者或缺数据等状态。
5. 检查覆盖与结果声明。最终文件为 `results/SAAP_evidence.tsv`，UTF-8 BOM，可在 Excel 以制表符导入。交付前核对各行 DP_hypotheses 之和与宽松表行数。需要审阅很长的逐谱注释时用 TSV 读取器读取原文件，Excel 单元格有显示长度限制。

## 运行

在同一工作目录下准备 workflow.json。例如：

```json
{
  "design":"design.json",
  "reference":"/data/ref/search.fasta",
  "isoforms":"/data/ref/expanded.fasta.gz",
  "contaminants":"/opt/MaxQuant/bin/conf/contaminants.fasta",
  "maxquant_command":["dotnet","/opt/MaxQuant/bin/MaxQuantCmd.dll"],
  "evidence":{
    "mzml":{"P01":"/data/peaks/P01.mzML"},
    "model":"Prosit_2024_intensity_PTMs_gl",
    "nce_grid":[24,28,32,36]
  }
}
```

- 仅生成两轮参数：`python scripts/pipeline.py workflow.json --parameters-only --out planning-01`。
- 接现有导出：`python scripts/pipeline.py workflow.json --dp-txt /data/r01/combined/txt --closed-txt /data/r03/combined/txt --out replay-01`。检查实际 XML 的 RAW、酶切与化学设置。第二轮历史 FASTA 若与本次库不同，其匹配状态单独写入表中。
- 实际检索：用户要求启动搜索时使用 `--execute-searches`；已有任一轮可传对应 txt。只在新的专用搜索目录建立原 RAW 同名链接。搜索日志、两轮 XML、返回码、阶段 job 和哈希保存在本次运行目录。搜索已完成而下游失败时，下一次直接接搜索导出。
- Prosit 新推断：在 evidence 明确填写学校 Koina 地址或获准使用的端点 `endpoint`，运行时加 `--predict`。有可信校准时提供 `calibrated_nce` 和 `calibration_source`。否则程序从普通 HCD PSM 选对照、剔除候选及基底序列、按 I/L 序列分出留出组，在同一可评分集合上按中位谱角选择能量。选择后用相同模型与能量预测候选及竞争解释。远端仅需模型输入，不传实测峰。缓存按序列、化学形式、电荷、能量和模型复用，支持项目既有 rows/response 格式。没有网络授权时在本地生成请求或使用缓存，将缺口如实写入表。
- 只重做评分或汇总：复制上次 `evidence_inputs.json`，更新输入，运行 `python scripts/saap_table.py evidence_inputs.json --out review-02`。先导出缺失请求可加 `--requests-only`；学校可用 `prosit_batch.py prediction_requests.tsv --model 模型名 --nce 已选能量 --cache 缓存目录` 生成原生模型请求，提供 `--endpoint` 时才提交。

### 竞争范围与输出口径

每行由基底肽、1 基位点及候选肽共同定义；同一替换的不同 RAW、扫描、来源蛋白和化学形式保留在该行。I/L 合并事件另有单独列。不同位点和不同基底肽各自成行。范围由本次 DP 筛选的质量窗口、参考、修饰、酶切和可用实测谱限定。

新项目的常规比较包含可恢复的基底肽修饰态、同扫描其他单替换位置、基底单 PTM、指定参考和污染库的固定化学酶切肽，以及已导出的闭合搜索 PSM。修饰破坏原位点或基底 PSM 未知时保留化学状态待复核。旧项目的扩展竞争单元含更多 PTM 组合与半特异解释；有这些资料时在 evidence 指定 `units`、`links`，按原比较单元的前体、质量与上下文资格重放，避免把不同化学状态的竞争者混在一起。可按模型允许的标记与电荷生成 Prosit 预测。强度评分覆盖 b/y 离子及最高三价碎片，默认 20 ppm，一对一配峰。实测竞争按解释强度比例排序，Prosit 竞争按平方根强度的归一化谱角排序。两类最优者可不同。

表中的直读列给出假设/谱图数、严格通过数、最佳 DP PEP、位点概率、可评分数、代表扫描的两种最优竞争解释及分差。代表扫描优先有完整竞争评分，其后按 RAW 与扫描排序。结果依据全体扫描，不能只引用代表扫描。详列包含 DP 注释与簇信息、闭合 PSM/PIF、观测与预测竞争者、双方分数、可区分的位点离子、尚未获得的预测和普通肽同化学身份歧义。原始强度差或谱角差是相对证据。

同一批次的既有附加证据可按需导入：`prior_unit_evidence` 接 primary_unit_evidence.tsv，附加其他扫描、RT 残差及诊断离子；`class_q` 接 candidate_134_FDR_results.tsv 格式的序列层条件 q；`paired_events` 与 `abundance` 成对接事件及 XIC/配对丰度。没有对应观测时保留不可判定。q 的样本、库、层级和打分策略必须与来源匹配，端到端新估计需要对应 target/decoy 审计；自发现同数据重搜有选择偏倚。TMT 的报告离子若存在会随闭合 PSM 保留，分析供体差异时应结合通道表、干扰、杂质校正及跨批归一化。

