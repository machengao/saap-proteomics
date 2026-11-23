# 将项目说明和口述实验设计写入 mqpar

收到 RAW 文件夹或口述方法时，由当前助手完成设计解析。优先查看相邻 README、实验设计表、原作者搜索参数和制样/采集章节。已有 mzML、导出扫描表或仪器方法可给出仪器、激活方式、质量分析器、隔离窗口、MS 层级和实际采集方式。核对目录内实际项目及 RAW 文件名，区分同一文件夹中的不同队列。按来源、原文、字段值和决策理由填写 evidence_sources 与每组 basis。向用户展示会改变结果的冲突。对当前资料已经支持的字段直接选择，关键化学或采集方式无法推断时仅询问缺失项，先保留草案。

| 实验线索 | 选择与核对 |
|---|---|
| DDA、DIA、HCD、CID、ETD、FTMS 或 ITMS | 分别记录采集策略、激活方式、检测器与仪器家族。优先采用实际文件头与设计的交叉证据。HCD 本身不能判定 DDA。常规 DDA 选 DP，DIA 转入相应分析设计；ETD 等保留未实现的碎片评分状态。混合仪器依实际参数组组织。 |
| 仅需鉴定、未标记强度、明确要求 MaxLFQ | 单通道 Standard，可按设计选 unlabelled 或 MaxLFQ；后者写入 lfqMode=1 与最小比值数。疾病名称、组织类型和物种均不足以推断定量方法。 |
| TMT6，MS2 reporter 或 SPS-MS3 | 写入 Reporter MS2 或 Reporter MS3，六组 N 端/K 标签、对应批次的四项杂质系数，MS3 时导出 MS3 表。文献只写 TMT 时，继续从方法或文件确认 plex、读出层级及通道。明确记录系数来源；只在项目允许无校正时使用四个零。混样谱用于鉴定，通道到供体的定量映射单独核查。其他 plex、SILAC 等可用同版同类型原生模板生成参数，当前下游可执行的化学范围见下。 |
| 胰酶或其他酶 | 区分 Trypsin 与 Trypsin/P、LysC、混合酶、漏切数、全酶切或半特异。按协议和原作者参数决定。不同蛋白酶切点会改变替换的可检验性。 |
| 烷基化、富集、其他制样条件 | 根据实际试剂选固定修饰。按研究目的选择可变修饰和准确残基/端基；有脱酰胺证据时纳入与 N→D、Q→E 的竞争。TMT 通过标签配置编码一次。酶与修饰名称用本机 MaxQuant 安装的 conf 目录核对。 |
| 污染物、decoy、样本/分级 | 明确内置污染物开关与自建库是否已含污染物，保留 decoy 规则。将 RAW 对到实验、分级、参数组、TMT 参考通道；用独立运行目录承载两次搜索。 |

### design.json

下例为口述“Q Exactive Plus，人源 DDA/HCD，未标记、需要 MaxLFQ，Trypsin/P 两处漏切，碘乙酰胺烷基化，使用常见污染库”的参数结构。填写实际路径和实验依据。项目说明不支持的修饰应从列表调整。相对路径均从启动命令时的工作目录解析。对远端主机应使用其文件系统路径。

```json
{
  "evidence_sources": ["实验方法第 2 节及文件头；样本清单第 1 行"],
  "unresolved": [],
  "include_contaminants": true,
  "threads": 8,
  "maxquant_conf": "/opt/MaxQuant/bin/conf",
  "dp_fasta": "/data/ref/human.fasta",
  "closed_fasta": "/data/project/run02/database/competition.fasta",
  "files": [{"raw_source":"/data/raw/P01.raw","dp_path":"/data/project/r01/P01.raw","closed_path":"/data/project/r03/P01.raw","experiment":"P01","group":0}],
  "groups": [{
    "instrument":"Q Exactive Plus", "instrument_family":"Orbitrap", "acquisition":"DDA", "fragmentation":"HCD",
    "quantification":"MaxLFQ", "enzymes":["Trypsin/P"], "enzyme_mode":"Specific", "missed_cleavages":2,
    "fixed_modifications":["Carbamidomethyl (C)"], "variable_modifications":["Oxidation (M)","Acetyl (Protein N-term)"],
    "basis":"制样记录及用户要求 MaxLFQ。氧化和蛋白 N 端乙酰化沿用已核实的搜索设置。"
  }]
}
```

运行 `python scripts/mq_design.py design.json --out parameters-01`。得到首轮 DP、第二轮闭合竞争的两份 XML 及参数决策记录。随包模板来自已实跑的 MaxQuant 2.8.1.0 Orbitrap DDA 参数，已移除旧数据路径、样本编号。其他仪器用组字段 `template` 指向该版本相匹配的原生 XML；按 `template_group` 选择组。检测器、碎裂容差和激活类型自动识别设置沿用核对后的原生模板。混合组须有一致的全局检测器配置。TMT6 时改用 `TMT6-MS2` 或 `TMT6-MS3`，以 `channel_corrections` 映射 126–131 各通道的 [M−2,M−1,M+1,M+2]。特定实验可用 `reporter_tolerance` 写入已核实的原生报告离子容差。

完整计算当前支持同一化学配置的人类 DDA、无标记或完整 TMT6、固定 CAM(C) 或无固定修饰、已登记的全特异酶切。蛋白末端参与酶切判断；起始 Met 移除等须按引擎另行复核。TMT 的编码按完整标记处理，K 获得/丢失会改变质量差。正式 TMT 数据先以普通标记肽和少量已知谱验收标记效率、DP 质量差、电荷与模型适用性。全流程统计分析及供体归一化应按具体队列设计执行。

原生参数组及命令行说明见 [MaxQuant 实验设计](https://cox-labs.github.io/coxdocs/rawFilesTab.html) 和 [MaxQuant 安装与参数配置](https://cox-labs.github.io/coxdocs/Download_Installation.html)。

