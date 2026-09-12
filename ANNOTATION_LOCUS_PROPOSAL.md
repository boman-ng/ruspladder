# Ruspladder 多位置注释处理与性能优化建议

日期：2026-09-12。状态：原研究建议；实现和验收记录另见本分支运行报告。已有实测见 [超时根因报告](/home/wubw/data/ruspladder/p0-timeout-root-cause/TIMEOUT_ROOT_CAUSE.md)。方案原件在 `perf/p0-timeout-root-cause`，本实现分支为 `perf/p0-locus-lifecycle`。

## 1. 问题归属与依据

不能把本次超时概括为“NCBI 注释有错，SplAdder 没问题”。证据支持三层结论：

- 当前 GTF 有 714 个 `gene_id` 跨 seqname 或 strand 使用。文件名和 RefSeq 标记不足以证明它是未经转换的 NCBI 原始下载；工作流的参考 manifest 记录了 DVC 来源及文件摘要，没有记录原始 NCBI 下载地址、注释版本和转换命令。
- GTF 2.2 将 `gene_id` 定义为转录本基因组位点的全局唯一标识。当前跨坐标系的用法与这一约定不一致。不能将 GFF3 中允许重复的生物学 GeneID 直接等同于 GTF 的唯一位点 ID。[GTF 2.2 原始定义](https://mblab.wustl.edu/GTF22.html)，[UCSC 格式说明](https://genome.ucsc.edu/goldenPath/help/GTF.html)。
- SplAdder v3.1.1 按 `gene_id` 分组本身有格式约定依据，但没有检查分组中的坐标系，使用首次出现的染色体/链和全组坐标极值构造区间，导致无意义的大跨度图和 BAM 查询。Ruspladder 为兼容复现了这一行为。这是应在注释导入边界修正的建模/鲁棒性问题；现有证据不意味着剪接事件算法整体错误。

NCBI 明确区分 GFF3 的 feature `ID`、`Parent` 与生物学数据库交叉引用：feature ID 在文件内唯一，gene symbol、GeneID 和转录本 accession 可以出现在多个 feature 上。同一 assembly 可存在共享 GeneID 的不同基因部分或等位位置。NCBI 还说明部分 feature ID 随文件生成而变化，不能承诺跨注释版本稳定。[NCBI GFF3 官方说明](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/)。

已有实现提供了可直接采用的思路：

- Bioconductor GenomicFeatures 将不能表示为单一 seqname/strand 区间的基因保存在 `GRangesList`（`single.strand.genes.only=FALSE`）。其 `genes()` 先按 gene ID 分组，再求保留坐标系的 ranges。借鉴多位置表达，不采用默认排除这些基因的策略。[维护者源码](https://github.com/Bioconductor/GenomicFeatures/blob/RELEASE_3_23/R/transcripts.R#L224)，[接口定义](https://github.com/Bioconductor/GenomicFeatures/blob/RELEASE_3_23/man/transcripts.Rd#L81)。
- GENCODE 对 chrY PAR 区域的 gene/transcript ID 使用位置相关后缀来消除重复。这支持“保留基因关联、区分位置实例”的导出方式，但不意味着 Ruspladder 应将其专用命名规则套用于所有物种。[GENCODE 数据格式](https://www.gencodegenes.org/pages/data_format.html)。

## 2. 最小实现范围

原则：保存原始基因身份，以单一参考序列、单一链上的位置实例构建普通剪接图。无需创建新的剪接算法或通用注释框架。

1. **当前 GTF 导入**：在现有两遍解析中，以 `(原 gene_id, seqname, strand)` 确定位置实例；gene 行与 exon 行使用同一作用域。转录本挂在所属实例下，重复 transcript ID 不跨实例拼接。保留所有原 exon 坐标、方向与归属记录。
2. **现有 GFF3 导入**：按显式 feature `ID/Parent` 维护关系；不用 gene symbol 或 `Dbxref=GeneID` 替换 feature ID。GFF3 允许一个不连续 feature 使用多行同 ID，不能把每一行都生成为独立基因。[Sequence Ontology GFF3 规范](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)。
3. **沿用现有下游模型**：每个位置实例继续使用当前 `Gene` 和剪接图算法。当前 `Gene.name` 也参与排除、跨样本合并与缓存；只改解析 HashMap 键不够。最小方案为受影响实例生成确定性唯一名称，并输出名称到原 gene ID、seqname、strand、位置范围的映射表；不冲突的名称保留。导出规范化 GTF 时同时消歧跨实例重复的 transcript ID，并保留原 ID 映射。
4. **保持实例身份贯穿全流程**：多样本必须使用同一注释及消歧规则；合并、排序、报告和缓存都使用一致身份。原模式与位置模式的缓存路径或命名必须隔离，避免旧图被当作修正后的图复用。第一版不需要改造全部缓存结构。
5. **明确推断边界**：这个三元组足以处理本次跨染色体/链问题，不是所有 locus 的完整定义。同一 contig、同一链的不同位置，优先利用显式 feature/parent 信息；缺少信息时不凭任意距离阈值拆分真实长基因。明确标注的跨位置转录、trans-splicing 等超出现有单位置剪接图契约的输入，应给出具体记录及不支持原因，不悄悄当作普通转录本处理。

规范化不删除 alt/haplotype contig，不修改 BAM 的比对归属或多重比对过滤，也不把不同位置的计数自动相加为基因级结果。当前 BAM 缺少某个 contig 时沿用已定义的处理，并在注释映射审计中保留该位置。

该处理可以消除如 ZNF84 的 chr12 与 Un_gl000223 坐标混合，但不承诺恢复源文件可能已经丢失的生物学身份或 feature 层级。

## 3. 输出一致性的两个验收对象

| 路径 | 比较基准 | 要求 |
| --- | --- | --- |
| 原 SplAdder 兼容路径 | Python 使用同一原始输入 | 保留原有结果契约；调度、读取、内存优化必须等价 |
| 按位置处理路径 | Python 与 Rust 使用同一确定性规范化注释 | 比较完整图、事件坐标、计数、PSI、顺序及公共文件数据；Rust 原始输入加位置处理须与其规范化输入等价 |

正式默认路径暂不变；旁路验证位置处理，后续可通过单一显式注释策略选择启用。修复跨坐标系合并与保留旧错误图的输出不可能同时成立，不能将修正结果宣称为旧结果的无差异加速。

没有冲突的回归注释要求原输出一致；对本次有冲突的整份注释，不能承诺其余基因结果都不变，因为 SplAdder 的自适应过滤参数存在全局传播。对照必须检查整套输出，并记录最终过滤参数，不能仅比较 714 个冲突 ID。

Python 在规范化完整输入上的比较也受每次 1 小时硬限制（30 分钟软边界）；若仍超时，只能报告已完成阶段和小型用例的一致性，不能宣称完整三样本已经通过。

## 4. 性能实现顺序

每项单独对照，上一项有效后再决定下一项；不预先建设统一调度器、全 BAM 内存缓存或新存储格式。

| 顺序 | 现有证据 | 最小候选改动 |
| --- | --- | --- |
| A：注释位置归属 | 冲突关联查询占 S026 初始 BAM 查询时间 75.4% | 上述位置实例导入与身份映射，先验证科学语义 |
| B：覆盖度内存 | 单 BAM 仍有两份完整跨度覆盖度数组 | 消除多余副本；需要窗口处理时复用项目已有差分/窗口能力，验证原区间边界与 CIGAR 计数语义 |
| C：并行任务粒度 | 初始内含子读取串行；cassette 约 346.6 秒只有一个工作线程继续查询 | 复用已有旁路并行提取实现，评估 Rayon `with_max_len`；保持归并顺序，检查 reader 初始化次数及内存，不能假定 `map_init` 就是一线程一个 reader |
| D：重复解压 | 初始阶段采样约 59% 用户态周期在解压 | 优先小规模验证现有 HTSlib 有界缓存；仍受重复读取限制时，再评估多区域读取和按原查询分发 |
| E：缓存生命周期 | 位置拆分 S026 中图缓存读写累计 156.9 秒，占构建 44.6% | 先消除单样本流程中已在内存的图被重复读回；不先设计新的 HDF5 布局或并发 HDF5 写入 |

Rayon 已提供任务粒度控制：[1.12.0 `with_max_len`](https://docs.rs/rayon/1.12.0/rayon/iter/trait.IndexedParallelIterator.html#method.with_max_len)。小任务与读取器初始化/局部性有成本权衡，不预设 `max_len=1` 最优。

HTSlib 多区域迭代能减少重复读取，但同时会减少跨查询重复返回的记录。若用于 Ruspladder，必须将每条比对重新分发到原本所有应接收它的基因查询，保留原查询边界、过滤和计数，不能直接全局去重。[samtools/HTSlib 多区域读取说明](https://www.htslib.org/doc/samtools-view.html)。缓存命中率和收益尚未测得。

增加并行之前先控制覆盖度工作内存；否则多个旧模型超长基因同时分配数组可能超过 8 GiB。消除冗余计算可能降低 CPU 总时间，平均占用率应与总耗时、CPU 秒数和各阶段空闲时间一起判断，不以占满 64 核作为成功标准。

## 5. 已验证可行性与下一步验收

S026 的同一完整 BAM、同一诊断二进制，在 4 核/8 GiB 下：原注释 1,200.2 秒超时、峰值 RSS 2,137.5 MiB；按 `(gene_id, seqname, strand)` 拆分的诊断副本 351.8 秒完成、峰值 RSS 438.1 MiB。副本生成另耗时 28.3 秒。所有 1,490,438 行保留，仅改动冲突行的 gene ID。这是根因及方案方向的实证，不是完整科学验证，也不是正式优化加速比。

实现验收依次为：

1. 小型真实注释片段覆盖 ZNF84、chr6 haplotype、同位置多转录本、跨位置重复 transcript ID，以及正常长基因；独立验证坐标、归属、行保留和名称映射，再与 Python 规范化基线比较。
2. 跑既有无冲突端到端回归，验证模式切换、缓存重用和跨样本图合并不会串用身份。
3. 使用 S025/S026/S027 原 BAM，每样本 4 核、8 GiB、全流程 1 小时硬截止、30 分钟软边界，无可视化。记录注释准备、图构建、定量、事件输出、缓存各阶段，以及 CPU 秒数、平均用核、进程 RSS 和 cgroup 内存。首次准备与缓存复用分开报告，避免将输入准备时间隐藏在构建之外。
4. 比较完整最终输出，临时文件不算成功；每一项优化均保留基线差异报告。需要扩展性测试时再使用 8/16/32/64 核，总并发核数不超过 64；4 核基准仍是主要验收依据。

以上是待执行方案。本轮只补充规范与实现依据，不改变生产代码、原 GTF、参考库或现有样本结果。
