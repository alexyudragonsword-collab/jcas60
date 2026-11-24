# jcas60 协作规则

> 本项目使用 Project Cairn：`AGENTS.md` 是规则与导航入口，`cairn/` 是项目知识/状态层。同目录 `CLAUDE.md` 只有一行 `@AGENTS.md`，使 Claude Code 读到同一份规则；Codex 直接读本文件。

## 项目一句话

60 GHz 通感一体：在同一套 60 GHz 前端数字上，把“用通信波形（802.11ay 的 Golay 训练场）做感知”和专用的 FMCW 雷达放进一张对照表——各自能做到什么、代价是什么。mmw60_trx 的边界是通信链路，fmcw_radar 的边界是 FMCW 雷达，用通信波形做感知两边都不收，所以单开一库；两半都从它们 vendor，外部的数从 spec_sources 进来。

## Init 配置

- 毕业 provider：暂缓对接（首次毕业时再接知识库）
- 知识库索引 / 毕业目标：尚未配置

## 进入项目后的阅读顺序

1. 先读本文件。
2. 读根目录 `ROADMAP.md`：只记没做的，每条附证据（复现方式与数字），含“报回上游”节。
3. 读 `cairn/LOG.md` 最新几条（新条目在顶部）：每轮的数字原值与被推翻的判断。
4. 按任务读 `cairn/comparison-table.md`（三个问题的结论，先结论后前提）、`cairn/architecture-invariants.md`（铁律 ↔ 实现 ↔ 测试）。

## 文档职责

| 文件 | 角色 | 维护方式 |
|---|---|---|
| `AGENTS.md` / `CLAUDE.md`（根） | 规则与导航 / 单行存根 | 很少变，≤ 65 行 |
| `ROADMAP.md`（根） | 待办与优先级 | 只写没做的；完成即删并在 LOG 记一笔；发现的问题带复现方式和数字 |
| `CHANGELOG.md`（根） | 已交付的演进 | 按里程碑追加 |
| `cairn/LOG.md` | 时间序日志 | 新条目加在顶部，每条 ≤ 20 行：数字原值、被推翻的判断、下一步 |
| `cairn/<topic>.md` | 知识专题文档（当前真相） | 就地更新；修正过去的判断时追加更正说明 |
| `cairn/Reference/` | 外部原始输入 | 按需创建，只追加 |

> 工程资产（来源表 `src/jcas60/data/`、vendored 的表、期望输出 `examples/expected/`、指纹记录 `tools/fingerprint.json`）留在代码树里；`cairn/` 只放关于它们的知识。

## 项目铁律（每条都有测试守着；论证见 `cairn/architecture-invariants.md`）

1. **外部数字只从 spec_sources 进来，带出处标签**（spec / survey / assumed / derived）。本库自己的取舍在 `data/assumed.yaml`（`jcas60.*`），全部标 assumed，note 写明“换了会动什么”；代码里不写物理数字，只按键读（`tests/test_laws.py::test_law1_*`）。
2. **两种感知方式用同一套前端数字**：发射功率、PA 回退（一条规则：回退到该行的 EIRP 上限）、接收噪声系数、天线增益、带内压缩点、收发隔离、本振相噪，只在 `frontend.py` 一处配置，一个实例交给每一行；A 行的输入与它逐项相等。唯一不能相同的是本振（FMCW 的本振就是出 chirp 的 ADPLL），表里标明（`test_law2_*`）。
3. **各自按自己适用的法规上限比**，两个上限都在表里写明；不许拿一方的满功率比另一方的合规功率。通信设备顺带感知适用哪一款没有定论之前，ay 的每一行出两个版本（`test_law3_*`）。
4. **一个结论至少有两条独立的路**：闭式预算（`budget.py`）与波形级仿真（`wavesim.py`），差多少进表；差超过 1 dB 先查原因、写进 LOG，不调参数凑（`test_law4_*`）。只有一条路的结论在表里写明“只有闭式”。
5. **数值指纹**：重构用 `tools/fingerprint.py` 证明逐位不变；记录入库（`tools/fingerprint.json`），CI 比较 base 与 head，`(numbers-move: 为什么)` 是逃生口（`tests/test_fingerprint.py`）。
6. **vendor 纪律**：`src/jcas60/vendor/` 是可校验的副本——pin 完整 SHA，`tools/vendor_check.py --strict` 全过，只拷用到的文件及其依赖，正文不改（只许改 import 行）；A 行的报告由 `tools/upstream_report.py` 从 pin 的提交重新生成后核对（`tests/test_vendor.py`、`tests/test_table.py`）。
7. **措辞**：在 spec_sources 有 802.11bf 的条目之前，本库的说法是“用 ay 波形做感知”，不作与 802.11bf 一致性有关的声称（`test_law7_*`）。

## 本仓库特有的硬规矩

- **边界**：里程碑 0 不 vendor mmw60_trx 的 `Frontend60`。“同一套前端”指同一组数字（从 spec_sources 和 mmw60_trx 的预设条目读），不是同一个对象。A 行的数直接取 fmcw_radar 的报告，不重算；报告里没有的列由报告里的数推出并标 derived。
- **不改 mmw60_trx、fmcw_radar、spec_sources**：需要它们做的，在那个库开 issue（标题以“jcas60 请求：”开头），记进 ROADMAP 的“报回上游”；对方合并后本库重新 pin、验收、关 issue。IEEE 标准正文与 Xplore 的许可限制输入 AI：不读，也不请用户发来。
- **预期先写下来，再跑**：每轮在 LOG 里先记预期，量出来不一样就按量到的写，并说明差在哪。
- **假设要可见**：每个取舍是一条 assumed，示例末尾整表打印，承重的有灵敏度行。
- 在 `feat/*` 分支上开 PR 合 main；提交信息写为什么，带实测数字，至少两条 break-check；模型名不进提交、PR 和文档。
- 界面（图形界面、Android）不在本库当前范围内。

## 冲突仲裁

- 优先级：**代码与实测数字 > 知识专题文档 > LOG 历史**；规则层面的冲突由本文件裁定。

## 知识沉淀规则

- 每次实质性推进后，在 `cairn/LOG.md` 顶部加一条（摘要 + 指针，带数字）；结论沉淀进 `cairn/` 专题文档，先写结论再写前提。修正过去的判断时追加更正说明，不要静默覆盖。
- 每个里程碑结束：CHANGELOG 一条、LOG 一条、ROADMAP 删掉做完的条目并加上过程中发现的问题（带复现方式和数字）。
