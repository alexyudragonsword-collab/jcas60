# 变更记录 / Changelog

按里程碑记录。每条写**做了什么**和**为什么**；数字原值在 `cairn/LOG.md`，结论在 `cairn/comparison-table.md`。

## 0.1.0 —— 里程碑 0：一张对照表（2026-11-24）

开库。mmw60_trx 的边界是通信链路，fmcw_radar 的边界是 FMCW 雷达，“用通信波形做感知”两边都不收；本库把两半 vendor 进来，回答同一套 60 GHz 前端上两种做法各自能做到什么、代价是什么。

- **对照表由代码生成**（`jcas60.table`，`examples/00_comparison_table.py` 打印，期望输出入库并有漂移检查）：5 m 处的人体，三种做法五行——A（fmcw_radar 默认配置，报告原样）、B（ay 单载波 L-CEF 单站，按通信限值 / 按雷达限值各一行）、C（双站，同样两行）。每个数带出处标签；`[derived/assumed]` 表示推导里用了至少一个假设。
- **两条路**：闭式（`budget.py`：雷达方程 + 处理增益）与波形级（`wavesim.py`：vendored 的 Golay 序列拼成 L-CEF，成形、点目标信道、相关接收、距离-多普勒图）。四个 ay 行差 −0.05 … −0.14 dB。单站的三项代价（带内压缩、ADC 量化、泄漏上的本振相噪）里后两项两条路都有。
- **共用的前端数字**（`frontend.py`，一个实例）取 fmcw_radar 默认配置用的那一组，因为 A 行不重算；A 行报告的 9 个输入与它逐项相等。本振是唯一的例外，表里标明。
- **法规两种读法都是行**：通信限值 (c)(1) 与雷达限值 (c)(2)；雷达读法直接跑 fmcw_radar 的检查（vendored）。哪一种适用在 spec_sources#54 里问了。
- **vendor**：spec_sources@`f9ce3f2`（读取器 + 5 张表）、mmw60_trx@`cfb401a`（`ay/golay.py`、它的本地表、同台表的期望输出）、fmcw_radar@`fb0fc16`（检测、CFAR、窗、读图、合规检查、示例 00 的期望输出）；15 个文件，13 verbatim、2 只改 import 行。A 行的报告由 `tools/upstream_report.py` 从 pin 的提交重新生成后核对。
- **没有按原计划 vendor 的**：`ay/training.py`（是 OFDM 速率的 EDMG-STF / EDMG-CEF，且拷不出来；L-CEF 由本库照 spec_sources 的两条拼，mmw60_trx#89）；`Frontend60`（里程碑 0 的边界：同一组数字，不是同一个对象）。
- **工具与 CI**：`vendor_check.py`（从 fmcw_radar 的改过来，多核 `.txt`）、`upstream_report.py`、`fingerprint.py`（504 个值，记录入库）、`example_drift.py`；CI 跑 ruff、四个 Python 上的 pytest、把三个上游检出到 pin 上的 vendor 核对、base 对 head 的指纹。
- 七条铁律在 `AGENTS.md`，各有测试（30 个）；Project Cairn 的文档齐。
