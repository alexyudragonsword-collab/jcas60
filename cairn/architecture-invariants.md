---
type: project_topic
status: active
summary: "七条铁律各自的论证、实现位置与守护测试；里程碑 0 的边界（同一组数字而不是同一个对象；A 行不重算）。"
tags: [jcas60, invariants, tests, vendor]
contains: [decision, scope]
created: "2026-11-24"
updated: "2026-11-25"
related: [comparison-table.md]
authoring_mode: ai_generated
---

# 铁律 ↔ 实现位置 ↔ 守护测试

规矩照两个上游的写法定（mmw60_trx 与 fmcw_radar 的 `AGENTS.md`、`cairn/architecture-invariants.md`）：每条有测试，数字只按键读，vendor 是可校验的副本，指纹证明重构。

| # | 铁律 | 实现 | 守护测试 | 状态 |
|---|---|---|---|---|
| 1 | 外部数字只从 spec_sources 进来，带标签 | `sources.py`：vendored spec_sources 五张表 + 本库 `data/assumed.yaml`（`jcas60.*`，19 条，全 assumed）+ mmw60_trx 自己的表，用 spec_sources 的读取器合起来、同键报错；数以 `Num`（值、标签、引用）流动，`derive()` 标 derived 并带上输入里最弱的标签（`[derived/assumed]`）；每次读取进 `TRACE`，示例 §10 整份打印 | `test_laws.py::test_law1_every_reference_read_resolves_…`、`_the_three_tables_are_kept_apart_by_prefix`、`_own_entries_are_all_assumed_and_say_what_would_move`、`_grep_count_of_each_tag_equals_the_entries`、`_no_stray_numbers_in_the_code`（AST 扫本库代码里的数字字面量：`physconst.py` 之外只许一位有效数字的算术常数与列明的几个） | 守住 |
| 2 | 同一套前端数字 | `frontend.py::shared()` 是唯一读前端引用的地方，一个 frozen `FrontEnd` 实例经 `config.setup()` 交给预算、仿真、法规与表；PA 回退是一条规则 `max(0, P + G − 上限)` | `test_law2_one_front_end_instance_serves_every_row`、`_row_a_ran_on_the_same_numbers`（A 行报告的 9 个输入逐项相等；规则给出的回退 = 报告的 8 dB）、`_front_end_numbers_are_read_in_one_module_only`、`_the_lo_is_named_as_the_one_number_not_shared` | 守住；本振例外，表里标明 |
| 3 | 各按各的法规上限 | `regulatory.py`：通信读法读 `fcc15255.eirp_avg_dbm`，雷达读法跑 vendored 的 `link/compliance.py`；ay 每行两个版本 | `test_law3_every_row_states_its_rule_and_limit_…`、`_no_row_runs_at_another_rows_limit`、`_the_vendored_check_reproduces_fmcw_radars_own_report`（逐字的副本经 `check_emission`、限值由本库的 `fcc_number` 从 spec_sources 取，在 fmcw_radar 默认帧上逐档等于它的报告） | 守住 |
| 4 | 两条独立的路 | `budget.py`（雷达方程 + 处理增益；压缩 / ADC / 相噪三项判据）对 `wavesim.py`（vendored 的构造 → L-CEF 或 TRN 场 → 成形 → 点目标信道 → 接收机的限幅器与 ADC → 相关接收 → 距离-多普勒图，用 fmcw_radar 的读图函数读峰与底） | `test_law4_the_two_roads_agree_within_1_db_on_every_ay_row`（八行，现状 ≤ 0.15 dB）、`_the_difference_is_a_column_of_the_table`、`_the_receivers_cost_has_both_roads_where_the_closed_form_claims_to_hold` | 守住。压缩判据 2026-11-25 起也有两条路（闭式：fmcw_radar 的小信号理论；仿真：接收端的限幅器），第一轮只有闭式 |
| 5 | 数值指纹 | `tools/fingerprint.py`：表的每个数、两条路的中间量、两种发射场与四种仿真帧的图，1717 个值（对照表 993 个、示例 01 的研究 724 个；标量 repr、数组 SHA-256）；记录 `tools/fingerprint.json` 入库并带环境（换了环境标量比到 1e-6，所以舍入残差不当读数记：低于 `Simulator.rounding_dbc` 的电平按界记）；CI 的 fingerprint 作业并排记录 base 与 head | `tests/test_fingerprint.py` | 守住 |
| 6 | vendor 纪律 | `src/jcas60/vendor/`：18 个文件，17 verbatim、1 adapted（`processing/rd.py` 的一行 import）；`tools/vendor_check.py --strict`；A 行的报告不是上游仓库里的文件，由 `tools/upstream_report.py` 从 pin 的提交（`git archive` 导出到临时目录）重新生成后核对 | `tests/test_vendor.py`（完整 SHA、三处写的 pin 一致、strict 通过、只拷用到的）、`test_table.py::test_the_recorded_report_is_what_the_pinned_fmcw_radar_writes`；CI 的 vendor-drift 作业把三个上游检出到 pin 上真跑 | 守住 |
| 7 | 措辞 | `jcas60.WORDING`，示例第 2 行；`bf.py` 把感知帧逐项对 spec_sources 的 802.11bf 条目（示例 §11），每项给出条目、本库的做法与“对得上 / 对不上 / 条目没有定 / 没有条目” | `test_law7_no_conformance_claim_anywhere`（全仓库没有一句一致性声称）、`_what_is_said_is_checked_item_by_item_against_the_entries`（“对得上”是由条目算出来的；放不进字段的帧会被报出来） | 守住。2026-11-25 之前的写法是“有条目之前不说”：条目到了（spec_sources#54），措辞没有放宽，多了逐项核对 |

## 里程碑 0 的边界

- **同一组数字，不是同一个对象。** 不 vendor mmw60_trx 的 `Frontend60`：它站在整条 wifitrx 链上（依赖私有的 wifi_trx_cal），fmcw_radar 0.2.8 做过同一个决定。代价：不是“同一个实例”跨库，而是本库里一个实例、数与出处同一份。
- **A 行不重算。** 表里 A 行的每个格子要么是报告里的一个值（`test_table.py::ROW_A` 逐项对），要么由报告里的数推出并标 derived（Pd 0.9 的距离、处理增益、空口占比）。
- **拷不出来的不硬拷，请上游拆。** `ay/training.py` 带进 mmw60 的垫片、它的本地表和第二份 spec_sources，没有拷（也用不上：它是 OFDM 速率的场）；单载波的训练场请 mmw60_trx 拆成了只依赖 numpy 的 `ay/preamble.py`（mmw60_trx#89）。`link/compliance.py` 原来的两行 import 指向 fmcw_radar 的配置层与来源表，第一版只改 import 行、给 `check()` 喂一个只带那几个属性的对象；fmcw_radar#42 之后它加载时不导入自己库的任何东西，取限值的函数是参数，现在逐字 vendor、调 `check_emission`。A 行“Pd 0.9 的距离”也从本库的外推改成读报告的 `detection.range_at_pd`。

## 分层

```
sources (三张表, Num, TRACE)        bf (感知帧 ↔ 802.11bf 条目, 只报告不出数)
   └─ frontend.shared() ── config.setup() ──┬─ budget   (闭式)
field.l_cef() (vendor: preamble, golay) ────┤─ wavesim  (波形级; vendor: rd, cfar)
fmcw_row (报告 JSON, 示例 00 文本) ─────────┤─ regulatory (vendor: compliance)
                                            └─ table.build() → render() → examples/00
studies (表之外的研究: 同一个 Setup, 每项只改一样, 两条路) ── build() → render() → examples/01
```
