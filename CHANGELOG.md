# 变更记录 / Changelog

按里程碑记录。每条写**做了什么**和**为什么**；数字原值在 `cairn/LOG.md`，结论在 `cairn/comparison-table.md`。

## 未发布

### 用 TRN 场做感知；压缩判据有了第二条路（2026-11-25，用户定“都做”）
- **对照表多了四行**：B、C 在 EDMG TRN-R 场上感知（同一个最短的包后面接 4 个单元 = 30 720 码片，每个子字段相关第 2–5 段，每包 20 480 码片），各出通信限值与雷达限值两个版本。通信限值下 SNR 43.75 dB，比 A 行的 43.34 dB 高 0.41 dB（L-CEF 行是 30.74 dB）；雷达限值下 41.75 dB。空口时间 4.51 %（L-CEF 1.16 %，A 行 19.2 %）。两条路差 −0.02 … −0.05 dB。问题 1 的结论因此改写：SNR 上追得平，输的只有分辨率。
- **接收机的压缩进了波形级仿真**（fmcw_radar 的 `link/compression.py`，vendored；闭式用它的小信号理论）：灵敏度压低 1 dB 的隔离是 43.3 dB（两条路都是），不是原来留 10 dB 余量的 56 dB；余量这条 assumed 删掉，换成限幅器的膝部 p（`jcas60.rx_rapp_p`）。对查到的 40 dB 差 3.3 dB（原 16 dB）。新发现一项灵敏度之外的代价：压缩在目标的镜像多普勒上造一个假目标，40 dB 隔离时 −11.4 dBc，要压在噪声底下得 48.8 dB。问题 2 的结论改写。
- **读峰修了一处**：静止径旁边几个多普勒格的慢目标（双站行在第 4.47 格）以前被“绕开零多普勒脊”的排除带挡住，读到相邻的格，低 0.32 dB。双站 L-CEF 行两条路的差从 −0.14 / −0.10 变成 +0.12 / +0.15 dB；单站行逐位不变。
- **预期先于代码提交**（`91549b5`）：TRN 的六条都对上；压缩的三条里“约 46 dB”方向对、偏保守，“像在深压缩之前很弱”错了。
- vendor 17 个文件（16 verbatim、1 adapted）；感知帧对 802.11bf 条目的核对从 10 项变 13 项（8 对得上、4 条目没有定、1 没有条目）；spec_sources#57 请它补 TRN 场的三处细节与双站的 PPDU。

### spec_sources 答完后的重新 pin（2026-11-25）
- **spec_sources → `0d17de6`，spec_sources#54 验收**：只换 pin 时指纹 507 个值逐位不变，合并表 89 → 96 条（802.11bf 在 60 GHz 的 5 条、DMG / EDMG 的 TRN 场 2 条）。
- **感知帧逐项对 802.11bf 的条目**（`jcas60.bf`，示例新增 §11）：十项里 6 项对得上（单站的定义、单站可用任何 DMG / EDMG PPDU、双站的收发关系、包间隔 / 帧间隔 / 每帧包数放得进字段），3 项条目没有定（要不要边发边收、信道估计用哪一段、带宽），1 项没有条目（双站用什么 PPDU：C 行的普通 SC 包不能说对上）。铁律 7 的措辞没有放宽：仍是“用 ay 波形做感知”，多了“能说的只有逐项对过的”。
- **氧气衰减改读信道中心 60.48 GHz 的点**（14.993 dB/km，原来借 60 GHz 的 14.651；数字动了）：ay 四行的 SNR −0.0034 dB（30.75 → 30.74、28.75 → 28.74 dB），Pd 0.9 的距离短 7 mm 左右；`jcas60.oxygen_frequency_ghz` 这条 assumed 删掉（19 条）。A 行不动（它的中心是 60 GHz）。

### 两个上游答复后的重新 pin（2026-11-25）
- **fmcw_radar → `c303566`（0.2.33），fmcw_radar#42 验收**：`link/compliance.py` 逐字 vendor（它加了不读配置的 `check_emission`，加载时不导入自己库的东西），本库改的两行 import 与喂给 `check()` 的那个对象删掉；vendor 16 个文件，15 verbatim、1 adapted。它算的 ay 帧与本库原来的结果逐项一致（20 dBm 档、关断 26.70 ms、最短帧周期 14.553 ms）。
- **mmw60_trx → `442fabd`，mmw60_trx#89 验收**：L-CEF 改由它新出的 `ay/preamble.py` 建，本库自己拼的那份删掉。码片的 π/2 旋转从 exp 换成精确的 1 / j / −1 / −j：指纹 504 个值里 64 个在末位动（最大 3.2e-10 dB），静止时 Golay 的数值零从 −281.5 到 −293.0 dBc；其余打印出来的数不变。
- **Pd 0.9 的距离改读报告、三行同一个门限**（数字动了）：A 行读报告新增的 `detection.range_at_pd`（28.78 m，原来是本库外推的 29.56 m）；三行的门限从噪声电平已知的 11.75 dB 改成报告的 12.24 dB（图上 CA-CFAR 要的），B 14.68 → 14.28 m、12.76 m（原 13.12 m），C 相同。
- **spec_sources#54 答了法规一半，没有定论**：两种读法继续都是行。它的 pin 不动（`f9ce3f2`），等 802.11bf 的条目合并后一次重新 pin。

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
