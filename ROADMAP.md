# ROADMAP

**这份文件只写还没做的事。** 做完即删，并在 `CHANGELOG.md` 记一笔、`cairn/LOG.md` 记数字。里程碑 0（对照表）的结论在 `cairn/comparison-table.md`。

复现所有数字：`python examples/00_comparison_table.py`（约 30 s；§n 指它输出的第 n 节）。

## 报回上游（本库不改别的库；要它们做的在那个库开 issue，对方合并后本库重新 pin、验收、关 issue）

### spec_sources —— pin 在 `f9ce3f2`；[spec_sources#54](https://github.com/alexyudragonsword-collab/spec_sources/issues/54) 待它答复
- **802.11bf 在 60 GHz 的感知条目**：用哪些字段、感知包的结构与长度、重复方式、单站 / 双站 / 多站的定义。到之前本库的感知场（L-CEF）、感知包（最短的 SC 包）、一帧 64 个包都是 assumed，措辞守铁律 7。到了之后：按条目重建感知帧，B、C 两行重出，再议措辞。
- **法规：通信设备顺带感知适用 (c)(1) 还是 (c)(2)**，以及关断时间管不管通信发射。到之前 ay 每行两个版本（EIRP 差 2 dB），雷达读法的空口占比是 1.16 % … 50 % 的区间。到了之后删掉不适用的版本，或把区间收成一个数。
- 顺带：`itu_p676.gamma_o_dry_sea_level` 补 60.48 GHz 的点（现取 60 GHz 的 14.651 dB/km，`jcas60.oxygen_frequency_ghz`，差 < 0.01 dB）。
- 里程碑 0 合并后在 #54 回报本库的提交与 vendor 的文件清单（`reader.py` 与 5 张表：ieee80211、regulatory、frontend60、radar_targets、atmosphere），请它登记进 `consumers.yaml`。

### mmw60_trx —— pin 在 `cfb401a`；[mmw60_trx#89](https://github.com/alexyudragonsword-collab/mmw60_trx/issues/89) 待它答复
- **标准的 SC 训练场做成能单独拷的模块**。`ay/training.py` 是 OFDM 包的 EDMG-STF / EDMG-CEF，而且拷不出来（带进它的 `spec` 垫片、本地表与 pin 在 `65ff5e4` 的第二份 spec_sources）；单载波的 L-CEF 它没有按标准拼法建。本库现在自己拼（`field.py`）。它建了之后：vendor 它的、删掉 `field.py` 里的拼法，指纹应逐位不变（同样的码片）。
- 以后做 OFDM 感知时要 `ay/training.py` 本身：请它把读条目与建波形分开。

### fmcw_radar —— pin 在 `fb0fc16`；[fmcw_radar#42](https://github.com/alexyudragonsword-collab/fmcw_radar/issues/42) 待它答复
- **合规检查留一个不读 `RadarConfig` 的入口**。现在 vendor 的 `link/compliance.py` 改了两行 import，`check()` 读配置的那几个属性由 `regulatory._view` 给；它多读一个属性本库要到运行时才知道（`test_law3_the_vendored_check_reproduces_fmcw_radars_own_report` 会红）。
- **报告里给“按 SNR 的最远探测距离”**。A 行的 29.56 m 是本库由报告的 43.34 dB 按 R⁴ 与氧气项外推的（derived），不含 IF 高通与相噪底随距离的变化。
- A 行的泄漏：默认配置的报告里 `leakage` 是 null（模型关着），所以表的 §5 对 A 只写“没有建模”。要比两边的隔离需求，得请它在默认配置上开泄漏模型或给一份开着的报告（它的示例 08）。

## 里程碑 1：阵列与波束扫描下的感知

单天线够不到通信限值：满功率 22 dBm EIRP，离 40 dBm 差 18 dB（§4）。N 个单元相干合成 EIRP +20·log10(N)：8 个单元 +18 dB 正好顶到上限，届时通信读法比雷达读法多 20 dB 而不是 2 dB，B 行 SNR 按闭式是 30.75 + 18 = 48.8 dB（发射阵列，接收仍单天线）。要做的：
- vendor mmw60_trx 的阵列层（它 vendor 的 arraytx `geometry` / `beam`）与扇区扫描的时长（`assumed.ay_sector_sweep`、`edmg.short_ssw`、`dmg.ssw_frame`）；波束扫描把一帧分给多个方向，空口时间 × 扇区数。
- 那时再决定要不要 vendor `Frontend60`（它依赖 wifitrx 整包；fmcw_radar 0.2.8 的理由仍成立）。

## 这一轮发现、没有做的（按对结论的分量排）

1. **压缩判据只有闭式**（铁律 4 的缺口，表里写明）。单站“差 16 dB”这个结论压在它上面，而 10 dB 的余量是本库取的（`jcas60.compression_margin_db`；取 0 时差 6 dB）。做法：vendor fmcw_radar 的 `link/compression.py`（Rapp 限幅器与小信号增益，只依赖 numpy），在 `wavesim` 的接收端加一级，量泄漏从 −44 到 −24 dBm 时回波 SNR 掉多少，把余量从假设换成读数。
2. **20 dB 隔离下的结论要改口径**：ADC 抖不开时（热噪声 0.06 LSB），量化误差的 69 % 落在零多普勒的脊上，与 5 m 处的目标等高（−87.1 对 −86.1 dBm，§5）。闭式的“底抬高 14.18 dB”不成立（仿真 9.24 dB）。静止或慢速目标（|v| ≤ 4 个多普勒格 ≈ 1.5 m/s）的可见性没有量；加抖动（dither）能不能把脊打散也没有试。
3. **泄漏只有一条 0.1 ns 的路径**（`jcas60.leakage_delay_s`）。相噪判据对它敏感：0.2 ns 需要的隔离多 5.5 dB。天线罩或近处物体的反射是更长时延的“泄漏”，相关抵消更差；没有建。
4. **A 行的本振不是同一个**（铁律 2 的例外，表里标明）。fmcw_radar 有一份用 mmw60_trx ×10 本振剖面的配置（`radar_60g_x10lo.yaml`）；要让三行同一个本振，A 行换成那份配置的报告，ay 行的相噪从它的剖面读。这一轮按用户定的默认配置没有动。
5. **双站的同步没有仿真**：每个包的相位取自直达径（已仿真），定时与载波频偏假设已由前导对齐。±40 ppm 是 2.4 MHz，在 0.58 µs 的相关场里转 1.4 圈，不补偿时互补性全毁；mmw60_trx 的 STF 估计残差 < 5 kHz（它的铁律 3），量级上够，但没有接进来量。直达径上的相噪：闭式 1.26 dB 对仿真 0.89 dB（差 0.37 dB，闭式的“减掉场内平均相位”是近似）。
6. **感知场只用了 L-CEF**：L-STF 的组成公开材料里没有（spec_sources 的 ROADMAP），用上是 +4.95 dB；TRN 字段没有建（DMG 的单元长度有，结构没有）。追平 A 行的 SNR 要 ×18.2 的相关码片（§2），最自然的来源就是 TRN。
7. **目标是一个点**：没有杂波、多目标、多径；RCS 是 77 GHz 射线追踪范围的下端（−10 dBsm），上端 +5 dBsm 时三行一起 +15 dB（§8）。循环前缀之外（> 10.9 m）的目标落进来会怎样（窗里不再是循环移位）没有量。
8. **距离窗**：相关接收没有加窗，旁瓣 −14.6 dB；加频域窗换旁瓣的取舍（与 FMCW 的 Hann 同一种）没有做，所以“处理后的宽度 7.43 对 5.34 cm”两边不是同一种取舍。
9. **信道绑定**（2 / 4 信道：3.52 / 7.04 GHz，c/2B 4.26 / 2.13 cm）：带宽上赢过 4 GHz FMCW 的唯一途径，mmw60_trx 第一轮也定了不做。

## 不在范围内（用户定）

图形界面、Android；OFDM 波形做感知、成像。
