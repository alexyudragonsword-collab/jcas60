# ROADMAP

**这份文件只写还没做的事。** 做完即删，并在 `CHANGELOG.md` 记一笔、`cairn/LOG.md` 记数字。里程碑 0（对照表）的结论在 `cairn/comparison-table.md`。

复现所有数字：`python examples/00_comparison_table.py`（约 30 s；§n 指它输出的第 n 节）。

## 报回上游（本库不改别的库；要它们做的在那个库开 issue，对方合并后本库重新 pin、验收、关 issue）

### spec_sources —— pin 在 `f9ce3f2`；[spec_sources#54](https://github.com/alexyudragonsword-collab/spec_sources/issues/54) 答了一半，issue 开着
- **802.11bf 在 60 GHz 的感知条目**：它读完了公开文件，等用户对照标准核对后另开 PR。到之前本库的感知场（L-CEF）、感知包（最短的 SC 包）、一帧 64 个包都是 assumed，措辞守铁律 7。到了之后：重新 pin，按条目重建感知帧，B、C 两行重出，再议措辞。
- **那次重新 pin 一并带进来的**（它建议不为这一点单独 pin）：`itu_p676.gamma_o_dry_sea_level` 在 `14e84fd` 已有 60.48 GHz 的点，14.993 dB/km。接入时删掉 `jcas60.oxygen_frequency_ghz`、改读 `@f_ghz=60.48`；5 m 处往返 10 m 动 0.0034 dB，是一次 numbers-move。
- **法规已答复（2026-11-25），没有定论**：五问里四问“原文没有”，FCC 23-35 没有讨论通信设备用自己的信号做感知；**KDB 没有查到**（检索页对它返回 503）。所以 ay 每行两个版本、雷达读法的空口占比 1.16 % … 50 % 的区间都留着。以后谁读到 KDB 或 FCC 对 WLAN sensing 的说法，在 #54 里补。
- 本库已登记在它的 `consumers.yaml`（at `8cb40e0`、pin `f9ce3f2`）；重新 pin 后在 #54 回报新的提交。

### mmw60_trx —— pin 在 `442fabd`；mmw60_trx#89 已验收
- 以后做 OFDM 感知时要 `ay/training.py` 本身：它记在自己的 ROADMAP（里程碑 3“按需”），到时在 #89 或新 issue 里说，照 `preamble.py` 的分法拆。
- L-STF 与 TRN 的组成它也没有公开出处（只有码片数）；有了、经 spec_sources 进来之后它再加，本库再 vendor。

### fmcw_radar —— pin 在 `c303566`；fmcw_radar#42 已验收
- A 行的泄漏：默认配置的报告里 `leakage` 是 null（模型关着），所以表的 §5 对 A 只写“没有建模”。要比两边的隔离需求，得请它在默认配置上开泄漏模型或给一份开着的报告（它的示例 08）。还没有开 issue。
- B、C 的 CFAR 损失是借 A 行的：三行的 Pd 0.9 距离都用报告的门限 12.24 dB（比噪声电平已知的 11.75 dB 高 0.49 dB，那是 A 行那张加 Hann 窗、4 倍补零的图上的数）。ay 的图不加距离窗、不补零，自己的损失应当更小；要算得 vendor 它的 `processing/cfarstats.py`，而那个文件顶上 import `RadarConfig`（`cfar_statistics_2d` 本身不用它）。换成 11.75 dB 三行都远 2.7–2.8 %。

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
