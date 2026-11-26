# jcas60

60 GHz 通感一体：在同一套 60 GHz 前端数字上，用通信波形（802.11ay 单载波的 Golay 训练场）做感知，和专用的 FMCW 雷达比——各自能做到什么、代价是什么。

本库的说法是“用 ay 波形做感知”，不作与 802.11bf 一致性有关的声称；感知帧逐项对过 spec_sources 里 802.11bf 的条目（取自工作组的公开文件），结果在示例输出的 §11。

```
pip install -e ".[test]"
python examples/00_comparison_table.py     # 对照表，约 70 s
python examples/01_receiver_and_window.py  # 表引出的十一个问题：限幅器、增益档、窗外的目标、快目标、粗 ADC、泄漏的时延、同一个本振、双站的频偏、距离窗、第二个目标、为感知设计的 TRN 波形，约 145 s
pytest && ruff check src tests tools examples
```

结论（5 m 处的人体；前提与每个数的出处见 `cairn/comparison-table.md` 与示例输出）：

| | FMCW（fmcw_radar 默认） | ay，只用 L-CEF | ay，用 TRN 场 |
|---|---:|---:|---:|
| 距离分辨率 c/2B | 3.75 cm | 8.52 cm | 8.52 cm |
| EIRP（各按各的上限；ay 按通信限值） | 14 dBm | 22 dBm（单天线到不了 40 dBm） | 22 dBm |
| SNR | 43.3 dB | 30.7 dB | 43.8 dB |
| 占用的空口时间（30 帧/s） | 19.2 % | 1.16 % | 4.51 % |
| 接收 ADC | 10 MS/s × 12 bit | 3.52 GS/s × 2 路 | 同左 |
| 单站需要的收发隔离（灵敏度压低 1 dB） | —（泄漏在高通之外） | 43.2 dB，查到 40 / 20 dB | 同左 |

按雷达限值读，ay 的 EIRP 是 20 dBm，SNR 各低 2 dB，空口时间的上界是 50 %；通信设备顺带感知适用哪一款，法规原文没有定论，对照表两种读法都出。

- 它是 [mmw60_trx](https://github.com/alexyudragonsword-collab/mmw60_trx)（通信链路）与 [fmcw_radar](https://github.com/alexyudragonsword-collab/fmcw_radar)（FMCW 雷达）的姊妹库；外部的数只从 [spec_sources](https://github.com/alexyudragonsword-collab/spec_sources) 进来，三者按完整 SHA vendor 在 `src/jcas60/vendor/`。
- 规则在 `AGENTS.md`，没做的在 `ROADMAP.md`，每轮的数字在 `cairn/LOG.md`。
- 上游的检出放在 `/home/user/<库>` 或用 `JCAS60_SIBLING_<库>` 指路后：`python tools/vendor_check.py --strict`、`python tools/upstream_report.py check --strict`、`python tools/fingerprint.py check`。
