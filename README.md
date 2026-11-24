# jcas60

60 GHz 通感一体：在同一套 60 GHz 前端数字上，用通信波形（802.11ay 单载波的 Golay 训练场）做感知，和专用的 FMCW 雷达比——各自能做到什么、代价是什么。

本库的说法是“用 ay 波形做感知”；spec_sources 里还没有 802.11bf 的条目，这里不作与它一致性有关的声称。

```
pip install -e ".[test]"
python examples/00_comparison_table.py     # 里程碑 0 的对照表，约 30 s
pytest && ruff check src tests tools examples
```

里程碑 0 的结论（5 m 处的人体；前提与每个数的出处见 `cairn/comparison-table.md` 与示例输出）：

| | FMCW（fmcw_radar 默认） | ay 单站，通信限值 | ay 单站，雷达限值 |
|---|---:|---:|---:|
| 距离分辨率 c/2B | 3.75 cm | 8.52 cm | 8.52 cm |
| EIRP（各按各的上限） | 14 dBm | 22 dBm（单天线到不了 40 dBm） | 20 dBm |
| SNR | 43.3 dB | 30.8 dB | 28.8 dB |
| 占用的空口时间（30 帧/s） | 19.2 % | 1.16 % | 1.16 % … 50 % |
| 接收 ADC | 10 MS/s × 12 bit | 3.52 GS/s × 2 路 | 同左 |
| 单站需要的收发隔离 | —（泄漏在高通之外） | 56 dB，查到 40 / 20 dB | 54 dB |

- 它是 [mmw60_trx](https://github.com/alexyudragonsword-collab/mmw60_trx)（通信链路）与 [fmcw_radar](https://github.com/alexyudragonsword-collab/fmcw_radar)（FMCW 雷达）的姊妹库；外部的数只从 [spec_sources](https://github.com/alexyudragonsword-collab/spec_sources) 进来，三者按完整 SHA vendor 在 `src/jcas60/vendor/`。
- 规则在 `AGENTS.md`，没做的在 `ROADMAP.md`，每轮的数字在 `cairn/LOG.md`。
- 上游的检出放在 `/home/user/<库>` 或用 `JCAS60_SIBLING_<库>` 指路后：`python tools/vendor_check.py --strict`、`python tools/upstream_report.py check --strict`、`python tools/fingerprint.py check`。
