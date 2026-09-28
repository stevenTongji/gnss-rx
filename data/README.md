# 数据目录

> ⚠️ **本目录里的 `*.bin` 已被 `.gitignore` 排除，不会进 Git。**
> 中频采样数据动辄几百 MB，GitHub 单文件上限 100 MB，硬传一定失败。
> 换电脑/换人协作时，按下面的链接重新下载即可。

## data/raw/ 放什么

放**原始中频采样文件**（8-bit 实数或 I/Q 交替的二进制）。

| 文件 | 采样率 | 中频 | 时长 |
|---|---|---|---|
| （等待下载） | | | |

## 去哪找数据

**① 先跑仿真（零下载，推荐第一步）**
`gnssrx.sim.synthesize()` 直接合成，答案已知，适合验证算法。
仓库里的 `scripts/00_smoke_test.py` 就是这么做的。

**② 真实数据集**

- **SoftGNSS 配套数据**（Borre 教材，最经典）
  <https://github.com/TMBOC/SoftGNSS>
  README 里给了 `gnss0.bin`（60.5 s）和 `gnsa14.bin`（40 s）的网盘/Google Drive 链接，
  SiGe GN3S 采集卡，8-bit I/Q，38.192 MHz 采样，9.548 MHz 中频。
  同名 MATLAB 源码可以直接当"参考答案"对着读。

- **Darius Plaušinaitis 的样点**（1 s 短数据，专测捕获/跟踪）
  <http://gfix.dk/?p=188>
  16.3676 MHz 采样，4.1304 MHz 中频，int8。含 GPS L1 与 Galileo GIOVE 信号。

- **GNSS-SDR 官方示例**
  <https://gnss-sdr.org/> · <https://github.com/gnss-sdr/gnss-sdr>
  C++ 全栈实现，适合做工程参照（不建议照抄，会失去"自己写"的价值）。

- **上海交大 GNSS 数据集**（含 IMU，留给项目 B）
  <https://bat.sjtu.edu.cn/?p=1239>
  注意：原始 IF 数据 70 GB 级别，**别下整包**，只看它的测量值部分。

**③ 自己采（加分项，但要买硬件）**
RTL-SDR（约几十元）+ 有源 GPS 天线，采 GPS L1 / 北斗 B1I。
这一步写在简历里很值钱，但**放到跟踪调通之后再做**，别一上来就卡在硬件上。
