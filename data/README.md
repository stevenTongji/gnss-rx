# 数据目录

> ⚠️ **本目录里的 `*.bin` 已被 `.gitignore` 排除，不会进 Git。**
> 中频采样数据动辄几百 MB（GitHub 单文件上限 100 MB），硬传一定失败。
> 换机器时按下面的链接重新下载即可。

## data/raw/ 放什么

放**原始中频采样文件**。

| 文件 | 采样率 | 中频 | 格式 | 时长 |
|---|---|---|---|---|
| `nottingham_gps_l1_1bit.bin` | 5.456 MHz | 1.364 MHz（混叠后） | 1-bit 实数 | ~80 s |
| `sige_gps_l1_8bit.dat` | 16.368 MHz | 4.092 MHz | 8-bit 有符号实数 IF | 原文件 ~120 s（本仓库只取前 40 s） |
| （合成数据不需要下载） | 4.092 MHz | 1.023 MHz | float64 | 由 `sim.py` 生成 |

## 数据集

### ① Nottingham GPS L1（真实信号，推荐从这里开始）

`scripts/02_real_data_test.py` 会自动使用它。

```bash
cd data/raw
curl -L -o nottingham_gps_l1_1bit.bin \
  "https://api.github.com/repos/JiaoXianjun/GNSS-GPS-SDR/contents/gps.samples.1bit.I.fs5456.if4092.bin" \
  -H "Accept: application/vnd.github.v3.raw"
```

- 约 55.8 MB，1-bit 量化（每字节 8 个采样点），**字节内 LSB 优先**
- 原始来源：jks.com/gps/gps.html（Nottingham GPS 数据集），
  本链接转自 [JiaoXianjun/GNSS-GPS-SDR](https://github.com/JiaoXianjun/GNSS-GPS-SDR)
- 1-bit 量化损失约 2 dB 信噪比，但对捕获/跟踪影响很小
- **欠采样注意**：模拟中频 4.092 MHz 高于 fs/2，混叠后数字中频为 1.364 MHz，
  且多普勒符号与真实值相反（详见 `scripts/02_real_data_test.py` 的注释）

### ①b SiGe GN3S（真实信号，卫星多、信号强，PVT 推荐）

来自 ION GNSS SDR 元数据标准的 **SiGe GN3S v3** 前端样本。这是经典的实数 IF 数据
（与本项目采集/跟踪链路原生兼容，无需改代码），可见 **PRN 1,7,8,9,11,17,28 共 7 颗**，
足以支撑 ≥4 星的 3D 定位。

- `fs = 16.368 MHz`，`IF = 4.092 MHz`（= fs/4），**8-bit 有符号实数 IF**，大端无关（8-bit）
- 编码 `TCA`（Two's Complement，按 `np.int8` 直接读取即可）
- 原文件 1.96 GB / ~120 s；本仓库只需前 40 s（约 654 MB），用 HTTP Range 分段取：

```bash
cd data/raw
# 只取前 40 秒（字节 0 .. 654719999），避免下整包 1.96 GB
curl -L -r 0-654719999 -o sige_gps_l1_8bit.dat \
  "https://sdr.ion.org/SiGe/SiGe_Bands-L1.dat"
```

- 元数据（采样率/中频/编码/采集位置）：<https://sdr.ion.org/SiGe/SiGe_Bands-L1.sdrx>
- **采集地：德国慕尼黑**（元数据 `lat=48.1715°, lon=11.8087°, h≈577 m`，2013-05-23 采集）——
  可用作 PVT 的真值参考（`scripts/07_sige_pvt_test.py` 即用慕尼黑先验）。
- 验证脚本：`scripts/05_sige_ephemeris_test.py`（星历）、PVT 脚本复用同一份数据

### ② 其他可选数据集

- **SoftGNSS 配套数据**（Borre 教材，60.5 s，8-bit I/Q，38.192 MHz）：
  <https://github.com/TMBOC/SoftGNSS>，README 里有网盘链接；
  同名 MATLAB 源码可直接当"参考答案"对着读。
- **Darius Plaušinaitis 样点**（1 s 短数据，测捕获/跟踪）：<http://gfix.dk/?p=188>
- **上海交大 GNSS 数据集**（含 IMU，留给项目 B）：<https://bat.sjtu.edu.cn/?p=1239>
  原始 IF 数据 70 GB 级别，**别下整包**。

### ③ 自己采（加分项，需要买硬件）

RTL-SDR（约几十元）+ 有源 GPS 天线，采 GPS L1 / 北斗 B1I。
建议放到跟踪调通之后再做，别一上来就卡在硬件上。
