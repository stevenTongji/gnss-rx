# gnss-rx

**中文** | [English](./README_EN.md)

从零实现的 GPS L1 C/A 软件接收机。不依赖任何现成 GNSS 库，从 C/A 码生成、并行码相位捕获、DLL/PLL 跟踪环路，一路手写到位同步、星历解码与最小二乘定位。

![跟踪环路验证结果](docs/figures/tracking.png)

## 为什么做这个项目

面向北斗 / GNSS 与无线通信算法岗位的能力准备。招聘 JD 里反复出现的能力项——捕获与跟踪、环路参数权衡、RTK/PPP、状态估计、C++ 工程化——在这里逐项落地为可运行、可验证的代码。

设计上遵循一条原则：**所有项目共用一条数据链和一个代码仓库**，前一级的输出是后一级的输入。这样两年里只做一件事，但简历上能写出三层深度，面试时能顺着一条线讲到底。

## 当前进展

| 阶段 | 模块 | 文件 | 状态 |
|---|---|---|---|
| 基础 | C/A 码生成（PRN 1–32） | `src/gnssrx/ca_code.py` | ✅ |
| 基础 | 中频信号仿真 / 读写 | `src/gnssrx/sim.py`、`io_if.py` | ✅ |
| ① 捕获 | 并行码相位 FFT + 平方环精频估计 | `src/gnssrx/acquisition.py` | ✅ |
| ② 跟踪 | 码环 DLL + 载波环 PLL + 载波辅助 | `src/gnssrx/tracking.py` | ✅ |
| ③ 同步 | 位同步、帧同步、导航电文 | `src/gnssrx/nav_msg.py` | ✅ |
| ④ 星历 | 星历解码 + 卫星位置（Sagnac / 相对论修正） | `src/gnssrx/ephemeris.py` | ✅ |
| ⑤ 定位 | 伪距提取 + 最小二乘 SPS（≥4 星）+ 3 星地球约束 | `src/gnssrx/pvt.py` | ✅ |
| ⑥ 工程化 | 核心模块 C++ 移植（Eigen + CMake） | `cpp/` | ⬜ |
| 进阶 | GNSS/INS 紧耦合、5G NR 定位与 NTN | — | ⬜ |

## 验证结果

所有指标均来自**真值已知的仿真数据**（`scripts/01_tracking_test.py` 可复现），不是估计值。

| 指标 | 结果 | 说明 |
|---|---|---|
| 多普勒估计残差 | **0.0 Hz** | 粗捕获 ±125 Hz 经平方环精频后收敛 |
| C/N₀ 估计 | 44.6 dB-Hz vs 真值 45.0 | 绝对误差 **−0.37 dB** |
| 载波鉴相器抖动 | 8.2°–10.8° | 测量噪声理论底约 7.2° |
| 平滑后真实相位抖动 | **1.31°** | 理论值 1.6° |
| 码环锁定保持 | \|I_P\| 保持 103%–125% | 1 秒内码多普勒漂移约 2 码片 |
| 捕获灵敏度 | 约 37 dB-Hz | 45 dB-Hz 需 5 ms 积分，37 dB-Hz 需 20 ms |

对照实验（同一颗星，代码路径开关）：关闭载波辅助后 C/N₀ 从 44.2 掉到 25.4 dB-Hz，
证明码环确实在跟踪码多普勒，而非依赖载波环兜底。

捕获结果（32 个 PRN × 全部多普勒频点，5 ms 非相干累积）：

![捕获结果](docs/figures/acquisition.png)

## 真实信号验证

在公开的 **Nottingham GPS L1 采集数据集**（1-bit，fs 5.456 MHz）上运行同一套代码：

![真实信号验证](docs/figures/real_data.png)

| 指标 | 结果 |
|---|---|
| 检出卫星数 | 10 颗（32 个 PRN 全搜） |
| 稳定锁定 | **9 / 10** |
| 最强卫星 C/N₀ | 47.2 dB-Hz（PRN 30） |
| 星座图 | 清晰可见导航电文的 ±1 跳变 |

### 星历解码（真实数据）

`scripts/04_ephemeris_test.py` 在 40 秒真实数据上完成位同步 → 帧同步 → 剥奇偶 → 解码子帧 1/2/3。
剥奇偶采用 IS-GPS-200 的 D30\* 规则（第 k 字数据位由其「前一个字的 D30」反相），子帧号由
TOW 结构性推出（`SFID = (TOW/6 mod 5)`，TOW=0→5），不依赖 HOW 内的 SFID 位。

判据全部不靠外部真值：`IODE2 == IODE3 == (IODC & 0xFF)`、`sqrt(A)` ≈ 5153.6（半长轴 ≈ 26560 km）、
`e` ≈ 0.01、`i0` ≈ 55°、卫星 ECEF 半径 ≈ 26560 km。

| 数据集 | 干净卫星（自检通过） | 说明 |
|---|---|---|
| SiGe GN3S v3（8-bit） | **PRN 1, 7, 8, 9, 11, 17, 28（7 颗）** | 与 ION 元数据标注的可见卫星完全一致 |
| Nottingham（1-bit） | **PRN 30, 29, 31, 1, 21, 5, 13, 23, 16（9 颗）** | sqrtA≈5153.6、e≈0.002–0.017、i0≈+55° 全部合理 |

> **关键 bug 修复（`ephemeris.py::strip_parity`）**：IS-GPS-200 的奇偶规则为
> `d_i = r_i ⊕ D30*(前字)`，而 D30\* 是**含 Costas 反相 c 的接收比特**，正确推导得
> `d_i = r_i ⊕ r_prevD30`（**c 自动抵消**）。原实现写成 `data ^ D30* ^ c`，对首字正确、
> 对后续字**多异或了一次 c**：c=0 的通道恰好正确，c=1 的通道整段数据被反相 → 解出"垃圾轨道"。
> 1-bit 数据里 c=1 比例高，导致大量真实卫星被误判成"C/A 互相关鬼影"而被丢弃。
> 修正（D30\* 同样带上 c）后，SiGe 由 3 颗恢复到 **7 颗**、Nottingham 由 2 颗恢复到 **9 颗**。

**奇偶校验**（`ephemeris.py::parity_failures`）：实现了 GPS LNAV 的 6 位奇偶校验（DF(2) 系数
由真实数据 RANSAC 反推、与 IS-GPS-200 方程一致），可判断一个子帧是否含可检出的比特错误。
实测：SiGe 全部子帧 **0 失败**（干净）；Nottingham 干净通道全 0，而互相关假锁的垃圾通道 **9/9 失败**。
`decode_subframes` 会在同类子帧的多个副本中优先取「奇偶失败最少」的一份。

值得一提的两个坑，都写在 `scripts/02_real_data_test.py` 注释里：

- **欠采样混叠**：模拟中频 4.092 MHz 高于 fs/2，混叠后的数字中频是 1.364 MHz，且多普勒符号相反
- **1-bit 数据的字节内比特序**：该数据集是 LSB 优先，排错的话相关峰会掉到 1/4

### 导航电文同步（位同步 → 帧同步）

在 30 秒真实信号上，三颗星全部完成同步：

| 指标 | 结果 |
|---|---|
| 位同步 | 能量法与跳变统计法**两种独立方法结果一致**，定位到 1 ms 精度 |
| 跳变法置信度 | 11.5–20.0（均匀分布应为 1，实测 92% 的比特跳变落在同一余数上） |
| 帧同步 | 每颗星找到 **5 个子帧**，间隔严格等于 300 比特 = 6.000 秒 |
| 子帧起始时刻 | 4.6 / 10.6 / 16.6 / 22.6 / 28.6 s（三颗星完全一致） |

"间隔必须恰好 300 比特"是一条极硬的约束——8 比特随机匹配的概率是 1/128，
但要连续多次间隔 300 命中几乎不可能，所以它同时验证了捕获、跟踪、位同步整条链路。

![导航电文同步](docs/figures/nav_msg.png)

### 单点定位（PVT：真实数据 + 合成验证）

**真实数据（Nottingham 1-bit，fs 5.456 MHz）**：完整链路 捕获 → 48 s 跟踪 → 星历解码
→ 伪距 → 加权最小二乘定位 已跑通，**9 颗干净卫星**（PRN 30/29/31/1/21/5/13/23/16）做完整 SPS，
定位结果 **52.93°N / 1.17°W / 高程 ≈ 0**，与公开采集点偏差 **1.96 km**，各星伪距残差 **±2 m**。
（SiGe 8-bit 数据集解出 **7 颗**干净卫星 PRN 1/7/8/9/11/17/28，采集地德国慕尼黑
48.1715°N/11.8087°E；其 PVT 残差约 ±40 km（≈0.1 ms）——奇偶校验确认其子帧 0 失败、星历无比特错误，
故该残余来自量测/数据侧，属该数据集固有水平；见 `scripts/07_sige_pvt_test.py`。Nottingham 则达米级。）

> **两个符号 bug（本轮定位精度从"数百 km"降到"米级"的关键）**，均由合成自检脚本
> `scripts/15_timing_selfcheck.py` 严格定位（该脚本用已知码相位/已知子帧时刻的合成信号，
> 逐项核对 `t_user` 与真实到达时刻之差）：
> 1. **码相位符号**：`t_user = m_blk/1000 − τ/1.023e6`（原为 `+`）。本跟踪器的码相位 τ 与
>    「信号码历元在块内的偏移」反号，用加号会引入随卫星变化的 ~1 ms 误差。
> 2. **卫星钟差符号**：`pseudorange = c·(t_user − t_sat) + c·δt_sat`（原为 `−`）。
>    物理推导 `c(t_user − t_sat) = ρ + c·b − c·δt_sat`，故还原到「已钟差改正的伪距」需取加号。
>
> 修正前 Nottingham 偏差 ~298 km、SiGe ~数百 km；修正后 Nottingham **1.96 km / 残差 ±2 m**。

**冷启动（无需任何先验位置）**：`pvt.solve_cold_start` 在全球经纬网格上搜索，代价函数为
「先用候选点消整数毫秒模糊（逐星取整）→ 残差散布」——由于 `ρ_raw` 的跨星差含整数毫秒的
精确台阶（±k·c ms），取整会把它**精确抵消**，故真位置处残差≈0、其余位置骤增，构成尖锐极小值；
粗搜（1°）→ 细搜 → 迭代最小二乘精化。Nottingham **不给任何先验**即得
52.93491°N / 1.16503°W（距采集点 1.96 km），与给先验的结果一致。

**合成验证（`scripts/13_synthetic_4sat_pvt.py`，已知真值）**：用已知真值场景走完整链路
（真实距离 → 注入整数毫秒模糊度 → `build_measurement` → 最小二乘）严格证明定位引擎正确：

| 场景 | 位置误差 | 说明 |
|---|---|---|
| 6 星无模糊，直接 WLS | **0.024 m** | 求解器数学正确 |
| 6 星 + 注入共同整毫秒歧义（K=434579） | **0.024 m** | 整数毫秒消模糊正确（GPS 周量级） |
| 4 星子集 | **0.024 m** | 完整 SPS 正确 |
| 码相位噪声 σ = 2 / 5 / 10 m | 6.3 / 15.8 / 31.5 m | 误差 ≈ σ·GDOP，符合理论 |
| 3 星 + 地球约束 | ~11 km | 两交点歧义，印证需 ≥4 星 |

**修复**：GPS 周号在电文里仅 10 比特（WN mod 1024），每 1024 周回卷；本数据解出 717，
真实扩展周号为 1741（2013），已在 `ephemeris.py` 按参考周号就近还原（该回卷不影响卫星
位置——位置只依赖周内 tk——但修正后日期一致）。

## 数据集

| 数据 | 来源 | 用途 |
|---|---|---|
| 合成中频信号 | 内置 `sim.py`，含 50 bps 导航电文与码多普勒 | 算法验证（真值已知） |
| 真实 GPS L1 采集 | SiGe GN3S v3（8-bit，fs 16.368 MHz，PRN 1/7/8/9/11/17/28） | 星历（7 干净星）/ PVT |
| 真实 GPS L1 采集 | Nottingham 数据集（1-bit，fs 5.456 MHz） | 捕获/跟踪/星历（9 干净星）/ PVT |

真实数据文件体积较大，不纳入版本管理；下载方式见 `data/README.md`。

## 仓库结构

```
gnss-rx/
├── src/gnssrx/
│   ├── ca_code.py       # C/A 码生成（含 GPS 三值特性自检）
│   ├── sim.py           # 中频信号仿真：导航电文、码多普勒、C/N₀ 可调
│   ├── io_if.py         # 中频数据读写（实数 / I-Q 交替 / 1-bit）
│   ├── acquisition.py   # 并行码相位 FFT 捕获 + 平方环精频估计
│   ├── tracking.py      # DLL + PLL 二阶环路 + 载波辅助
│   └── plotting.py      # 作图公共配置
├── scripts/
│   ├── 00_smoke_test.py       # 环境 + C/A 码 + 捕获验证
│   ├── 01_tracking_test.py    # 跟踪环路验证（真值对照）
│   ├── 02_real_data_test.py   # 真实 GPS L1 信号验证（Nottingham）
│   ├── 03_nav_msg_test.py     # 导航电文同步验证
│   ├── 04_ephemeris_test.py   # 星历解码验证（Nottingham）
│   ├── 05_sige_ephemeris_test.py  # SiGe 数据多通道星历解码
│   ├── 06_hackrf_ephemeris_test.py # HackRF 数据星历解码
│   ├── 07_sige_pvt_test.py    # SiGe 真实数据端到端 PVT
│   ├── 08_pvt_inspect.py      # PVT 诊断（复用缓存）
│   ├── 09_pvt_ambiguity_validation.py # 整毫秒歧义/3 星歧义验证
│   ├── 10_ephemeris_audit.py  # 全通道星历审计（找第 4 颗星）
│   ├── 11_3sat_geometry.py    # 3 球面交点几何定界
│   ├── 12_time_and_geometry.py # GPS 时基 / 卫星星下点核验
│   ├── 13_synthetic_4sat_pvt.py  # 合成 4 星 PVT 验证（已知真值）
│   ├── 14_nottingham_pvt.py   # Nottingham 端到端 PVT（含 --diag 逐通道诊断）
│   └── 15_timing_selfcheck.py  # 时基自检：合成信号核对 t_user 口径（定位两个符号 bug）
├── docs/figures/        # 结果图
└── data/                # 数据目录（*.bin 不入库）
```

## 快速开始

```bash
git clone git@github.com:stevenTongji/gnss-rx.git
cd gnss-rx
uv sync                                      # 安装依赖（需 Python ≥ 3.11）

uv run python scripts/00_smoke_test.py       # 捕获验证
uv run python scripts/01_tracking_test.py    # 跟踪验证（约 4 秒）
uv run python scripts/02_real_data_test.py   # 真实信号验证（需先下载数据，见 data/README.md）
```

## 算法说明

各模块的实现细节与公式推导写在源码注释里，这里只列关键点：

- **捕获**：一次 FFT 同时搜索全部 1023 个码相位，配合多普勒串行搜索；判据用峰/次峰比。
  45 dB-Hz 下 1 ms 相干积分余量不足，需 5 ms 非相干累积。
- **精频估计**：码剥离 → 逐毫秒复相关值 → 平方消除导航电文 → 补零 FFT。
  这一步不可省略——二阶 PLL 无法在合理时间内自行牵引上百赫兹频差。
- **码环**：归一化、非相干超前减滞后鉴相器，增益 −2 /码片，线性区 ±(1 − d/2) 码片。
- **载波环**：Costas 鉴相器 `atan(Q/I)`，对导航电文 180° 翻转不敏感，代价是留有 180° 相位模糊。
- **环路滤波器**：二阶 PI，系数取 Kaplan 公式 ω_n = 8ζB_n/(4ζ²+1)、τ₂ = 2ζ/ω_n、τ₁ = 1/ω_n²。
  鉴相器输出已归一化到被控量本身，故环路总增益为 1。
- **载波辅助**：码速率由载波多普勒推算 f_code = 1.023 MHz × (1 + f_d / f_L1)，
  码环只处理残差。3000 Hz 多普勒对应码速率偏 1.95 码片/秒，无辅助时 2 Hz 码环会失锁。

## 参考文献

- K. Borre, D. Akos, et al., *A Software-Defined GPS and Galileo Receiver: A Single-Frequency Approach*, Birkhäuser, 2007.
- E. D. Kaplan, C. Hegarty, *Understanding GPS/GNSS: Principles and Applications*, 3rd ed., Artech House, 2017.
- IS-GPS-200, *NAVSTAR GPS Space Segment / User Segment Interfaces*.
- T. D. Barfoot, *State Estimation for Robotics*, Cambridge, 2017.
- 3GPP TS 38.305（NR 定位）、TS 38.821（NTN 非地面网络）。

## 许可

MIT
