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
| ③ 同步 | 位同步、帧同步、导航电文 | `src/gnssrx/nav_msg.py` | ⬜ |
| ④ 星历 | 星历解码与卫星位置计算 | `src/gnssrx/ephemeris.py` | ⬜ |
| ⑤ 定位 | 伪距提取 + 最小二乘单点定位 | `src/gnssrx/pvt.py` | ⬜ |
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

值得一提的两个坑，都写在 `scripts/02_real_data_test.py` 注释里：

- **欠采样混叠**：模拟中频 4.092 MHz 高于 fs/2，混叠后的数字中频是 1.364 MHz，且多普勒符号相反
- **1-bit 数据的字节内比特序**：该数据集是 LSB 优先，排错的话相关峰会掉到 1/4

## 数据集

| 数据 | 来源 | 用途 |
|---|---|---|
| 合成中频信号 | 内置 `sim.py`，含 50 bps 导航电文与码多普勒 | 算法验证（真值已知） |
| 真实 GPS L1 采集 | Nottingham 数据集（1-bit，fs 5.456 MHz） | 真实信号验证 |

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
│   └── 02_real_data_test.py   # 真实 GPS L1 信号验证
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
