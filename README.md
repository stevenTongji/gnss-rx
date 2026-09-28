# gnss-rx · 从零写一个 GNSS 软件接收机

> 对应《导航×通信 岗位清单与学习地图》**Step 3 · 项目 A**
> 目标岗位：北斗基带 / 芯片算法、通信物理层、航天院所、智驾融合定位

## 这条链路在做什么

```
天线 → 射频前端 → [中频采样数据 .bin]
                          │
       ┌──────────────────┴──────────────────┐
       ▼                                     ▼
   ① 捕获 acquisition.py          ② 跟踪 tracking.py
   找出"天上有哪几颗星、           DLL + PLL 环路，逐毫秒
   码相位多少、多普勒多少"          咬住码和载波，输出 I/Q
       │                                     │
       └──────────────────┬──────────────────┘
                          ▼
              ③ 位同步 / 帧同步  nav_msg.py
                 解出导航电文里的星历
                          ▼
              ④ 伪距 + 最小二乘定位  pvt.py
                 输出经纬度与误差曲线
```

## 当前进度

| 模块 | 文件 | 状态 |
|---|---|---|
| C/A 码生成 | `src/gnssrx/ca_code.py` | ✅ 完成（含 GPS 三值特性自检） |
| 中频合成 / 读写 | `src/gnssrx/sim.py`、`io_if.py` | ✅ 完成 |
| 并行码相位 FFT 捕获 | `src/gnssrx/acquisition.py` | ✅ 完成（32 星搜索 0.15 s） |
| 跟踪环路 DLL + PLL | `src/gnssrx/tracking.py` | ⬜ 待做 |
| 位同步 / 帧同步 | `src/gnssrx/nav_msg.py` | ⬜ 待做 |
| 星历解码 | `src/gnssrx/ephemeris.py` | ⬜ 待做 |
| 伪距 + 最小二乘定位 | `src/gnssrx/pvt.py` | ⬜ 待做 |
| C++ 工程化移植 | `cpp/` | ⬜ 待做（第 6 周起） |

## 快速开始

```bash
cd ~/Developer/gnss-rx
uv sync                             # 安装依赖
uv run python scripts/00_smoke_test.py   # 冒烟测试：环境 + 算法自检
```

## 关键参数（改数据时要同步改）

| 参数 | 当前值 | 说明 |
|---|---|---|
| `fs` | 4.092 MHz | 采样率 = 4 采样点/码片 |
| `f_if` | 1.023 MHz | 中频 |
| `CODE_RATE_HZ` | 1.023 MHz | C/A 码速率，周期 1 ms |

换成真实数据集（如 SoftGNSS 的 38.192 MHz / 9.548 MHz）时，
**采样率和中频必须写成配置常量，不要在代码里硬编码**——这是后面所有模块的公共约定。
