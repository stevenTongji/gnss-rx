"""绘图公共配置：让 matplotlib 在 macOS 上正常显示中文。"""
from __future__ import annotations

import matplotlib as mpl

mpl.rcParams["font.sans-serif"] = [
    "PingFang SC", "Hiragino Sans GB", "Heiti SC", "Arial Unicode MS", "sans-serif",
]
mpl.rcParams["axes.unicode_minus"] = False   # 负号显示
