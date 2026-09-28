"""gnssrx —— 从零开始的 GNSS 软件接收机（项目 A）。

模块路线（对应《导航×通信_岗位清单与学习地图》Step 3 项目 A）：
    采集/读取 IF 数据 → 捕获 → 跟踪 → 位同步/帧同步 → 解星历 → 伪距 → 最小二乘定位
"""
from .ca_code import generate_ca, CODE_LENGTH, CODE_RATE_HZ, PRN_LIST

__all__ = ["generate_ca", "CODE_LENGTH", "CODE_RATE_HZ", "PRN_LIST"]
__version__ = "0.1.0"
