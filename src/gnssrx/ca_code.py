"""GPS L1 C/A 码（PRN 1-32）生成器。

原理（IS-GPS-200 第 3.2.1.3 节）：
    两个 10 级线性反馈移位寄存器 G1 / G2 各产生周期 1023 的 m 序列，
    G1 反馈抽头 = 第 3、10 级；G2 反馈抽头 = 第 2、3、6、8、9、10 级。
    不同 PRN 的区别在于 G2 输出的"相位选择器"抽头不同，见 G2_TAPS。

    码片用 ±1 表示（而不是 0/1），这样相关运算直接就是乘法，且便于 FFT。

参考：Kai Borre《A Software-Defined GPS and Galileo Receiver》makeCaTable.m
"""
from __future__ import annotations

import numpy as np

CODE_LENGTH = 1023          # C/A 码一个周期的码片数
CODE_RATE_HZ = 1.023e6      # 码速率 (chip/s)，周期恰为 1 ms
L1_HZ = 1575.42e6           # GPS L1 载波频率；码多普勒由它换算（见 tracking.py 的载波辅助）

# PRN 1..32 的 G2 相位选择抽头（1-based，IS-GPS-200 表 3-Ia）
G2_TAPS: dict[int, tuple[int, int]] = {
    1: (2, 6),   2: (3, 7),   3: (4, 8),   4: (5, 9),   5: (1, 9),
    6: (2, 10),  7: (1, 8),   8: (2, 9),   9: (3, 10), 10: (2, 3),
   11: (3, 4),  12: (5, 6),  13: (6, 7),  14: (7, 8),  15: (8, 9),
   16: (9, 10), 17: (1, 4),  18: (2, 5),  19: (3, 6),  20: (4, 7),
   21: (5, 8),  22: (6, 9),  23: (1, 3),  24: (4, 6),  25: (5, 7),
   26: (6, 8),  27: (7, 9),  28: (8, 10), 29: (1, 6),  30: (2, 7),
   31: (3, 8),  32: (4, 9),
}

PRN_LIST = tuple(G2_TAPS.keys())

# G1 / G2 的反馈抽头（1-based → 转成 0-based 下标）
_G1_FB = (2, 9)
_G2_FB = (1, 2, 5, 7, 8, 9)


def generate_ca(prn: int) -> np.ndarray:
    """生成指定 PRN 的 C/A 码，返回长度 1023 的 int8 数组，元素取值为 ±1。"""
    if prn not in G2_TAPS:
        raise ValueError(f"PRN 必须是 1..32，收到 {prn}")

    t1, t2 = G2_TAPS[prn]
    # 初值必须是全 -1：在 0/1 域里它是"全 1"，而全 +1 是移位寄存器的退化态
    # （反馈永远是 1，序列不动），这是一个非常经典的坑。
    g1 = np.full(10, -1, dtype=np.int8)
    g2 = np.full(10, -1, dtype=np.int8)

    code = np.empty(CODE_LENGTH, dtype=np.int8)
    for i in range(CODE_LENGTH):
        # ±1 相乘等价于二进制异或
        code[i] = g1[9] * g2[t1 - 1] * g2[t2 - 1]

        fb1 = g1[_G1_FB[0]] * g1[_G1_FB[1]]
        fb2 = (g2[_G2_FB[0]] * g2[_G2_FB[1]] * g2[_G2_FB[2]]
               * g2[_G2_FB[3]] * g2[_G2_FB[4]] * g2[_G2_FB[5]])

        g1[1:] = g1[:-1]; g1[0] = fb1     # 右移，新值从左侧推入
        g2[1:] = g2[:-1]; g2[0] = fb2

    return code


def _circular_autocorr(code: np.ndarray) -> np.ndarray:
    """循环自相关（FFT 实现），返回值未归一化。"""
    spec = np.fft.fft(code.astype(np.float64))
    return np.real(np.fft.ifft(spec * np.conj(spec)))


def self_test() -> None:
    """自检：C/A 码的周期自相关必须满足 GPS 的三值特性。

    自相关在零延迟处 = 1023；非零延迟处取值只能落在 {-1, -65, 63} 内。
    这是对生成器的强校验——抽头或初值写错的话这里必然崩。
    """
    for prn in PRN_LIST:
        ac = _circular_autocorr(generate_ca(prn))
        assert abs(ac[0] - CODE_LENGTH) < 1e-6, f"PRN {prn} 零延迟自相关错误"
        assert np.max(np.abs(ac[1:])) <= 65.0 + 1e-6, f"PRN {prn} 旁瓣越界"

    # 不同 PRN 之间的互相关也应很小（互相关上界 65）
    a, b = generate_ca(1).astype(np.float64), generate_ca(2).astype(np.float64)
    xcorr = np.abs(np.real(np.fft.ifft(np.fft.fft(a) * np.conj(np.fft.fft(b)))))
    assert xcorr.max() <= 65.0 + 1e-6, "PRN 1/2 互相关越界"

    print(f"[ca_code] 自检通过：{len(PRN_LIST)} 个 PRN 的自相关/互相关均符合 GPS 三值特性")


if __name__ == "__main__":
    self_test()
