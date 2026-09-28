"""合成 GPS L1 中频信号——用来在没有真实数据之前先把捕获/跟踪算法跑通。

为什么先做仿真：
    仿真里"答案"是已知的（哪几颗星、码相位多少、多普勒多少），
    所以算法写错会立刻暴露；拿到真实数据时出问题，你分不清是数据的问题还是代码的问题。

用法：
    from gnssrx.sim import synthesize
    sig, truth = synthesize(prns=[7, 11, 19, 24, 30], n_ms=100, cn0_dbhz=45)
"""
from __future__ import annotations

import numpy as np

from .ca_code import CODE_LENGTH, CODE_RATE_HZ, L1_HZ, generate_ca

DEFAULT_FS = 4.092e6      # 采样率：4 samples/chip，够用且文件不大
DEFAULT_IF = 1.023e6      # 中频
NAV_BIT_RATE = 50.0       # 导航电文速率 50 bps → 每个比特占 20 ms
__all__ = ["DEFAULT_FS", "DEFAULT_IF", "L1_HZ", "NAV_BIT_RATE",
           "synthesize", "to_int8"]


def synthesize(
    prns: list[int],
    *,
    fs: float = DEFAULT_FS,
    f_if: float = DEFAULT_IF,
    n_ms: int = 100,
    cn0_dbhz: float = 45.0,
    amplitude: float = 1.0,
    nav_bits: bool = True,
    code_doppler: bool = True,
    seed: int = 20260928,
) -> tuple[np.ndarray, list[dict]]:
    """合成多颗卫星的 GPS L1 实数中频信号。

    返回 (信号, 真值列表)。真值里给出每颗星的真实码延迟、多普勒与码速率，用于验证捕获/跟踪。

    两个可选的物理效应（默认都开，因为它们正是跟踪环存在的理由）：

    nav_bits
        50 bps 导航电文，每 20 ms 可能把载波相位翻转 180°。
        所以载波环必须用对 180° 翻转不敏感的鉴相器（Costas），否则一遇比特翻转就失锁。

    code_doppler
        码速率被载波多普勒按同一比例牵引：f_code = 1.023 MHz × (1 + f_d / f_L1)。
        典型 3000 Hz 多普勒下码速率偏差约 1.95 chip/s，1 秒就漂出 2 个码片。
        码环若不跟踪，相关峰会一路滑走 —— 这是检验 DLL 是否真在工作的关键。
    """
    rng = np.random.default_rng(seed)
    spc = int(round(fs / 1e3))          # 每个 C/A 周期（1 ms）的采样点数
    n = spc * n_ms
    t = np.arange(n) / fs
    samples_per_chip = fs / CODE_RATE_HZ

    sig = np.zeros(n)
    truth: list[dict] = []

    for prn in prns:
        code = generate_ca(prn).astype(np.float64)
        # 码延迟取小数（而不是整数）：若恰好是整数个采样点，标称码索引会落在整数上，
        # 此时码多普勒哪怕只有 1e-6 的偏移，floor() 也会立刻把整段码错位 1 个采样点，
        # 造成系统性偏置。真实信号的码相位本来就是任意小数。
        delay_samples = float(rng.uniform(0.0, spc))
        doppler = float(rng.uniform(-3500.0, 3500.0))
        phase = float(rng.uniform(0.0, 2 * np.pi))

        # 码相位（单位：码片）随采样点线性推进；码多普勒体现为推进速率的微小偏移
        code_rate_hz = CODE_RATE_HZ * (1.0 + doppler / L1_HZ) if code_doppler else CODE_RATE_HZ
        chip_idx = (np.floor(np.arange(n) * (code_rate_hz / fs)
                             - delay_samples / samples_per_chip)
                    .astype(np.int64) % CODE_LENGTH)
        chips = code[chip_idx]

        # 导航电文：每颗星独立的随机比特流，且比特边界各自随机偏移
        if nav_bits:
            bit_len = int(round(fs / NAV_BIT_RATE))
            n_bits = n // bit_len + 2
            bits = rng.integers(0, 2, size=n_bits).astype(np.float64) * 2.0 - 1.0
            offset = int(rng.integers(0, bit_len))
            nav = bits[(np.arange(n) + offset) // bit_len]
        else:
            nav = np.ones(n)

        sig += amplitude * nav * chips * np.cos(2 * np.pi * (f_if + doppler) * t + phase)

        truth.append({"prn": prn,
                      "code_delay_samples": delay_samples,
                      "doppler_hz": doppler,
                      "code_rate_hz": code_rate_hz,
                      "cn0_dbhz": cn0_dbhz})

    # 按 C/N0 反推噪声强度：N0 = Ps / 10^(C/N0/10)，噪声功率 = N0 * (fs/2)
    ps = amplitude ** 2 / 2.0
    n0 = ps / (10.0 ** (cn0_dbhz / 10.0))
    sigma = np.sqrt(n0 * fs / 2.0)
    sig += rng.normal(0.0, sigma, size=n)

    return sig, truth


def to_int8(sig: np.ndarray, scale: float = 32.0) -> np.ndarray:
    """量化成采集卡常见的 8-bit 实数格式。"""
    return np.clip(np.round(sig * scale), -128, 127).astype(np.int8)
