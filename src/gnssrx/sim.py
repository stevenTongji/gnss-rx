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

from .ca_code import CODE_LENGTH, CODE_RATE_HZ, generate_ca

DEFAULT_FS = 4.092e6      # 采样率：4 samples/chip，够用且文件不大
DEFAULT_IF = 1.023e6      # 中频


def synthesize(
    prns: list[int],
    *,
    fs: float = DEFAULT_FS,
    f_if: float = DEFAULT_IF,
    n_ms: int = 100,
    cn0_dbhz: float = 45.0,
    amplitude: float = 1.0,
    seed: int = 20260928,
) -> tuple[np.ndarray, list[dict]]:
    """合成多颗卫星的 GPS L1 实数中频信号。

    返回 (信号, 真值列表)。真值里给出每颗星的真实码延迟与多普勒，用于验证捕获结果。
    """
    rng = np.random.default_rng(seed)
    spc = int(round(fs / 1e3))          # 每个 C/A 周期（1 ms）的采样点数
    n = spc * n_ms
    t = np.arange(n) / fs

    sig = np.zeros(n)
    truth: list[dict] = []

    for prn in prns:
        code = generate_ca(prn).astype(np.float64)
        delay_samples = int(rng.integers(0, spc))
        doppler = float(rng.uniform(-3500.0, 3500.0))
        phase = float(rng.uniform(0.0, 2 * np.pi))

        # 每个采样点落在哪个码片上（含码延迟，对 1023 取模实现周期延拓）
        chip_idx = (np.floor((np.arange(n) - delay_samples) * (CODE_RATE_HZ / fs))
                    .astype(np.int64) % CODE_LENGTH)
        chips = code[chip_idx]
        sig += amplitude * chips * np.cos(2 * np.pi * (f_if + doppler) * t + phase)

        truth.append({"prn": prn,
                      "code_delay_samples": delay_samples,
                      "doppler_hz": doppler})

    # 按 C/N0 反推噪声强度：N0 = Ps / 10^(C/N0/10)，噪声功率 = N0 * (fs/2)
    ps = amplitude ** 2 / 2.0
    n0 = ps / (10.0 ** (cn0_dbhz / 10.0))
    sigma = np.sqrt(n0 * fs / 2.0)
    sig += rng.normal(0.0, sigma, size=n)

    return sig, truth


def to_int8(sig: np.ndarray, scale: float = 32.0) -> np.ndarray:
    """量化成采集卡常见的 8-bit 实数格式。"""
    return np.clip(np.round(sig * scale), -128, 127).astype(np.int8)
