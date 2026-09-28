"""中频（IF）采样数据的读写。

常见数据集格式对照：
    SoftGNSS / Borre 教材配套  —— 8-bit 实数，采样率 38.192 MHz，IF 9.548 MHz
    部分采集卡                 —— 8-bit I/Q 交替（I0,Q0,I1,Q1,...）
    RTL-SDR rtl_sdr 输出       —— 8-bit I/Q 交替（无符号偏移 127，需减 127）

本模块只处理"实数 IF"和"I/Q 交替"两种最常见的格式。
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


def read_real_int8(path: str | Path, n_samples: int = -1, offset: int = 0) -> np.ndarray:
    """读取 8-bit 实数中频采样（SoftGNSS 数据集就是这个格式）。"""
    data = np.fromfile(str(path), dtype=np.int8, count=n_samples,
                       offset=offset * np.dtype(np.int8).itemsize)
    return data.astype(np.float64)


def read_iq_int8(path: str | Path, n_samples: int = -1,
                 signed: bool = True) -> np.ndarray:
    """读取 8-bit I/Q 交替数据，返回复数数组。signed=False 用于 RTL-SDR 原始输出。"""
    raw = np.fromfile(str(path), dtype=np.uint8 if not signed else np.int8, count=n_samples)
    if not signed:
        raw = raw.astype(np.float64) - 127.5
    iq = raw.astype(np.float64).reshape(-1, 2)
    return iq[:, 0] + 1j * iq[:, 1]


def write_real_int8(path: str | Path, samples: np.ndarray, scale: float = 32.0) -> None:
    """把浮点信号量化成 8-bit 实数写入文件，模仿真实采集卡的输出。"""
    q = np.clip(np.round(samples * scale), -128, 127).astype(np.int8)
    q.tofile(str(path))


def read_1bit_i(path: str | Path, n_samples: int = -1, offset_samples: int = 0,
                msb_first: bool = True) -> np.ndarray:
    """读取 1-bit 量化的实数中频采样（每个字节打包 8 个采样点，每点 1 bit）。

    1-bit 量化会带来约 2 dB 信噪比损失，但对捕获/跟踪影响很小 —— 早期采集卡
    与不少公开数据集（如 Nottingham GPS 数据集）都是这个格式。

    msb_first 决定一个字节内 8 个点的排列顺序，两种都要试：排错了序列会乱、
    相关峰出不来。整体极性反了则无所谓（只是 I/Q 反号，相关值取模后一样）。
    """
    byte_offset = offset_samples // 8
    if n_samples < 0:
        raw = np.fromfile(str(path), dtype=np.uint8, offset=byte_offset)
    else:
        raw = np.fromfile(str(path), dtype=np.uint8,
                          count=(n_samples + 7) // 8, offset=byte_offset)
    bits = np.unpackbits(raw, bitorder="big" if msb_first else "little")
    if n_samples >= 0:
        bits = bits[:n_samples]
    return bits.astype(np.float64) * 2.0 - 1.0


def probe(path: str | Path, fs: float) -> dict:
    """快速体检一个数据文件：时长、幅度、频谱峰值位置（用来确认 IF 设对没有）。"""
    data = read_real_int8(path)
    n = data.size
    spec = np.abs(np.fft.rfft(data[: min(n, int(fs * 0.01))] * np.hanning(
        min(n, int(fs * 0.01)))))
    freqs = np.fft.rfftfreq(min(n, int(fs * 0.01)), d=1.0 / fs)
    peak_hz = float(freqs[np.argmax(spec)])
    return {
        "samples": n,
        "duration_s": n / fs,
        "mean": float(data.mean()),
        "std": float(data.std()),
        "spectral_peak_hz": peak_hz,
    }
