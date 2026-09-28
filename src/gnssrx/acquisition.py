"""并行码相位 FFT 捕获（Borre《软件定义的 GPS 与伽利略接收机》第 6 章算法）。

核心思想：
    穷举多普勒频点 × 一次 FFT 搞定全部 1023 个码相位。
    1 ms 数据做一次 FFT 相关，就能同时试出所有码相位——这是软件接收机之所以可行的关键。

    对某个多普勒候选 f：
        1) 把中频信号下变频到基带：x[n] * exp(-j*2π*f*n/fs)
        2) acc = |IFFT( FFT(基带信号) * conj(FFT(本地码)) )|²
        3) acc 的最大值位置 = 该卫星的码相位，最大值对应的 f = 多普勒

判据：
    峰值 / 次峰值 > 门限（通常 2.0~2.5）才算捕获成功。
    只看绝对峰值会被强信号漏进来的旁瓣骗到。
"""
from __future__ import annotations

import numpy as np

from .ca_code import CODE_LENGTH, CODE_RATE_HZ, PRN_LIST, generate_ca

TWO_PI = 2.0 * np.pi


def code_replica(prn: int, fs: float, tau_samples: float = 0.0) -> np.ndarray:
    """把 1023 个码片按采样率重采样成 1 ms 长的本地码副本。

    tau_samples 是本地码相位偏移（采样点，可为小数）。捕获给出的码延迟 D 对应
    tau = (-D) mod spc —— 见 tracking.py 的符号约定。
    """
    spc = int(round(fs / 1e3))
    idx = (np.floor((np.arange(spc) + tau_samples) * (CODE_RATE_HZ / fs))
           .astype(np.int64)) % CODE_LENGTH
    return generate_ca(prn).astype(np.float64)[idx]


def acquire_prn(
    data: np.ndarray,
    fs: float,
    f_if: float,
    prn: int,
    *,
    ms: int = 1,
    doppler_half_range: float = 6000.0,
    doppler_step: float = 500.0,
    replica: np.ndarray | None = None,
) -> dict:
    """对单颗卫星做捕获，返回最佳（多普勒, 码相位）及其置信度指标。"""
    spc = int(round(fs / 1e3))
    x = np.asarray(data, dtype=np.float64)
    if x.size < spc * ms:
        raise ValueError(f"数据不足：需要 {spc * ms} 点，实际 {x.size} 点")

    if replica is None:
        replica = code_replica(prn, fs)
    code_fft = np.conj(np.fft.fft(replica))

    ts = 1.0 / fs
    phase = 2.0 * np.pi * ts * np.arange(spc)

    # 相关峰主瓣宽度约 ±1 码片，判次峰时要挖掉它旁边的点
    guard = max(1, int(round(2.0 * fs / CODE_RATE_HZ)))
    all_idx = np.arange(spc)

    acc = np.empty(spc)
    best: dict | None = None

    for f_d in np.arange(-doppler_half_range, doppler_half_range + doppler_step, doppler_step):
        acc[:] = 0.0
        for m in range(ms):                       # 非相干累加，抗噪声
            block = x[m * spc:(m + 1) * spc]
            iq = block * np.exp(-1j * (f_if + f_d) * phase)
            acc += np.abs(np.fft.ifft(np.fft.fft(iq) * code_fft)) ** 2

        k = int(np.argmax(acc))
        peak = float(acc[k])
        circular = (all_idx - k + spc // 2) % spc - spc // 2
        second = float(acc[np.abs(circular) > guard].max())

        if best is None or peak > best["peak"]:
            best = {
                "prn": prn,
                "doppler_hz": float(f_d),
                "code_phase_samples": k,
                "code_phase_chips": float(k * CODE_RATE_HZ / fs),
                "peak": peak,
                "second_peak": second,
                "peak_ratio": float(peak / second) if second > 0 else np.inf,
                "peak_over_mean": float(peak / acc.mean()) if acc.mean() > 0 else np.inf,
            }

    assert best is not None
    return best


def acquire_all(
    data: np.ndarray,
    fs: float,
    f_if: float,
    prns=PRN_LIST,
    **kwargs,
) -> list[dict]:
    """对所有 PRN 依次捕获，按峰值/次峰值比从高到低返回。"""
    results = [acquire_prn(data, fs, f_if, p, **kwargs) for p in prns]
    return sorted(results, key=lambda r: r["peak_ratio"], reverse=True)


def detect(results: list[dict], threshold: float = 2.5) -> list[dict]:
    """按门限筛出"捕获成功"的卫星。"""
    return [r for r in results if r["peak_ratio"] >= threshold]


def refine_doppler(
    data: np.ndarray,
    fs: float,
    f_if: float,
    prn: int,
    *,
    code_phase_samples: int,
    coarse_doppler_hz: float,
    ms: int = 50,
    fft_size: int = 16384,
) -> float:
    """捕获后的精频估计 —— 跟踪环能不能牵进去，全看这一步。

    为什么必须做：
        捕获的多普勒分辨率就是搜索步长（典型 250-500 Hz），而二阶 PLL 靠 atan
        鉴相器自行牵引上百赫兹的频差要花好几秒（牵入过程会反复滑周）。
        真实接收机都是"粗捕获 + 精频估计"两步走，SoftGNSS 里同样有这一步。

    做法（BPSK 载波恢复里经典的"平方环"思路）：
        1) 用捕获到的码相位把码剥离，剩下的就是一个被 ±1 电文调制的中频单音
        2) 用粗多普勒下变频到基带，逐毫秒取一个复相关值
           z[m] ∝ d[m] · exp(j·2π·Δf·m·T)，Δf 是残余多普勒
        3) 取平方 z[m]² → (±1)² = 1，电文调制被消掉，频率翻倍成 2Δf
        4) 补零 FFT 找峰值 → Δf，分辨率可达零点几赫兹

    ms 越大分辨率越高，但码多普勒会让相关幅度缓慢下降；
    50 ms 内码只漂约 0.1 码片（幅度损失约 10%），足够。
    """
    spc = int(round(fs / 1e3))
    if data.size < spc * ms:
        raise ValueError(f"数据不足：需要 {spc * ms} 点")

    tau = float((-code_phase_samples) % spc)
    replica = code_replica(prn, fs, tau)      # 必须按捕获到的码相位对齐，否则剥离不掉码
    k = np.arange(spc)
    ts = 1.0 / fs
    f_lo = f_if + coarse_doppler_hz

    z = np.empty(ms, dtype=np.complex128)
    for m in range(ms):
        seg = data[m * spc:(m + 1) * spc]
        iq = seg * np.exp(-1j * TWO_PI * f_lo * (m * spc + k) * ts)
        z[m] = iq @ replica

    # 平方去掉电文调制，频率翻倍
    spec = np.abs(np.fft.fft(z ** 2, n=fft_size))
    freqs = np.fft.fftfreq(fft_size, d=1e-3)      # z 的采样率是 1 kHz（每毫秒一个）
    peak_hz = float(freqs[int(np.argmax(spec))])

    return coarse_doppler_hz + peak_hz / 2.0
