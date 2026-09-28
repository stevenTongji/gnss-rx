"""导航电文同步：位同步（bit sync）→ 帧同步（frame sync）。

================================================================================
一、这一步的输入和输出
================================================================================
输入：跟踪环每毫秒吐出的 I_P（一个实数序列，长度 = 处理的毫秒数）
输出：50 bps 的导航电文比特流，以及每个子帧的起始位置

    I_P（每毫秒一个值）
        │
        ▼ 位同步：找出比特边界在哪一毫秒
    比特流（每 20 毫秒一个比特，50 bps）
        │
        ▼ 帧同步：搜索同步头 10001011
    子帧边界（每 300 比特 = 6 秒一个子帧）
        │
        ▼ 下一步：星历解码
    卫星轨道参数

================================================================================
二、位同步在做什么
================================================================================
导航电文速率 50 bps，一个比特占 20 ms，正好是 20 个 C/A 码周期。
跟踪环的积分周期是 1 ms，所以我们拿到的是"每个 C/A 周期一个 I_P"，
但不知道哪 20 个连续的值属于同一个比特。

关键性质：**同一个比特内的 20 个 I_P 同号**（因为电文比特在这 20 ms 内不变）。
于是判定准则很直接 —— 遍历 20 种可能的划分起点，取让下式最大的那个：

    score(s) = Σ_k | Σ_{j=0..19} I_P[s + 20k + j] |

划分对了，20 个同号值相加是"相长"的，|和| 大；
划分错了，每个比特都被切成两半、两半异号，加起来互相抵消，|和| 小。

这个方法的优点是**不需要检测比特跳变**——弱信号下跳变检测很不靠谱，
而相长/相消的差别在低信噪比下依然存在。

================================================================================
三、帧同步在做什么
================================================================================
GPS L1 C/A 的导航电文格式（IS-GPS-200）：

    1 比特    = 20 ms
    1 字      = 30 比特
    1 子帧    = 10 字 = 300 比特 = 6 秒
    1 帧（主帧）= 5 子帧 = 1500 比特 = 30 秒

每个子帧的第 1 个字是遥测字（TLM），它的前 8 比特是固定的同步头：

     preamble = 1 0 0 0 1 0 1 1

所以帧同步就是：在比特流里搜这 8 个比特，并且要求**相邻的命中点间隔恰好
300 比特**。这个"间隔必须恰好 300"是极强的约束，随机碰撞几乎不可能满足，
所以它同时也是检验整个链路是否正确的判据。

关于极性：
    Costas 鉴相器留下 180° 相位模糊，整条比特流可能整体反号。
    所以搜索时正反两种极性都要试，并以"同一条链上极性一致"作为约束。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MS_PER_BIT = 20                                  # 50 bps → 1 比特 20 ms
BITS_PER_WORD = 30
WORDS_PER_SUBFRAME = 10
BITS_PER_SUBFRAME = BITS_PER_WORD * WORDS_PER_SUBFRAME   # 300
PREAMBLE = (1, 0, 0, 0, 1, 0, 1, 1)              # IS-GPS-200 规定的同步头
_PREAMBLE_PM1 = np.array([1 if b else -1 for b in PREAMBLE], dtype=np.int8)


@dataclass
class SyncResult:
    """一个通道的同步结果。"""
    prn: int
    bit_offset_ms: int          # 比特边界在毫秒序列中的偏移
    bits: np.ndarray            # 比特流，元素 ±1（int8）
    bit_strength: np.ndarray    # 每个比特的置信度（20 个 ms 之和的绝对值）
    subframes: list[dict]       # 每个元素含 start_bit / polarity / starts
    confidence: float = 0.0     # 位同步置信度，见 bit_sync

    @property
    def n_bits(self) -> int:
        return len(self.bits)

    def subframe_start_seconds(self, start_bit: int) -> float:
        """把子帧起始的比特序号换算成秒（相对本段数据的开头）。"""
        return start_bit * MS_PER_BIT / 1e3


def bit_sync(ip: np.ndarray, ms_per_bit: int = MS_PER_BIT) -> tuple[int, float]:
    """找出比特边界的毫秒偏移，返回 (offset, 置信度)。

    置信度 = 最佳划分得分 / 最差划分得分。参考值（实测标定）：

        ≈ 2.0   无噪声的理想比特流
        ≈ 1.7   噪声较大但仍能稳定判定
        ≈ 1.0   数据里根本没有 20 ms 的比特结构（判定不可信）

    为什么不用"最佳/均值"：正确划分与相邻划分只差 1 ms，得分几乎一样
    （20 vs 18），导致 best/mean 的理论上限只有约 1.33，区分度太差。
    而最差划分（错开 10 ms）会把每个比特劈成两半，得分掉一半，
    所以 best/min 能拉开到 2.0 左右。
    """
    ip = np.asarray(ip, dtype=np.float64)
    if len(ip) // ms_per_bit < 5:
        raise ValueError("数据太短，至少需要 5 个比特的长度才能做位同步")

    scores = np.empty(ms_per_bit)
    for s in range(ms_per_bit):
        n = (len(ip) - s) // ms_per_bit
        seg = ip[s:s + n * ms_per_bit].reshape(n, ms_per_bit)
        scores[s] = np.abs(seg.sum(axis=1)).sum()

    best = int(np.argmax(scores))
    confidence = float(scores[best] / max(scores.min(), 1e-9))
    return best, confidence


def bit_sync_transitions(ip: np.ndarray,
                         ms_per_bit: int = MS_PER_BIT) -> tuple[int, float] | None:
    """位同步的第二种独立方法：统计符号跳变位置对 20 取模的分布。

    同一个比特内的 20 个 ms 同号，比特之间的跳变只能发生在比特边界上。
    所以把所有跳变的毫秒序号对 20 取模，正确的边界会集中在一个余数上。

    返回 (offset, 置信度) 或 None（跳变太少，比如极弱的信号）。
    置信度 = 主余数计数 / 平均计数。均匀分布是 1，实测真实信号约 18（即 92%
    的跳变都落在同一个余数上）—— 这是最有区分度的一个指标。
    """
    ip = np.asarray(ip, dtype=np.float64)
    sgn = np.sign(ip)
    sgn[sgn == 0] = 1
    transitions = np.where(sgn[1:] != sgn[:-1])[0] + 1     # 新比特开始的毫秒序号
    if len(transitions) < 10:
        return None
    residues = transitions % ms_per_bit
    counts = np.bincount(residues, minlength=ms_per_bit)
    off = int(np.argmax(counts))
    return off, float(counts[off] / max(counts.mean(), 1e-9))


def extract_bits(ip: np.ndarray, offset_ms: int,
                 ms_per_bit: int = MS_PER_BIT) -> tuple[np.ndarray, np.ndarray]:
    """按给定的比特边界把每毫秒的 I_P 归并成比特流。

    返回 (bits, strength)：bits 元素取 ±1；strength 是每个比特内 20 个值之和的
    绝对值，可以当作该比特的置信度（信噪比低时它会变小）。
    """
    ip = np.asarray(ip, dtype=np.float64)
    n = (len(ip) - offset_ms) // ms_per_bit
    if n <= 0:
        raise ValueError("偏移超出数据长度")
    seg = ip[offset_ms:offset_ms + n * ms_per_bit].reshape(n, ms_per_bit)
    sums = seg.sum(axis=1)
    bits = np.where(sums >= 0, 1, -1).astype(np.int8)
    return bits, np.abs(sums)


def find_subframes(bits: np.ndarray,
                   spacing: int = BITS_PER_SUBFRAME,
                   min_chain: int = 2) -> list[dict]:
    """在比特流中搜索同步头，返回成链的子帧起始位置。

    判据：命中点之间必须间隔恰好 `spacing`（默认 300）比特，且同一条链上极性一致。
    随机 8 比特匹配的概率是 1/128，但要连续多次间隔 300 命中，概率低到可以忽略。
    """
    bits = np.asarray(bits, dtype=np.int8)
    n = len(bits)
    matches: dict[int, int] = {}
    for i in range(n - len(PREAMBLE) + 1):
        w = bits[i:i + len(PREAMBLE)]
        if np.array_equal(w, _PREAMBLE_PM1):
            matches[i] = 1
        elif np.array_equal(w, -_PREAMBLE_PM1):
            matches[i] = -1

    used: set[int] = set()
    chains: list[dict] = []
    for i in sorted(matches):
        if i in used:
            continue
        pol = matches[i]
        # 向前延伸成链
        chain = [i]
        j = i
        while (j + spacing) in matches and matches[j + spacing] == pol:
            j += spacing
            chain.append(j)
        used.update(chain)
        if len(chain) >= min_chain:
            chains.append({"start_bit": chain[0],
                           "polarity": pol,
                           "n_subframes": len(chain),
                           "starts": chain})

    return sorted(chains, key=lambda c: -c["n_subframes"])


def sync_channel(ip: np.ndarray, prn: int, ms_per_bit: int = MS_PER_BIT) -> SyncResult:
    """对一个通道做完整的位同步 + 帧同步。"""
    offset, confidence = bit_sync(ip, ms_per_bit)
    bits, strength = extract_bits(ip, offset, ms_per_bit)
    subframes = find_subframes(bits)
    return SyncResult(prn=prn, bit_offset_ms=offset, bits=bits,
                      bit_strength=strength, subframes=subframes,
                      confidence=confidence)
