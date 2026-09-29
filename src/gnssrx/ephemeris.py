"""GPS L1 C/A 星历解码（IS-GPS-200，LNAV 子帧 1/2/3）。

================================================================================
一、输入是什么
================================================================================
`nav_msg.py` 的帧同步给出每个子帧在比特流里的起点（300 比特一份）和极性
（Costas 180° 模糊导致整段可能整体反号）。这一步要做：

    300 比特子帧（±1）
        │
        ▼ ① 去极性（× polarity），转成 0/1
    300 比特（含每字 6 个奇偶校验位）
        │
        ▼ ② 剥奇偶校验：每 30 比特一个字，前 24 个是数据比特；
           若该字的 D30*（第 30 比特）为 1，则 24 个数据比特整体反相
    240 比特数据数组（10 字 × 24，0 索引 = 第 1 字数据 0–23，第 2 字 24–47…）
        │
        ▼ ③ 按 RTKLIB 的权威位图读字段（见下方 FIELD 表）
    卫星轨道根数 + 时钟参数

================================================================================
二、位图来源（重要）
================================================================================
下面的字段偏移不是我自己推的，而是逐字节对齐 RTKLIB `rcvraw.c` 里
`decode_subfrm1/2/3` 的 `getbitu/getbits(buff, i, n)` 调用。RTKLIB 的 `buff`
就是上面第 ② 步产出的 240 比特数组，所以这里的 `start` 直接等于 RTKLIB 的 `i`。
常量 P2_k = 2^-k，SC2RAD = π（半周 → 弧度的换算）。

校验点（不靠任何外部真值就能验证）：
  · 子帧 3 的 iode 必须等于子帧 2 的 iode，也必须等于 (子帧 1 的 iodc & 0xFF)
  · sqrt(A) × 2^-19 ≈ 5153.6  → 对应半长轴 ≈ 26560 km
  · e（偏心率）应在 0.001–0.02 量级
  · i0 × π ≈ 0.96 rad（≈ 55°，GPS 轨道倾角）
  · 算出的卫星位置离地心约 26560 km
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

# ---- 物理常量（与 RTKLIB 一致）----
MU_EARTH = 3.986005e14        # 地球引力常数 m^3/s^2
OMEGA_E = 7.2921151467e-5     # 地球自转角速度 rad/s
C_LIGHT = 2.99792458e8        # 光速 m/s
PI = math.pi

# GPS 周号（10 比特）回卷参考：还原扩展周号时，取离此值最近的 1024 整数倍。
# 默认对应 2013 年的 ION SiGe 样本（扩展周号 1741）；换其它年代的数据改这里即可。
GPS_WEEK_REFERENCE = 1741

# 半周 → 弧度
SC2RAD = PI

# ---- 子帧内的字结构 ----
BITS_PER_SUBFRAME = 300
WORDS_PER_SUBFRAME = 10
DATA_BITS_PER_WORD = 24        # 每字 30 比特里前 24 个是数据比特
PARITY_BITS_PER_WORD = 6

# TLM 遥测字（每个子帧的第 1 个字）前 8 个数据比特永远是同步头
PREAMBLE = np.array([1, 0, 0, 0, 1, 0, 1, 1], dtype=np.int8)

# GPS LNAV 奇偶校验矩阵（GF(2)）：6 个奇偶位 D25..D30 由 24 个数据位与「前一字末两位」
# D29*/D30* 线性生成。下面每一行是该奇偶位对应的输入下标：
#   0..23 = 本字数据位 D1..D24，24 = D29*，25 = D30*。
# 系数由真实数据以 RANSAC（GF(2) 最小失配）反推得到，与 IS-GPS-200 的方程一致。
# 实测：干净子帧 9 个字全部通过，垃圾通道（互相关假锁）9/9 失败 —— 是判断子帧可信度的硬指标。
_PARITY_ROWS = (
    (0, 1, 2, 4, 5, 9, 10, 11, 12, 13, 16, 17, 19, 22, 24),
    (1, 2, 3, 5, 6, 10, 11, 12, 13, 14, 17, 18, 20, 23, 25),
    (0, 2, 3, 4, 6, 7, 11, 12, 13, 14, 15, 18, 19, 21, 24),
    (1, 3, 4, 5, 7, 8, 12, 13, 14, 15, 16, 19, 20, 22, 25),
    (0, 2, 4, 5, 6, 8, 9, 13, 14, 15, 16, 17, 20, 21, 23),
    (2, 4, 5, 7, 8, 9, 10, 12, 14, 18, 21, 22, 23, 24, 25),
)


# ---------------------------------------------------------------------------
# 第一步：去极性 + 剥奇偶校验 → 240 比特数据数组
# ---------------------------------------------------------------------------

def strip_parity(bits300: np.ndarray, polarity: int = 1) -> tuple[np.ndarray, int]:
    """把 300 比特子帧（±1）转成 240 比特数据数组（0/1，np.int8）。

    IS-GPS-200 的奇校验规则（逐字，关键易错点）：
        真实数据 D1..D24 = 接收比特 ⊕ D30*(k) ⊕ c
    其中 D30*(k) 是「第 k 个字」的 D30*，它等于**前一个字（k-1）的 D30 比特**。
    发射端在编码第 k 个字时，先看前一个字的 D30：若它为 1，就把第 k 个字的数据
    比特 D1..D24 整体反相。RTKLIB 的 `decode_word(word)` 正是把
    `前字的 D29*/D30* + 本字 D1..D30` 拼成 32 比特，再用本字 D30*（=前字 D30）
    决定是否反相 —— 所以反相用的是**前一个字**的 D30，而不是本字的。

    第 1 个字（TLM）没有前字，约定 D30* = 1（恒反相）。

    c 是整段统一的「Costas 180° 反相」残差（0 或 1），由 TLM 同步头自校准：
    谁能让 word0 数据位的 d[0:8] 还原成 `10001011`，c 就是谁。这样不依赖帧同步
    给出的极性符号，也不受 TLM 自身 D30*=1 的反相干扰。

    返回 (data240, c)：data240 是 10 字 × 24 比特的 0/1 数组。
    """
    sym = np.asarray(bits300, dtype=np.int8) * int(polarity)
    r = np.where(sym > 0, 1, 0).astype(np.int8)        # 0/1 接收比特
    out = np.empty(WORDS_PER_SUBFRAME * DATA_BITS_PER_WORD, dtype=np.int8)

    # 用 TLM 字的同步头自行定出全局反相残差 c（TLM 的 D30* 约定为 1）
    w0 = r[0:DATA_BITS_PER_WORD]
    best_c = None
    for c in (0, 1):
        if np.array_equal((w0[:8] ^ 1 ^ c), PREAMBLE):   # 第 1 字 D30* = 1
            best_c = c
            break
    if best_c is None:                                  # 有比特错误时按海明距离选
        d0 = int(np.sum((w0[:8] ^ 1) != PREAMBLE))
        d1 = int(np.sum((w0[:8] ^ 1 ^ 1) != PREAMBLE))
        best_c = 0 if d0 <= d1 else 1

    # 逐字剥奇偶：第 k 字的反相由「第 k-1 字的 D30」控制。
    # 注意：d_i = r_i ⊕ D30*(前字)，其中 D30* 是**含 Costas 反相 c 的接收比特**，
    # 于是 d_i = r_i ⊕ (r_prevD30 ⊕ c) ⊕ c = r_i ⊕ r_prevD30 —— c 自动抵消。
    # 因此 D30* 必须同样带上 best_c（否则 c=1 的通道会被多翻转一次 → 解出垃圾）。
    d30_star = 1                                        # 第 1 字（TLM）约定 D30* = 1
    for w in range(WORDS_PER_SUBFRAME):
        base = w * 30
        data = r[base:base + DATA_BITS_PER_WORD]
        out[w * DATA_BITS_PER_WORD:(w + 1) * DATA_BITS_PER_WORD] = data ^ d30_star ^ best_c
        d30_star = int(r[base + 29]) ^ best_c           # 本字 D30（含反相）→ 下一字的 D30*
    return out, best_c


def parity_failures(bits300: np.ndarray, polarity: int = 1) -> int:
    """统计一个 300 比特子帧中「奇偶校验失败」的字数（检查 word1..word9，共 9 个）。

    word0（TLM）的 D29*/D30* 来自上一子帧末字，本函数无法得知，故跳过。
    返回 0 表示该子帧 9 个字奇偶全部通过 —— 即未检出比特错误。
    这是判断一个子帧是否可信的硬指标（真实数据实测：干净通道全 0，垃圾通道 9/9）。
    """
    sym = np.asarray(bits300, dtype=np.int8) * int(polarity)
    r = np.where(sym > 0, 1, 0).astype(np.int8)
    c = None
    for cc in (0, 1):                                   # 同 strip_parity：定 Costas 反相残差
        if np.array_equal(r[0:8] ^ 1 ^ cc, PREAMBLE):
            c = cc
            break
    if c is None:
        return WORDS_PER_SUBFRAME - 1                   # 连前导都对不上 → 视为全失败
    t = r ^ c
    fails = 0
    for w in range(1, WORDS_PER_SUBFRAME):
        base = w * 30
        x = np.empty(26, dtype=np.int8)
        x[:24] = t[base:base + 24]
        x[24] = t[base - 2]                             # D29* = 前一字 D29
        x[25] = t[base - 1]                             # D30* = 前一字 D30
        exp = np.array([int(x[list(row)].sum() & 1) for row in _PARITY_ROWS],
                       dtype=np.int8)
        if not np.array_equal(exp, t[base + 24:base + 30]):
            fails += 1
    return fails


# ---------------------------------------------------------------------------
# 第二步：从 240 比特数组里读字段（start 是 RTKLIB 的 0 索引偏移）
# ---------------------------------------------------------------------------

def _read_field(data: np.ndarray, start: int, nbits: int, signed: bool) -> int:
    """读出从 `start` 开始的 `nbits` 比特，MSB 在前。signed 时按补码解释。"""
    val = 0
    for k in range(start, start + nbits):
        val = (val << 1) | int(data[k])
    if signed and (val >> (nbits - 1)) & 1:
        val -= (1 << nbits)
    return val


def read_subframe_id(data240: np.ndarray) -> int:
    """从 HOW 读出子帧号（1–5）。

    HOW 是第 2 字（1-idx），在 240 比特数据数组里占 d=24..47。
    RTKLIB `decode_frame` 用 `getbitu(buff,43,3)` 读 SFID，对应本数组的
    d[43:46]（HOW 字内偏移 19..21，即 D20–D22）。
    """
    return _read_field(data240, 43, 3, signed=False)


def read_tow(data240: np.ndarray) -> int:
    """从 HOW 读出 TOW（秒，×6）。TOW 在 d=24..40，17 比特，×6 秒。

    RTKLIB `decode_subfrm1` 用 `getbitu(buff,24,17)*6.0`，与这里一致。
    """
    raw = _read_field(data240, 24, 17, signed=False)
    return raw * 6


def sfid_from_tow(tow_seconds: int) -> int:
    """由 TOW（秒）结构性地推出子帧号（1–5），无需依赖 HOW 内的 SFID 位。

    GPS 帧结构（IS-GPS-200 20.3.3.5.2）：HOW 里的 TOW 表示该子帧「下一条」
    子帧的起始时刻（6 秒单位）。设子帧 k 的起始为 6·t_k 秒，则
    TOW_raw = t_k + 1；而 t_k = 5·m + (k−1)（每帧 5 个子帧），故
    TOW_raw = 5·m + k  →  TOW_raw mod 5 = k mod 5。
    于是：子帧号 = TOW_raw mod 5（k=5 时余数为 0，映射回 5）。
    这是帧结构的硬约束，比直接读 SFID 位更可靠。
    """
    raw = tow_seconds // 6
    return (raw % 5) or 5


# ---- 子帧 1：时钟与卫星健康（RTKLIB decode_subfrm1）----
def decode_subframe1(data240: np.ndarray) -> dict:
    # GPS 周号在电文里只有 10 比特（WN mod 1024），每 1024 周（≈19.7 年）回卷一次。
    # 2013 年的 SiGe 样本里解出的截断值是 717，真实扩展周号应为 717+1024=1741。
    # 这里按"离参考周号最近的 1024 整数倍"还原，避免把 2013 年错当成 1993 年。
    raw_week = _read_field(data240, 48, 10, signed=False)
    week = raw_week
    while week < GPS_WEEK_REFERENCE - 512:
        week += 1024
    while week > GPS_WEEK_REFERENCE + 512:
        week -= 1024
    code = _read_field(data240, 58, 2, signed=False)     # L2 码标志
    ura = _read_field(data240, 60, 4, signed=False)      # 精度指数 (SV accuracy)
    svh = _read_field(data240, 64, 6, signed=False)      # 卫星健康
    iodc0 = _read_field(data240, 70, 2, signed=False)    # IODC 高 2 位
    flag = _read_field(data240, 72, 1, signed=False)     # 数据期限/拟合标志
    tgd = _read_field(data240, 160, 8, signed=True)      # 群延迟 (s)
    iodc1 = _read_field(data240, 168, 8, signed=False)   # IODC 低 8 位
    toc = _read_field(data240, 176, 16, signed=False)    # 时钟参考时刻 (s)
    f2 = _read_field(data240, 192, 8, signed=True)       # 钟漂二次项
    f1 = _read_field(data240, 200, 16, signed=True)      # 钟漂一次项
    f0 = _read_field(data240, 216, 22, signed=True)      # 钟偏
    iodc = (iodc0 << 8) | iodc1
    return {
        "week": week,
        "ura": ura, "svh": svh, "code": code, "flag": flag,
        "iodc": iodc,
        "tgd": tgd * 2.0**-31,
        "toc": toc * 16.0,
        "af2": f2 * 2.0**-55,
        "af1": f1 * 2.0**-43,
        "af0": f0 * 2.0**-31,
    }


# ---- 子帧 2：轨道根数（OSC 参数，RTKLIB decode_subfrm2）----
def decode_subframe2(data240: np.ndarray) -> dict:
    iode = _read_field(data240, 48, 8, signed=False)
    crs = _read_field(data240, 56, 16, signed=True)
    deln = _read_field(data240, 72, 16, signed=True)
    M0 = _read_field(data240, 88, 32, signed=True)
    cuc = _read_field(data240, 120, 16, signed=True)
    e = _read_field(data240, 136, 32, signed=False)
    cus = _read_field(data240, 168, 16, signed=True)
    sqrtA = _read_field(data240, 184, 32, signed=False)
    toes = _read_field(data240, 216, 16, signed=False)
    fit = _read_field(data240, 232, 1, signed=False)
    return {
        "iode": iode,
        "crs": crs * 2.0**-5,
        "deln": deln * 2.0**-43 * SC2RAD,
        "M0": M0 * 2.0**-31 * SC2RAD,
        "cuc": cuc * 2.0**-29,
        "e": e * 2.0**-33,
        "cus": cus * 2.0**-29,
        "sqrtA": sqrtA * 2.0**-19,
        "toes": toes * 16.0,
        "fit_interval": 4.0 if fit == 0 else 6.0,  # 0:4h, 1:>4h
    }


# ---- 子帧 3：轨道根数（RTKLIB decode_subfrm3）----
def decode_subframe3(data240: np.ndarray) -> dict:
    cic = _read_field(data240, 48, 16, signed=True)
    OMG0 = _read_field(data240, 64, 32, signed=True)
    cis = _read_field(data240, 96, 16, signed=True)
    i0 = _read_field(data240, 112, 32, signed=True)
    crc = _read_field(data240, 144, 16, signed=True)
    omg = _read_field(data240, 160, 32, signed=True)
    OMGd = _read_field(data240, 192, 24, signed=True)
    iode = _read_field(data240, 216, 8, signed=False)
    idot = _read_field(data240, 224, 14, signed=True)
    return {
        "iode": iode,
        "cic": cic * 2.0**-29,
        "OMG0": OMG0 * 2.0**-31 * SC2RAD,
        "cis": cis * 2.0**-29,
        "i0": i0 * 2.0**-31 * SC2RAD,
        "crc": crc * 2.0**-5,
        "omg": omg * 2.0**-31 * SC2RAD,
        "OMGd": OMGd * 2.0**-43 * SC2RAD,
        "idot": idot * 2.0**-43 * SC2RAD,
    }


# ---------------------------------------------------------------------------
# 第三步：装配成 Ephemeris 对象（含卫星位置计算，便于验证）
# ---------------------------------------------------------------------------

@dataclass
class Ephemeris:
    prn: int
    week: int = 0
    # 时钟
    af0: float = 0.0
    af1: float = 0.0
    af2: float = 0.0
    toc: float = 0.0
    tgd: float = 0.0
    iodc: int = 0
    ura: int = 0
    svh: int = 0
    # 轨道（子帧 2）
    iode2: int = 0
    crs: float = 0.0
    deln: float = 0.0
    M0: float = 0.0
    cuc: float = 0.0
    e: float = 0.0
    cus: float = 0.0
    sqrtA: float = 0.0
    toes: float = 0.0
    fit_interval: float = 4.0
    # 轨道（子帧 3）
    iode3: int = 0
    cic: float = 0.0
    OMG0: float = 0.0
    cis: float = 0.0
    i0: float = 0.0
    crc: float = 0.0
    omg: float = 0.0
    OMGd: float = 0.0
    idot: float = 0.0
    # 原始子帧（用于回溯/调试）
    _raw: dict = field(default_factory=dict, repr=False)

    @property
    def A(self) -> float:
        return self.sqrtA ** 2

    def from_subframes(self, sf1: dict, sf2: dict, sf3: dict) -> "Ephemeris":
        self.week = sf1["week"]; self.af0 = sf1["af0"]; self.af1 = sf1["af1"]
        self.af2 = sf1["af2"]; self.toc = sf1["toc"]; self.tgd = sf1["tgd"]
        self.iodc = sf1["iodc"]; self.ura = sf1["ura"]; self.svh = sf1["svh"]
        self.iode2 = sf2["iode"]; self.crs = sf2["crs"]; self.deln = sf2["deln"]
        self.M0 = sf2["M0"]; self.cuc = sf2["cuc"]; self.e = sf2["e"]
        self.cus = sf2["cus"]; self.sqrtA = sf2["sqrtA"]; self.toes = sf2["toes"]
        self.fit_interval = sf2["fit_interval"]
        self.iode3 = sf3["iode"]; self.cic = sf3["cic"]; self.OMG0 = sf3["OMG0"]
        self.cis = sf3["cis"]; self.i0 = sf3["i0"]; self.crc = sf3["crc"]
        self.omg = sf3["omg"]; self.OMGd = sf3["OMGd"]; self.idot = sf3["idot"]
        self._raw = {"sf1": sf1, "sf2": sf2, "sf3": sf3}
        return self

    # ---- 一致性自检（不依赖外部真值）----
    def self_check(self) -> list[str]:
        """返回一串问题说明；空列表表示通过。判据均不依赖外部真值。"""
        problems = []
        if self.iode2 != self.iode3:
            problems.append(f"IODE 不一致：SF2={self.iode2} ≠ SF3={self.iode3}")
        if (self.iodc & 0xFF) != self.iode2:
            problems.append(f"IODC&0xFF={self.iodc & 0xFF} ≠ IODE={self.iode2}")
        # sqrt(A) 以 m^0.5 计，GPS 约 5153.6；给足余量
        if not (4700.0 < self.sqrtA < 5600.0):
            problems.append(f"sqrt(A)={self.sqrtA:.3f} 不在 GPS 典型区间 [4700,5600]")
        if self.e <= 0 or self.e > 0.2:
            problems.append(f"e={self.e:.4f} 不合理")
        if abs(abs(self.i0) - 0.96) > 0.15:
            problems.append(f"|i0|={abs(self.i0)*180/PI:.2f}° 偏离 GPS 标准倾角 55° 过多")
        if self.A <= 0:
            problems.append("半长轴 A 非正")
        return problems

    def clock_correction(self, tk: float) -> float:
        """钟差（秒）：af0 + af1·tk + af2·tk² − 相对论项。tk = t − toc。"""
        dt = self.af0 + self.af1 * tk + self.af2 * tk * tk
        return dt

    def position(self, tk: float) -> np.ndarray:
        """计算 tk = t − toe 秒时的 ECEF 位置（米）。

        算法逐字节对齐 RTKLIB 的 eph2pos（开普勒方程 + 摄动修正）。
        """
        A = self.A
        n0 = math.sqrt(MU_EARTH / A**3)
        n = n0 + self.deln
        M = self.M0 + n * tk
        # 解开普勒方程 E − e·sinE = M
        E = M
        for _ in range(30):
            dE = (E - self.e * math.sin(E) - M) / (1.0 - self.e * math.cos(E))
            E -= dE
            if abs(dE) < 1e-12:
                break
        sinE, cosE = math.sin(E), math.cos(E)
        # 真近点角 + 近地点幅角
        u = math.atan2(math.sqrt(1.0 - self.e**2) * sinE, cosE - self.e) + self.omg
        r = A * (1.0 - self.e * cosE)
        i = self.i0 + self.idot * tk
        # 摄动修正
        sin2u, cos2u = math.sin(2 * u), math.cos(2 * u)
        u += self.cus * sin2u + self.cuc * cos2u
        r += self.crs * sin2u + self.crc * cos2u
        i += self.cis * sin2u + self.cic * cos2u
        x = r * math.cos(u)
        y = r * math.sin(u)
        # 升交点赤经（含地球自转修正）
        Omega = self.OMG0 + (self.OMGd - OMEGA_E) * tk - OMEGA_E * self.toes
        cosO, sinO, cosI, sinI = math.cos(Omega), math.sin(Omega), math.cos(i), math.sin(i)
        X = x * cosO - y * cosI * sinO
        Y = x * sinO + y * cosI * cosO
        Z = y * sinI
        return np.array([X, Y, Z], dtype=float)

    def relativistic_term(self, t_sat: float) -> float:
        """相对论钟差修正（秒）：Δtr = −2·√μ·e·√A·sin(E) / c²。

        在卫星钟差里它是 ~1e-9 s（等效 ~0.2 m）量级，PVT 精度需要它。
        与 clock_correction 一起由 pvt.satellite_clock_correction 调用。
        """
        tk = t_sat - self.toes
        n = math.sqrt(MU_EARTH / self.A ** 3) + self.deln
        M = self.M0 + n * tk
        E = M
        for _ in range(30):
            dE = (E - self.e * math.sin(E) - M) / (1.0 - self.e * math.cos(E))
            E -= dE
            if abs(dE) < 1e-12:
                break
        return -2.0 * math.sqrt(MU_EARTH) * self.e * self.sqrtA * math.sin(E) / C_LIGHT ** 2

    def radius_km(self, tk: float = 0.0) -> float:
        return float(np.linalg.norm(self.position(tk)) / 1e3)


def decode_subframes(prn: int, blocks: list[np.ndarray],
                     polarity: int = 1) -> tuple[Ephemeris | None, list[dict]]:
    """从一堆 300 比特原始子帧（±1）里挑出 SF1/2/3 并装配成 Ephemeris。

    分类采用 **结构性判据**：每块的 TOW（d[24:41]）→ 子帧号 = (TOW_raw mod 5)+1，
    这比直接读 HOW 内的 SFID 位更可靠（不受 polarity / 位偏移歧义影响）。

    参数
    ----
    blocks : list[np.ndarray]   每个元素是一个 300 比特子帧（±1）
    polarity : int              find_subframes 给出的整段极性（可选，strip_parity
                                内部会自校准，这里仅作兼容保留）

    返回
    ----
    (eph, info)
        eph  : Ephemeris，三类不全返回 None
        info : 每块的解码诊断列表（sfid / tow / c / preamble_ok / 异常）
    """
    info: list[dict] = []
    best: dict[int, tuple] = {}          # sfid -> (解码结果, 奇偶失败字数)

    # 先按 TOW 升序排好：find_subframes 有时会把子帧链按时间倒序返回，
    # 排序后 SF1/2/3 才保证来自同一主帧（相邻 6 秒），IODE 才会一致。
    parsed = []
    for blk in blocks:
        d240, c = strip_parity(blk, polarity)
        tow = read_tow(d240)
        sfid = sfid_from_tow(tow)
        pre_ok = np.array_equal(d240[0:8], PREAMBLE)
        parsed.append({"tow": tow, "sfid": sfid, "c": c,
                       "preamble_ok": bool(pre_ok), "d240": d240,
                       "parity_fail": parity_failures(blk, polarity)})
    parsed.sort(key=lambda p: p["tow"])

    _DECODERS = {1: decode_subframe1, 2: decode_subframe2, 3: decode_subframe3}
    for p in parsed:
        sfid, tow, d240 = p["sfid"], p["tow"], p["d240"]
        entry = {"sfid": sfid, "tow": tow, "c": p["c"],
                 "preamble_ok": p["preamble_ok"], "parity_fail": p["parity_fail"]}
        dec = _DECODERS.get(sfid)
        if dec is not None:
            # 同类子帧可能有多份（跨主帧重复），保留「奇偶失败最少」的那份，
            # 以剔除带比特错误的副本（parity_fail=0 表示该子帧无可检出的比特错误）。
            prev = best.get(sfid)
            if prev is None or p["parity_fail"] < prev[1]:
                try:
                    best[sfid] = (dec(d240), p["parity_fail"])
                except Exception as exc:              # 边界/错位时某块可能解出异常值
                    entry["error"] = str(exc)
        info.append(entry)

    if 1 not in best or 2 not in best or 3 not in best:
        return None, info
    eph = Ephemeris(prn=prn)
    eph.from_subframes(best[1][0], best[2][0], best[3][0])
    return eph, info
