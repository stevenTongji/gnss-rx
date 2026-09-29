# gnss-rx

[中文](./README.md) | **English**

A GPS L1 C/A software receiver written from scratch. No existing GNSS library is used — from C/A code generation, parallel code-phase acquisition, and DLL/PLL tracking loops all the way to bit synchronization, ephemeris decoding and least-squares positioning.

![Tracking loop verification](docs/figures/tracking.png)

## Motivation

Preparation for GNSS and wireless communication algorithm roles. The capabilities that show up repeatedly in job descriptions — acquisition and tracking, loop parameter trade-offs, RTK/PPP, state estimation, C++ engineering — are implemented here as runnable and verifiable code.

Design principle: **all stages share one data chain and one repository**, with each stage feeding the next. One line of work yields three layers of depth on a résumé, and one coherent story in an interview.

## Status

| Stage | Module | File | State |
|---|---|---|---|
| Base | C/A code generation (PRN 1–32) | `src/gnssrx/ca_code.py` | ✅ |
| Base | IF signal simulation / I/O | `src/gnssrx/sim.py`, `io_if.py` | ✅ |
| ① Acquisition | Parallel code-phase FFT + squaring-loop fine Doppler | `src/gnssrx/acquisition.py` | ✅ |
| ② Tracking | DLL + PLL + carrier aiding | `src/gnssrx/tracking.py` | ✅ |
| ③ Sync | Bit sync, frame sync, navigation message | `src/gnssrx/nav_msg.py` | ✅ |
| ④ Ephemeris | Ephemeris decoding, satellite position (Sagnac / relativistic) | `src/gnssrx/ephemeris.py` | ✅ |
| ⑤ PVT | Pseudorange + least-squares SPS (≥4 sats) + 3-sat earth constraint | `src/gnssrx/pvt.py` | ✅ |
| ⑥ C++ port | Core modules with Eigen + CMake | `cpp/` | ⬜ |

## Verified results

Measured against **synthetic data with known ground truth** (`scripts/01_tracking_test.py`, reproducible).

| Metric | Result | Note |
|---|---|---|
| Doppler estimation residual | **0.0 Hz** | after squaring-loop refinement from ±125 Hz |
| C/N₀ estimate | 44.6 dB-Hz vs 45.0 truth | absolute error **−0.37 dB** |
| Carrier discriminator jitter | 8.2°–10.8° | measurement noise floor ≈ 7.2° |
| Smoothed true phase jitter | **1.31°** | theory 1.6° |
| Code loop lock retention | \|I_P\| 103%–125% | code drifts ≈ 2 chips over 1 s |
| Acquisition sensitivity | ≈ 37 dB-Hz | 5 ms integration at 45 dB-Hz, 20 ms at 37 dB-Hz |

Controlled experiment on the same satellite: disabling carrier aiding drops C/N₀ from 44.2 to 25.4 dB-Hz,
confirming the code loop genuinely tracks code Doppler instead of relying on the carrier loop.

![Acquisition results](docs/figures/acquisition.png)

## Real-signal validation

The same code runs on the public **Nottingham GPS L1 dataset** (1-bit, fs 5.456 MHz):

![Real-signal validation](docs/figures/real_data.png)

| Metric | Result |
|---|---|
| Satellites detected | 10 (full 32-PRN search) |
| Stable locks | **9 / 10** |
| Strongest C/N₀ | 47.2 dB-Hz (PRN 30) |
| Constellation | navigation message ±1 transitions clearly visible |

Two data-specific traps, documented in `scripts/02_real_data_test.py`, `scripts/03_real_data_test.py`:

- **Bandpass-sampling aliasing**: the analog IF (4.092 MHz) exceeds fs/2, so the digital IF becomes 1.364 MHz and the Doppler sign flips.
- **Bit packing order** of the 1-bit data: LSB-first for this dataset; a wrong order reduces the correlation peak to about a quarter.

### Navigation message sync (bit sync → frame sync)

On 30 s of real signal, all three satellites synchronize:

| Metric | Result |
|---|---|
| Bit sync | Two independent methods **agree exactly**, 1 ms resolution |
| Transition confidence | 11.5–20.0 (uniform would be 1; 92% of bit transitions fall on one residue) |
| Frame sync | **5 subframes** per satellite, spaced exactly 300 bits = 6.000 s |
| Subframe start times | 4.6 / 10.6 / 16.6 / 22.6 / 28.6 s (identical across satellites) |

Requiring an exact 300-bit spacing is a very hard constraint — a random 8-bit match has
probability 1/128, and hitting that repeatedly at exactly 300-bit intervals is essentially
impossible, so this simultaneously validates acquisition, tracking and bit sync.

![Navigation message sync](docs/figures/nav_msg.png)

### Single-point positioning (PVT: real data + synthetic validation)

**Real data (SiGe GN3S v3, 8-bit, fs 16.368 MHz)**: the full chain — acquisition → 30 s tracking
→ ephemeris decode → pseudorange → weighted least squares — runs end to end. After fixing the parity
bug (below), this dataset (ION GNSS SDR metadata-standard sample, collected 2013-05-23) yields
**7 clean satellites** (PRN 1/7/8/9/11/17/28, C/N₀ 45–50 dB-Hz), which meets the full-SPS (≥4 sats)
requirement; the Nottingham 1-bit dataset likewise yields **9 clean satellites**
(`scripts/14_nottingham_pvt.py`).

> **Parity bug fixed (`ephemeris.py::strip_parity`)**: the IS-GPS-200 parity rule is
> `d_i = r_i ⊕ D30*(prev word)`, where D30\* is the *received* bit (including the Costas inversion c),
> so correctly `d_i = r_i ⊕ r_prevD30` (the c cancels). The old code wrote `data ^ D30* ^ c`, which is
> right for the first word but adds a spurious `c` on every later word: c=0 channels happened to be
> correct while c=1 channels were fully inverted → a "garbage orbit". With 1-bit data dominated by c=1,
> many real satellites were mislabelled as C/A cross-correlation ghosts. Fixing it restored SiGe from
> 3 → **7** and Nottingham from 2 → **9** clean satellites.

> **Two sign bugs fixed** — the key to going from "hundreds of km" to metre level, both pinned down by
> the synthetic self-check `scripts/15_timing_selfcheck.py`:
> (1) **code-phase sign**: `t_user = m_blk/1000 − τ/1.023e6` (was `+`); this tracker's τ has the
> opposite sign to the code-epoch offset within the block; and
> (2) **satellite-clock sign**: `pseudorange = c(t_user − t_sat) + c·δt_sat` (was `−`), since
> `c(t_user − t_sat) = ρ + c·b − c·δt_sat`.
> Before the fix Nottingham was ~298 km off; after it, the 9-satellite real-data fix lands at
> **52.93°N / 1.17°W / h ≈ 0** with **±2 m residuals**, i.e. **1.96 km** from the published collection
> point — metre-level real-data positioning.

**Synthetic 4-sat validation (`scripts/13_synthetic_4sat_pvt.py`, known truth)**: because the real
dataset has only 3 sats, a known-truth scenario exercises the full path (true range → injected
integer-millisecond ambiguity → `build_measurement` → least squares) to prove the engine correct:

| Scenario | Position error | Note |
|---|---|---|
| 6 sats, no ambiguity, direct WLS | **0.024 m** | solver math correct |
| 6 sats + injected common ms ambiguity (K=434579) | **0.024 m** | integer-ms disambiguation correct (GPS-week scale) |
| 4-sat subset | **0.024 m** | full SPS correct |
| code-phase noise σ = 2 / 5 / 10 m | 6.3 / 15.8 / 31.5 m | error ≈ σ·GDOP, as predicted |
| 3 sats + earth constraint | ~11 km | two-point ambiguity — motivates ≥4 sats |

**Fix**: the GPS week number is only 10 bits in the message (WN mod 1024), rolling over every 1024
weeks. This data decodes 717 but the true extended week is 1741 (2013); `ephemeris.py` now restores it
to the nearest 1024-multiple of a reference week. (The rollover does not affect satellite position,
which depends only on the within-week tk, but the corrected date is now consistent.)

## Quick start

```bash
git clone git@github.com:stevenTongji/gnss-rx.git
cd gnss-rx
uv sync                                      # requires Python >= 3.11

uv run python scripts/00_smoke_test.py       # acquisition check
uv run python scripts/01_tracking_test.py    # tracking check (~4 s)
uv run python scripts/02_real_data_test.py   # real-signal check (needs data)
uv run python scripts/03_nav_msg_test.py     # nav message sync (needs data)
uv run python scripts/13_synthetic_4sat_pvt.py  # synthetic 4-sat PVT verification (known truth)
uv run python scripts/15_timing_selfcheck.py # timing self-check (synthetic, pins the sign bugs)
uv run python scripts/07_sige_pvt_test.py    # real-data PVT on SiGe (needs data)
uv run python scripts/14_nottingham_pvt.py   # real-data PVT on Nottingham (needs data)
```

## Implementation notes

Details and derivations live in the source comments; key points:

- **Acquisition**: one FFT covers all 1023 code phases, with a serial Doppler search; detection uses peak-to-second-peak ratio. At 45 dB-Hz, 1 ms coherent integration leaves too little margin — 5 ms non-coherent accumulation is used.
- **Fine Doppler**: code wipe-off → per-ms complex correlation → squaring removes the navigation message → zero-padded FFT. Essential, since a second-order PLL cannot pull in hundreds of hertz on its own.
- **Code loop**: normalized non-coherent early-minus-late discriminator, gain −2 /chip, linear range ±(1 − d/2) chips.
- **Carrier loop**: Costas discriminator `atan(Q/I)`, insensitive to 180° navigation bit flips, at the cost of a 180° phase ambiguity.
- **Loop filter**: second-order PI with Kaplan's coefficients ω_n = 8ζB_n/(4ζ²+1), τ₂ = 2ζ/ω_n, τ₁ = 1/ω_n². Discriminator outputs are normalized to the controlled quantity itself, so total loop gain is 1.
- **Carrier aiding**: code rate derived from carrier Doppler, f_code = 1.023 MHz × (1 + f_d / f_L1); the code loop handles only the residual. A 3000 Hz Doppler shifts the code rate by 1.95 chips/s, which a 2 Hz unaided code loop cannot follow.

## References

- K. Borre, D. Akos, et al., *A Software-Defined GPS and Galileo Receiver: A Single-Frequency Approach*, Birkhäuser, 2007.
- E. D. Kaplan, C. Hegarty, *Understanding GPS/GNSS: Principles and Applications*, 3rd ed., Artech House, 2017.
- IS-GPS-200, *NAVSTAR GPS Space Segment / User Segment Interfaces*.
- T. D. Barfoot, *State Estimation for Robotics*, Cambridge, 2017.
- 3GPP TS 38.305 (NR positioning), TS 38.821 (NTN).

## License

MIT
