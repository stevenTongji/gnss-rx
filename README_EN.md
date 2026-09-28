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
| ③ Sync | Bit sync, frame sync, navigation message | `src/gnssrx/nav_msg.py` | ⬜ |
| ④ Ephemeris | Ephemeris decoding, satellite position | `src/gnssrx/ephemeris.py` | ⬜ |
| ⑤ PVT | Pseudorange + least-squares positioning | `src/gnssrx/pvt.py` | ⬜ |
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

Two data-specific traps, documented in `scripts/02_real_data_test.py`:

- **Bandpass-sampling aliasing**: the analog IF (4.092 MHz) exceeds fs/2, so the digital IF becomes 1.364 MHz and the Doppler sign flips.
- **Bit packing order** of the 1-bit data: LSB-first for this dataset; a wrong order reduces the correlation peak to about a quarter.

## Quick start

```bash
git clone git@github.com:stevenTongji/gnss-rx.git
cd gnss-rx
uv sync                                      # requires Python >= 3.11

uv run python scripts/00_smoke_test.py       # acquisition check
uv run python scripts/01_tracking_test.py    # tracking check (~4 s)
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
