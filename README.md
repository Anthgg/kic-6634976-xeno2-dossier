# Cured Multi-Wavelength Dossier: KIC 6634976 (`bnty_0268` — XENO-2 Revise Edition)

**Target:** `KIC 6634976` (`Gaia DR3 2077624676762285440` / `AllWISE J195516.58+420147.5` / `2MASS J19551658+4201475`)  
**Coordinates (ICRS J2000):** `RA = 298.819090 deg` (`19h 55m 16.58s`), `Dec = +42.029850 deg` (`+42d 01m 47.5s`)  
**Galactic Coordinates:** `l = 76.8407 deg`, `b = +7.1103 deg` (Cygnus ISM / Cirrus Complex)  
**Repository & Reproducible Pipeline:** [`kic6634976_dossier_pipeline.py`](./kic6634976_dossier_pipeline.py) | [`results.json`](./results.json) | [`run_output.txt`](./run_output.txt)

---

## 1. Executive Summary & Final Disposition

This repository contains the complete, zero-dependency, end-to-end Python verification pipeline (`kic6634976_dossier_pipeline.py`) and cured multi-wavelength dossier for **`KIC 6634976`**, originally flagged in `XENO-2` (`bnty_0246`) as an optically quiet main-sequence `F8V/G0V` star (`Teff = 6029 K`, `d ~ 1030 pc`) with a dramatic mid-infrared excess (`6.82x` at `12 um W3` and `75.62x` at `22 um W4`, `f_IR = L_IR / L_* = 1.74%`, `Delta_BIC = +512.46`).

Following independent replication (`bnty_0264`) and the extended multi-wavelength checks executed in this repository (adding **NVSS 1.4 GHz**, **VLASS 3.0 GHz**, **IRSA SFD98/SF11 Dust & 100 um DIRBE/IRAS**, and **60-arcsec AllWISE/unWISE PSF profile morphology**), **we formally incorporate the Galactic Cirrus Filament explanation as the definitive astrophysical resolution of `KIC 6634976`**:

1. **Optical Photometric Stability (Confirmed):** Across all `17` STScI MAST Kepler Long-Cadence FITS files (`Q1-Q17`, `69,302` raw cadences, `55,909` valid `SAP_QUALITY == 0` cadences spanning `1,459.49 days`), `KIC 6634976` exhibits a median quarterly RMS of **`192.43 ppm`** (`162.90 ppm` to `222.75 ppm`) and a maximum downward dip of **`0.1252%`** (`min_norm_pdc = 0.998748`). There are no transits, eclipses, or occultation megastructures along the line of sight.
2. **NVSS (1.4 GHz) & VLASS (3.0 GHz) Radio Counterpart Check Within 5 Arcsec (Confirmed Null):**
   - **NVSS (`VIII/65/nvss`, 1.4 GHz, `5-sigma ~ 2.5 mJy/beam`):** **`0` matches within `5.0 arcsec`** and **`0` matches within `60.0 arcsec`**. The nearest NVSS radio source is `NVSS 195518+415628` (`RA = 298.829000, Dec = +41.941361`) at an angular separation of **`319.66 arcsec`** (`5.33 arcmin`) with flux density `S_1.4 = 3.8 +/- 1.0 mJy`.
   - **VLASS (`J/ApJS/255/30/comp`, CIRADA VLASS Epoch 1 QL, 3.0 GHz, `5-sigma ~ 0.60 mJy/beam`):** **`0` matches within `5.0 arcsec`** and **`0` matches within `60.0 arcsec`**. The nearest VLASS radio component is `VLASS1QLCIR J195524.12+420133.6` (`RA = 298.850505, Dec = +42.026011`) at an angular separation of **`85.14 arcsec`** (`1.42 arcmin`) with integrated flux `F_tot = 4.083 +/- 0.211 mJy` (`F_peak = 4.354 +/- 0.126 mJy/beam`).
   - **Implication:** This rules out radio-loud background AGN contamination, compact H II regions, and gyrosynchrotron/maser emission within the WISE `W3`/`W4` beam (`< 2.5 mJy` at 1.4 GHz and `< 0.60 mJy` at 3.0 GHz).
3. **Galactic Cirrus Filament Resolution (Incorporated Head-On with Quantitative Proof):**
   - **High Line-of-Sight Dust Column & Far-IR Cirrus Brightness (`IRSA Dust Service`):** At `(298.819090, +42.029850)` (`b = +7.11 deg` in Cygnus), the Schlegel, Finkbeiner & Davis (1998) reddening is **`E(B-V)_SFD = 0.3338 mag`** (`5 deg` local mean `0.3294 +/- 0.0162 mag`; Schlafly & Finkbeiner 2011 recalibrated **`E(B-V)_SandF = 0.2871 mag`**), well above the `E(B-V) < 0.10 mag` clean-sky threshold for mid-IR technosignature sweeps. The co-spatial IRAS/DIRBE `100 um` surface brightness is **`21.3249 MJy/sr`** (`mean = 21.0427 +/- 1.0972 MJy/sr`) at an equilibrium ISM dust temperature of **`T_dust,ISM = 18.5543 K`**.
   - **Shared `W3`/`W4` Excess Floor Across Unrelated Field Stars (`60"` AllWISE Cone):** Querying `II/328/allwise` within `60.0 arcsec` recovers **`22` sources**, of which **`9` sources (`7` within `41.46 arcsec`)** are strongly detected in `W3` (`snr3 = 9.7` to `23.6`) along a narrow NE-SW filamentary ridge. Even though these neighboring stars span **`6.37 magnitudes` (`280x` in optical/near-IR brightness, `W1 = 11.656` to `18.028 mag`)** at unrelated distances, **all 7 inner ridge sources share virtually identical `W3` (`10.133` to `11.412 mag`, `~0.87-2.80 mJy`) and `W4` (`7.476` to `9.208 mag`, `~1.74-8.55 mJy`) flux floors**, whereas `13` off-ridge stars located `28"-59"` away are `W3`/`W4` upper limits (`snr3 = -8.0` to `+4.1`, `qph = AAUU`/`ABUU`).
   - **Extended PSF Profile Morphology (`chi2W3` Diagnostic):** Because AllWISE list-driven photometry forces PSF fitting at `W1`/`W2` stellar centroids, an underlying diffuse `12 um` PAH / `22 um` dust filament produces systematic PSF wing residuals in the `6.5"` `W3` beam. On `KIC 6634976`, the reduced PSF profile fit statistic is **`chi2W3 = 1.701`** (and `chi2W3 = 1.575` on its `15.98"` neighbor `J195515.64+420159.5`; ridge mean **`chi2W3 = 1.299`**), compared to **`chi2W3 = 0.957`** for the off-ridge upper-limit stars.

---

## 2. Radio Counterpart Check (`NVSS` 1.4 GHz & `VLASS` 3.0 GHz)

Per Requirement (2) of `bnty_0268`, `kic6634976_dossier_pipeline.py` queries both the NRAO VLA Sky Survey (`VIII/65/nvss`, Condon et al. 1998, `1.4 GHz`, `45"` FWHM, `rms ~ 0.45 mJy/beam`) and the Very Large Array Sky Survey (`J/ApJS/255/30/comp`, Gordon et al. 2021 CIRADA VLASS Epoch 1 Quick Look, `2-4 GHz`, `2.5"` FWHM, `rms ~ 0.12 mJy/beam`) via VizieR ADQL TAP around `RA = 298.819090 deg, Dec = +42.029850 deg`:

| Survey / Catalog | Frequency & Beam | Search Radius | Matches (`<= 5"`) | Matches (`<= 60"`) | Nearest Detected Source (`<= 600"`) | Nearest Source Coordinates (J2000) | Angular Separation | Integrated / Peak Flux Density |
| :--- | :--- | :---: | :---: | :---: | :--- | :--- | :---: | :---: |
| **NVSS** (`VIII/65/nvss`) | `1.4 GHz` (`45"` FWHM) | `5.0"` & `600.0"` | **`0`** | **`0`** | `NVSS 195518+415628` | `298.829000, +41.941361` | **`319.66"`** (`5.33'`) | `S_1.4 = 3.8 +/- 1.0 mJy` |
| **VLASS QL** (`J/ApJS/255/30/comp`) | `3.0 GHz` (`2.5"` FWHM) | `5.0"` & `600.0"` | **`0`** | **`0`** | `VLASS1QLCIR J195524.12+420133.6` | `298.850505, +42.026011` | **`85.14"`** (`1.42'`) | `F_tot = 4.083 +/- 0.211 mJy` (`F_peak = 4.354 +/- 0.126 mJy/beam`) |

All 5 nearest VLASS components within `600"` (`10 arcmin`):
1. `VLASS1QLCIR J195524.12+420133.6` — `sep = 85.14"`, `F_tot = 4.083 +/- 0.211 mJy`, `F_peak = 4.354 +/- 0.126 mJy/beam`
2. `VLASS1QLCIR J195503.50+415938.2` — `sep = 194.88"`, `F_tot = 1.653 +/- 0.223 mJy`, `F_peak = 1.487 +/- 0.119 mJy/beam`
3. `VLASS1QLCIR J195534.34+415808.2` — `sep = 296.45"`, `F_tot = 3.243 +/- 0.366 mJy`, `F_peak = 2.094 +/- 0.145 mJy/beam`
4. `VLASS1QLCIR J195518.69+415628.0` — `sep = 320.37"` (VLASS counterpart to `NVSS 195518+415628`), `F_tot = 2.866 +/- 0.251 mJy`, `F_peak = 2.660 +/- 0.135 mJy/beam`
5. `VLASS1QLCIR J195550.25+420229.1` — `sep = 377.41"`, `F_tot = 1.240 +/- 0.212 mJy`, `F_peak = 1.213 +/- 0.120 mJy/beam`

---

## 3. Addressing the Galactic Cirrus Finding Head-On (`Requirement 3`)

### 3.1 IRSA Dust Reddening & 100 um Cirrus Emission at `KIC 6634976`
Querying `https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr=298.819090+42.029850+equ+j2000` returns:
- **Schlegel, Finkbeiner & Davis (1998) Reddening (`E(B-V)_SFD`):**
  - Reference Pixel Value: **`0.3338 mag`**
  - `5 deg` Mean / Std / Min / Max: **`0.3294 +/- 0.0162 mag`** (`min = 0.3005 mag`, `max = 0.3574 mag`)
- **Schlafly & Finkbeiner (2011) Reddening (`E(B-V)_SandF`):**
  - Reference Pixel Value: **`0.2871 mag`** (`mean = 0.2833 +/- 0.0139 mag`)
- **IRAS/DIRBE `100 um` Surface Brightness:**
  - Reference Pixel Value: **`21.3249 MJy/sr`** (`mean = 21.0427 +/- 1.0972 MJy/sr`, `max = 22.8957 MJy/sr`)
- **Thermal Dust Temperature (`T_dust,ISM`):**
  - Reference Pixel Value: **`18.5543 K`** (`mean = 18.5470 K`)

### 3.2 Spatial Field Analysis (`II/328/allwise` within `60.0 arcsec`)
Below is the complete table of all `22` `AllWISE` sources within `60.0 arcsec` of `KIC 6634976`, sorted by angular separation:

| # | AllWISE Designation | Sep (`"`) | `W1` (mag) | `W2` (mag) | `W3` (mag) | `snr3` | `chi2W3` | `W4` (mag) | `snr4` | `chi2W4` | `qph` | Zone / Classification |
| :-: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **1** | **`J195516.58+420147.5`** | **`0.13"`** | **`12.335`** | **`12.388`** | **`10.133`** | **`23.6`** | **`1.701`** | **`7.476`** | **`10.4`** | **`1.164`** | **`AAAA`** | **Target (`KIC 6634976`, On-Ridge)** |
| **2** | **`J195515.64+420159.5`** | **`15.98"`** | **`15.408`** | **`15.906`** | **`10.376`** | **`20.1`** | **`1.575`** | **`7.960`** | **`6.4`** | **`0.959`** | **`ABAB`** | **On-Ridge Field Star (`3.07 mag` fainter in `W1`)** |
| **3** | **`J195516.36+420211.7`** | **`24.36"`** | **`15.038`** | **`15.208`** | **`10.986`** | **`13.3`** | **`1.360`** | **`8.534`** | **`3.7`** | **`1.147`** | **`AAAB`** | **On-Ridge Field Star** |
| **4** | **`J195518.94+420146.7`** | **`26.31"`** | **`11.656`** | **`11.691`** | **`10.344`** | **`20.6`** | **`1.177`** | **`7.882`** | **`7.0`** | **`1.027`** | **`AAAB`** | **On-Ridge SIMBAD `dS*` (`2MASS J19551889+4201467`)** |
| 5 | `J195516.93+420119.6` | `28.06"` | `15.297` | `15.735` | `12.595` | `1.1` | `0.954` | `9.131` | `2.3` | `0.986` | `AAUC` | Off-Ridge (`W3` Upper Limit) |
| 6 | `J195513.88+420155.8` | `31.19"` | `15.498` | `15.842` | `12.516` | `1.0` | `0.817` | `9.166` | `0.2` | `0.867` | `ABUU` | Off-Ridge (`W3/W4` Upper Limit) |
| **7** | **`J195518.29+420213.1`** | **`32.02"`** | **`14.724`** | **`14.844`** | **`11.132`** | **`11.5`** | **`1.069`** | **`9.208`** | **`2.2`** | **`0.980`** | **`AAAC`** | **On-Ridge Field Star** |
| 8 | `J195517.04+420219.8` | `32.77"` | `15.150` | `15.532` | `12.283` | `4.1` | `0.952` | `9.391` | `-0.9` | `1.084` | `AABU` | Ridge Margin (`W4` Upper Limit) |
| 9 | `J195518.99+420127.9` | `33.21"` | `15.567` | `15.738` | `12.541` | `1.5` | `1.128` | `8.602` | `1.9` | `1.105` | `ABUU` | Off-Ridge (`W3/W4` Upper Limit) |
| 10 | `J195513.57+420139.7` | `34.36"` | `13.041` | `13.159` | `12.230` | `1.7` | `0.892` | `8.819` | `1.5` | `0.972` | `AAUU` | Off-Ridge (`W3/W4` Upper Limit) |
| 11 | `J195514.50+420216.8` | `37.44"` | `15.446` | `15.547` | `13.086` | `-2.7` | `0.916` | `9.119` | `-0.3` | `0.972` | `AAUU` | Off-Ridge (`W3/W4` Upper Limit) |
| **12** | **`J195519.69+420209.4`** | **`41.09"`** | **`16.779`** | **`17.480`** | **`10.761`** | **`15.5`** | **`1.246`** | **`8.937`** | **`2.9`** | **`1.085`** | **`BUAC`** | **On-Ridge Field Star (`4.44 mag` fainter in `W1`)** |
| **13** | **`J195520.09+420133.7`** | **`41.46"`** | **`18.028`** | **`17.516`** | **`11.412`** | **`9.7`** | **`1.218`** | **`9.023`** | **`2.5`** | **`0.884`** | **`BUBC`** | **On-Ridge Field Star (`5.69 mag` fainter in `W1`)** |
| 14 | `J195517.15+420234.2` | `47.23"` | `13.916` | `14.097` | `13.134` | `-7.0` | `0.959` | `9.487` | `-3.1` | `0.970` | `AAUU` | Off-Ridge (`W3/W4` Upper Limit) |
| 15 | `J195518.31+420103.4` | `48.07"` | `12.514` | `12.670` | `13.154` | `-2.7` | `1.024` | `9.151` | `-2.9` | `0.915` | `AAUU` | Off-Ridge (`W3/W4` Upper Limit) |
| 16 | `J195514.39+420103.3` | `50.39"` | `15.335` | `15.851` | `12.821` | `-1.5` | `0.884` | `9.172` | `-0.1` | `1.002` | `ABUU` | Off-Ridge (`W3/W4` Upper Limit) |
| **17** | **`J195521.12+420159.2`** | **`51.94"`** | **`11.012`** | **`11.081`** | **`10.464`** | **`20.6`** | **`1.073`** | **`8.388`** | **`4.7`** | **`1.020`** | **`AAAB`** | **On-Ridge Bright Star** |
| 18 | `J195512.22+420126.2` | `52.93"` | `12.097` | `12.151` | `12.219` | `3.6` | `0.901` | `9.304` | `-0.8` | `0.956` | `AABU` | Off-Ridge (`W4` Upper Limit) |
| **19** | **`J195520.78+420214.8`** | **`54.26"`** | **`15.684`** | **`15.330`** | **`10.873`** | **`14.4`** | **`1.268`** | **`8.240`** | **`4.8`** | **`1.020`** | **`AAAB`** | **On-Ridge Field Star** |
| 20 | `J195521.05+420123.9` | `55.08"` | `14.183` | `14.499` | `13.183` | `-1.6` | `1.059` | `9.200` | `-0.6` | `0.932` | `AAUU` | Off-Ridge (`W3/W4` Upper Limit) |
| 21 | `J195516.12+420244.8` | `57.66"` | `15.508` | `15.520` | `13.188` | `-8.0` | `0.943` | `9.313` | `-2.8` | `0.912` | `AAUU` | Off-Ridge (`W3/W4` Upper Limit) |
| 22 | `J195511.35+420154.4` | `58.61"` | `15.268` | `16.022` | `12.912` | `-4.5` | `1.032` | `9.102` | `-1.5` | `0.931` | `ABUU` | Off-Ridge (`W3/W4` Upper Limit) |

### 3.3 Why the Data Conlusively Favors Galactic Cirrus Over a Circumstellar Technosignature
1. **Flux Floor Invariance Across Stellar Luminosity:** A circumstellar Dyson swarm or warm debris disk scales with the host star's bolometric flux (`W1`/`W2`). Instead, `J195515.64+420159.5` (`15.98"` away, `W1 = 15.408 mag`, `17x` fainter than `KIC 6634976`) has **`W3 = 10.376 mag` (`2.24 mJy`) and `W4 = 7.960 mag` (`5.48 mJy`)** — nearly identical to `KIC 6634976` (`W3 = 10.133 mag`, `2.80 mJy`; `W4 = 7.476 mag`, `8.55 mJy`). Similarly, `J195519.69+420209.4` (`41.09"` away, `W1 = 16.779 mag`, `60x` fainter) has **`W3 = 10.761 mag` (`1.57 mJy`)**.
2. **PSF Wing Broadening (`chi2W3 = 1.701`):** A true unresolved circumstellar shell at `1030 pc` (`< 1 AU` angular radius `< 0.001"`) is a pure point source in the `6.5"` `W3` beam (`chi2W3 ~ 1.0`). The measured `chi2W3 = 1.701` on `KIC 6634976` and `1.575` on `J195515.64+420159.5` prove that the `12 um` emission is spatially extended across the PSF wings.
3. **Simulated Color Temperature (`T_eff,dust ~ 189 K`):** Non-equilibrium stochastic heating of interstellar PAH molecules (`7.7, 8.6, 11.3, 12.7 um` bands in `W3`) and very small grains (`22 um` `W4`) by the interstellar radiation field produces an apparent `W3 - W4 = 2.657 mag` color that mimics a `~189 K` blackbody when fitted with a single-temperature Planck curve (`Delta_BIC = +512.46`), even though the bulk large-grain cirrus equilibrium temperature measured by IRAS/DIRBE `100 um` is `18.55 K`.

---

## 4. Kepler Q1–Q17 Long-Cadence Photometry Manifest (`STScI MAST`)

Archive URL: `https://archive.stsci.edu/missions/kepler/lightcurves/0066/006634976/`

| Quarter | FITS Filename | Raw Cadences | Valid (`SAP_QUALITY=0`) | Quarterly RMS (`ppm`) | Max Dip (`%`) | SHA-256 Prefix |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: |
| Q1 | `kplr006634976-2009166043257_llc.fits` | 1,639 | 1,432 | 162.90 | 0.0449% | `4e2500c408979d42` |
| Q2 | `kplr006634976-2009259160929_llc.fits` | 4,354 | 3,669 | 197.16 | 0.0651% | `81d402f707d8dee1` |
| Q3 | `kplr006634976-2009350155506_llc.fits` | 4,370 | 3,927 | 192.43 | 0.0626% | `e5fe18c0d1d8ef16` |
| Q4 | `kplr006634976-2010078095331_llc.fits` | 4,397 | 3,549 | 222.75 | 0.0740% | `d0a60b019e043a37` |
| Q5 | `kplr006634976-2010174085026_llc.fits` | 4,634 | 3,389 | 176.49 | 0.0672% | `7f42ca80b9f344ab` |
| Q6 | `kplr006634976-2010265121752_llc.fits` | 4,398 | 3,676 | 173.45 | 0.0566% | `98ee3f45fb4537c2` |
| Q7 | `kplr006634976-2010355172524_llc.fits` | 4,365 | 3,646 | 184.08 | 0.0609% | `47b604e05a21d082` |
| Q8 | `kplr006634976-2011073133259_llc.fits` | 3,279 | 2,535 | 189.20 | 0.0637% | `706f10145158c3ec` |
| Q9 | `kplr006634976-2011177032512_llc.fits` | 4,768 | 3,955 | 200.51 | 0.0679% | `b4958ea1490a39c6` |
| Q10 | `kplr006634976-2011271113734_llc.fits` | 4,581 | 3,841 | 187.58 | 0.0671% | `19f47b78d39f59d7` |
| Q11 | `kplr006634976-2012004120508_llc.fits` | 4,757 | 3,520 | 208.13 | 0.0725% | `10275916420e4a0c` |
| Q12 | `kplr006634976-2012088054726_llc.fits` | 4,057 | 3,130 | 207.79 | 0.0713% | `93f4b32dc7e7b31c` |
| Q13 | `kplr006634976-2012179063303_llc.fits` | 4,418 | 3,670 | 188.62 | 0.0645% | `b6434848f2e846a8` |
| Q14 | `kplr006634976-2012277125453_llc.fits` | 4,766 | 3,477 | 209.68 | 0.0766% | `a7d074e6c7b3095b` |
| Q15 | `kplr006634976-2013011073258_llc.fits` | 4,782 | 3,479 | 204.24 | 0.1252% | `f622959bc2d91746` |
| Q16 | `kplr006634976-2013098041711_llc.fits` | 4,203 | 3,234 | 190.41 | 0.0684% | `6df32d9c92a55f5d` |
| Q17 | `kplr006634976-2013131215648_llc.fits` | 1,534 | 1,071 | 196.49 | 0.0641% | `8f5e507880fa35b0` |
| **Total / Summary** | **`17` FITS Files (`1,459.49 d`)** | **`69,302`** | **`55,909`** | **Median: `192.43 ppm`** | **Max: `0.1252%`** | — |

---

## 5. Reproduction Instructions

Requires only Python 3.8+ standard library (`urllib`, `json`, `struct`, `math`, `hashlib`, `statistics`, `xml.etree.ElementTree`):

```bash
python kic6634976_dossier_pipeline.py
```

This script fetches all 17 FITS files directly from STScI MAST, executes live ADQL TAP queries against VizieR (`VIII/65/nvss`, `J/ApJS/255/30/comp`, `II/328/allwise`, `I/355/gaiadr3`) and the IRSA Dust Reddening XML service, prints the summary to stdout (`run_output.txt`), and writes the full structured dataset to `results.json`.
