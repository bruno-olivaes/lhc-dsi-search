# lhc-dsi-search: A Certified Search for Discrete Scale Invariance in LHC Dimuon Mass Spectra

[![Paper PDF](https://img.shields.io/badge/Paper-PDF_Manuscript-red.svg)](manuscript/main.pdf)
[![Data: CERN Open Data](https://img.shields.io/badge/Data-CERN%20Open%20Data%20CMS-blue.svg)](https://doi.org/10.7483/OPENDATA.CMS.7001.303)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python: 3.12](https://img.shields.io/badge/Python-3.12+-green.svg)](https://www.python.org/)
[![CDD: 539.7 & 006.3](https://img.shields.io/badge/CDD-539.7%20%7C%20006.3-purple.svg)](https://github.com/bruno-olivaes/dsi-cosmo-search/blob/main/CDD_PORTFOLIO_MAP.md)

Open-source scientific pipeline adapting the pre-registered three-stage Discrete Scale Invariance (DSI) metrological protocol (from Paper 1, [`dsi-cosmo-search`](https://github.com/bruno-olivaes/dsi-cosmo-search)) to invariant mass spectra of elementary particles recorded at the Large Hadron Collider (LHC).

---

## 🔬 Summary of Scientific Findings

We analyze **62,214 physical dimuon ($\mu^+\mu^-$) collision events** from the CMS experiment at CERN ($\sqrt{s} = 7\,\text{TeV}$) retrieved from the CERN Open Data Portal, transforming invariant mass coordinates to logarithmic space $u = \ln(m / 1\,\text{GeV})$.

| Kinematic Window | Candidate Mode | Preferred Ratio $\hat\lambda$ | Significance ($p_{\text{corrected}}$) | Failed Gates | Final Protocol Verdict |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Drell–Yan Continuum ($11 \le m \le 75\,\text{GeV}$)** | $\hat\omega = 5.898\,\text{rad}$ | $\hat\lambda = 2.902 \pm 0.075$ | $p = 0.324$ | Gate 2 (Amp), Gate 4 (Bootstrap) | **NULL RESULT (REJECTED)** |
| **Full Spectrum ($2 \le m \le 110\,\text{GeV}$)** | $\hat\omega = 3.804\,\text{rad}$ | $\hat\lambda = 5.217 \pm 0.152$ | $p = 0.039$ | Gate 2 (Amplitude < $2\sigma$ RMS) | **NULL RESULT (REJECTED)** |

### Key Conclusions:
1. **Zero False-Positive Hallucination:** The 5-gate protocol successfully rejects spurious log-periodic candidate modulations in both kinematic regimes.
2. **Absence of Cross-Domain Universality:** The LHC continuum candidate ($\hat\lambda_{\text{LHC}} = 2.902 \pm 0.075$) is inconsistent with cosmic microwave background bounds from Planck 2018 ($\hat\lambda_{\text{Planck}} = 4.1208 \pm 0.3092$, $p=0.975$). Both represent independent realizations of Poissonian background noise.
3. **Modulo-9 Exploratory Analysis:** While standard vector meson masses ($J/\psi$, $\Upsilon$, $Z^0$) exhibit phenomenological power-of-3 ratio clustering ($m_\Upsilon / m_{J/\psi} \approx 3^1$, $m_Z / m_\Upsilon \approx 3^2$, $m_Z / m_{J/\psi} \approx 3^3$), this does not correspond to a continuous log-periodic DSI modulation in LHC collision residuals.

---

## 🛠️ Quickstart & Reproduction

```bash
# Clone the repository
git clone https://github.com/bruno-olivaes/lhc-dsi-search.git
cd lhc-dsi-search

# Install dependencies
pip install numpy scipy matplotlib

# Run full pipeline and produce publication figure
python scripts/run_analysis.py

# Compile manuscript
cd manuscript
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

---

## 📜 Citation & License

```bibtex
@article{olivaes2026lhcdsi,
  title={A Certified Search for Discrete Scale Invariance in LHC Dimuon Mass Spectra},
  author={Olivaes, Bruno Lu\c{c}ardo},
  journal={Preprint},
  year={2026},
  url={https://github.com/bruno-olivaes/lhc-dsi-search}
}
```

Licensed under the [MIT License](LICENSE).
