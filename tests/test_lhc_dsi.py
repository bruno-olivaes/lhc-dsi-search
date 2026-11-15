"""Unit and physics validation tests for lhc-dsi-search (Paper 6)."""
import os
import sys

import numpy as np

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)

from src.protocol_lhc import (
    detrend_continuum,
    stage1_detection,
    stage2_refinement,
    stage3_certify,
)


def test_white_noise_null_rejected():
    """Pure Gaussian white noise must fail the 5-gate DSI certification protocol."""
    rng = np.random.default_rng(2026)
    u = np.linspace(2.4, 4.3, 120)
    y = rng.normal(0.0, 1.0, size=len(u))
    omega_peak, _, _, _ = stage1_detection(u, y)
    s2 = stage2_refinement(u, y, omega_peak)
    assert s2 is not None
    cert = stage3_certify(u, y, s2, sigma=np.ones_like(y), n_bootstrap=100, rng=rng)
    assert cert["certified"] is False
    assert cert["delta_bic"] < 10.0


def test_synthetic_dsi_signal_certified():
    """Strong log-periodic harmonic injection must be recovered and pass all 5 gates."""
    rng = np.random.default_rng(77)
    u = np.linspace(2.4, 4.3, 180)
    lam_true = 1.65
    omega_true = 2.0 * np.pi / np.log(lam_true)
    smooth_bg = 500.0 * np.exp(-1.5 * (u - 2.4))
    y = smooth_bg + 8.0 * np.cos(omega_true * (u - u.mean()) + 0.4) + rng.normal(0.0, 1.0, size=len(u))

    _, residuals = detrend_continuum(u, y, deg=5)
    omega_peak, _, _, _ = stage1_detection(u, residuals)
    s2 = stage2_refinement(u, residuals, omega_peak)
    assert s2 is not None
    assert abs(s2["lambda_hat"] - lam_true) / lam_true < 0.05

    cert = stage3_certify(u, residuals, s2, sigma=np.ones_like(residuals), n_bootstrap=150, rng=rng)
    assert cert["gate_resolvable"] is True
    assert cert["gate_amplitude"] is True
    assert cert["gate_precision"] is True
    assert cert["gate_plevel"] is True
    assert cert["gate_bic"] is True
    assert cert["certified"] is True

