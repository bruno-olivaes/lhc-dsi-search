"""Protocolo de Certificação de DSI em 3 Estágios para Espectros de Massa do LHC.

Adaptação do Protocolo Metrológico (Paper 1 / The Matrioska Protocol) para
espectros de massa invariante dimuon m_μμ.
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import lombscargle


def detrend_continuum(u, y, deg=4):
    """Ajusta e subtrai o background contínuo suave em u = ln(m/m0)."""
    coeffs = np.polyfit(u, y, deg)
    baseline = np.polyval(coeffs, u)
    residuals = y - baseline
    return baseline, residuals


def stage1_detection(u, residuals, n_omegas=5000):
    """Estágio 1: Varredura espectral Lomb-Scargle no domínio logarítmico u."""
    u_work = u - u.mean()
    y_work = residuals - residuals.mean()
    u_span = u_work.max() - u_work.min()
    omega_min = 2.0 * np.pi / u_span
    du_med = np.median(np.diff(np.sort(u_work)))
    omega_max = 2.0 * np.pi / (4.0 * du_med)  # Nyquist operacional
    omega_grid = np.logspace(np.log10(omega_min), np.log10(omega_max), n_omegas)
    power = lombscargle(u_work, y_work, omega_grid, normalize=True)
    idx = int(np.argmax(power))
    return float(omega_grid[idx]), float(power[idx]), omega_grid, power


def _sinusoid(u, amp, omega, phi):
    return amp * np.cos(omega * u + phi)


def stage2_refinement(u, residuals, omega_peak):
    """Estágio 2: Refinamento não-linear via Levenberg-Marquardt com matriz de covariância."""
    u_work = u - u.mean()
    y_work = residuals - residuals.mean()
    p0 = [np.std(y_work), omega_peak, 0.0]
    bounds = (
        [0.0, omega_peak * 0.5, -np.pi],
        [np.inf, omega_peak * 1.5, np.pi],
    )
    try:
        popt, pcov = curve_fit(_sinusoid, u_work, y_work, p0=p0, bounds=bounds, maxfev=20000)
    except Exception:
        return None

    amp, omega, phi = popt
    perr = np.sqrt(np.diag(pcov))
    sigma_omega = float(perr[1])
    ln_lambda_hat = 2.0 * np.pi / omega
    sigma_ln_lambda = ln_lambda_hat * sigma_omega / omega
    lambda_hat = float(np.exp(ln_lambda_hat))

    return {
        "amp": float(amp),
        "omega": float(omega),
        "phi": float(phi),
        "sigma_omega": sigma_omega,
        "ln_lambda_hat": float(ln_lambda_hat),
        "lambda_hat": lambda_hat,
        "sigma_lambda": float(lambda_hat * sigma_ln_lambda),
        "rel_err_lambda": float(sigma_ln_lambda),
    }


def _scramble_phases(u, y, rng):
    u_grid = np.linspace(u.min(), u.max(), len(u))
    y_uniform = np.interp(u_grid, u, y)
    F = np.fft.rfft(y_uniform)
    phases = rng.uniform(0.0, 2.0 * np.pi, len(F))
    F_scram = np.abs(F) * np.exp(1j * phases)
    y_scram_uniform = np.fft.irfft(F_scram, n=len(u_grid))
    return np.interp(u, u_grid, y_scram_uniform)


def stage3_certify(u, residuals, stage2, sigma=None, n_bootstrap=1000, rng=None):
    """Estágio 3: Auditoria dos 5 Gates Mandatórios de Integridade Metrológica."""
    if rng is None:
        rng = np.random.default_rng(42)

    u_work = u - u.mean()
    u_span = u_work.max() - u_work.min()
    du_med = np.median(np.diff(np.sort(u_work)))
    omega_min = 2.0 * np.pi / u_span
    omega_max = 2.0 * np.pi / (4.0 * du_med)
    omega = stage2["omega"]

    # Gate 1: Resolubilidade espectral (Rayleigh & Nyquist)
    gate_resolvable = bool(omega_min <= omega <= omega_max)

    # Gate 2: Significância de Amplitude (> 2.0 * RMS)
    y_work = residuals - residuals.mean()
    rms_noise = float(np.sqrt(np.mean(y_work**2)))
    gate_amplitude = bool(abs(stage2["amp"]) > 2.0 * rms_noise)

    # Gate 3: Precisão Paramétrica (sigma_lambda / lambda_hat < 0.20)
    gate_precision = bool(stage2["rel_err_lambda"] < 0.20)

    # Gate 4: Look-Elsewhere Effect (p_corrected < 0.05 via Bootstrap de Máximo)
    max_powers = np.zeros(n_bootstrap)
    for i in range(n_bootstrap):
        y_scram = _scramble_phases(u_work, y_work, rng)
        omega_grid = np.linspace(omega * 0.5, omega * 1.5, 500)
        max_powers[i] = float(
            np.max(
                lombscargle(
                    u_work, y_scram - y_scram.mean(), omega_grid, normalize=True
                )
            )
        )
    power_obs = float(
        np.atleast_1d(lombscargle(u_work, y_work, np.array([omega]), normalize=True))[0]
    )
    p_corrected = float(np.mean(max_powers >= power_obs))
    gate_plevel = bool(p_corrected < 0.05)

    # Gate 5: Critério de Informação Bayesiana (Delta BIC = BIC_null - BIC_DSI > 10.0)
    err = np.full_like(y_work, rms_noise) if sigma is None else np.asarray(sigma)
    err = np.where(err <= 0, rms_noise, err)
    y_model = _sinusoid(u_work, stage2["amp"], omega, stage2["phi"])

    # Log-likelihood gaussiana
    loglike_null = -0.5 * np.sum((y_work / err)**2 + np.log(2.0 * np.pi * err**2))
    loglike_dsi = -0.5 * np.sum(((y_work - y_model) / err)**2 + np.log(2.0 * np.pi * err**2))
    n_pts = len(y_work)
    bic_null = -2.0 * loglike_null
    bic_dsi = 3.0 * np.log(n_pts) - 2.0 * loglike_dsi
    delta_bic = float(bic_null - bic_dsi)
    gate_bic = bool(delta_bic > 10.0)

    certified = all(
        [
            gate_resolvable,
            gate_amplitude,
            gate_precision,
            gate_plevel,
            gate_bic,
        ]
    )

    return {
        "gate_resolvable": gate_resolvable,
        "gate_amplitude": gate_amplitude,
        "gate_precision": gate_precision,
        "gate_plevel": gate_plevel,
        "gate_bic": gate_bic,
        "power_observed": power_obs,
        "p_corrected": p_corrected,
        "delta_bic": delta_bic,
        "certified": certified,
        "rms_noise": rms_noise,
    }
