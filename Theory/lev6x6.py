"""Core model of LevitationSim6x6.ipynb as an importable module.

Extracted verbatim (definitions only, no tests or example runs) from the code cells of
Theory/LevitationSim6x6.ipynb so that other studies can reuse the simulator without
modifying that notebook. See the notebook for the physics derivations.
"""
from dataclasses import dataclass, field, replace

from pathlib import Path

from typing import Callable, Optional

import numpy as np

import matplotlib.pyplot as plt

from matplotlib import animation, colors

from scipy.integrate import solve_ivp, quad

from scipy.special import elliprf, elliprd

from IPython.display import display, HTML, SVG, Markdown

plt.rcParams.update({"figure.dpi": 90, "axes.grid": True, "grid.alpha": 0.3, "font.size": 9,
                     "axes.formatter.useoffset": False})

NB_DIR = Path(__file__).resolve().parent

FIG_DIR = NB_DIR / "figures"

TEST_RESULTS = {}

def _report(label, ok, got, expected, err, kind):
    TEST_RESULTS[label] = bool(ok)
    tag = "PASS" if ok else "FAIL"
    print(f"{tag}  {label}\n      got {np.array2string(np.asarray(got), precision=8)}   "
          f"expected {np.array2string(np.asarray(expected), precision=8)}   {kind} err {err:.2e}")

def check_close(label, got, expected, rtol=1e-3):
    """Relative check on the 2-norm (works for scalars, vectors, complex numbers)."""
    g, e = np.atleast_1d(np.asarray(got, dtype=complex)), np.atleast_1d(np.asarray(expected, dtype=complex))
    err = np.linalg.norm(g - e) / max(np.linalg.norm(e), 1e-300)
    _report(label, err <= rtol, got, expected, err, "rel.")

def check_abs(label, got, expected, atol):
    err = np.linalg.norm(np.atleast_1d(np.asarray(got, dtype=complex) - np.asarray(expected, dtype=complex)))
    _report(label, err <= atol, got, expected, err, "abs.")

def quat_to_R(q):
    q0, q1, q2, q3 = q
    return np.array([[1 - 2*(q2*q2 + q3*q3), 2*(q1*q2 - q0*q3),     2*(q1*q3 + q0*q2)],
                     [2*(q1*q2 + q0*q3),     1 - 2*(q1*q1 + q3*q3), 2*(q2*q3 - q0*q1)],
                     [2*(q1*q3 - q0*q2),     2*(q2*q3 + q0*q1),     1 - 2*(q1*q1 + q2*q2)]])

def q_dot(q, w):
    q0, q1, q2, q3 = q
    return 0.5*np.array([-q1*w[0] - q2*w[1] - q3*w[2],
                          q0*w[0] - q3*w[1] + q2*w[2],
                          q3*w[0] + q0*w[1] - q1*w[2],
                         -q2*w[0] + q1*w[1] + q0*w[2]])

def q_mul(a, b):
    a0, av, b0, bv = a[0], np.asarray(a[1:]), b[0], np.asarray(b[1:])
    return np.concatenate([[a0*b0 - av @ bv], a0*bv + b0*av + np.cross(av, bv)])

def normalize_q(q):
    q = np.asarray(q)
    return q / np.sqrt(q @ q)

def axis_angle_q(axis, deg):
    k = np.asarray(axis, float) / np.linalg.norm(axis)
    h = np.deg2rad(deg) / 2
    return np.concatenate([[np.cos(h)], np.sin(h)*k])

def skew(a):
    return np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])

def tilt_deg(q):
    return np.degrees(np.arccos(np.clip(quat_to_R(normalize_q(q))[2, 2], -1, 1)))

I3 = np.eye(3)

Z3 = np.zeros((3, 3))

def iso3(s):           return s * np.eye(3)

def diag3(a, b, c):    return np.diag([a, b, c]).astype(float)

def from_packed(t):    # the JS storage format [t11, t12, t13, t22, t23, t33]
    t11, t12, t13, t22, t23, t33 = t
    return np.array([[t11, t12, t13], [t12, t22, t23], [t13, t23, t33]], float)

def symmetrize(m):     return (m + m.T) / 2

def resistance_from_blocks(K, B, Om):  return np.block([[K, B], [B.T, Om]]).astype(float)

def thermo_from_blocks(Gam, Lam):      return np.vstack([Gam, Lam]).astype(float)

def shift_matrix(r):                   return np.block([[I3, -skew(r)], [Z3, I3]])

@dataclass
class Particle:
    name: str
    mass: float
    inertia: np.ndarray                      # 3x3 about COM, body frame
    res: np.ndarray                          # 6x6 resistance about H
    thermo: np.ndarray                       # 6x3 thermophoretic matrix about H
    rH: np.ndarray = field(default_factory=lambda: np.zeros(3))   # COM -> H, body frame
    shape: dict = field(default_factory=lambda: {"type": "cylinder", "half_length": 0.10, "radius": 0.034})
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.inertia = np.asarray(self.inertia, float)
        self.res = np.asarray(self.res, float)
        self.thermo = np.asarray(self.thermo, float)
        self.rH = np.asarray(self.rH, float)

    def at_com(self):
        """Return (R_COM, T_COM): the tensors referred to the centre of mass."""
        P = shift_matrix(self.rH)
        return P.T @ self.res @ P, P.T @ self.thermo

    @property
    def Gamma(self):  return self.thermo[:3]
    @property
    def length_scale(self):
        if "length_scale" in self.meta: return self.meta["length_scale"]
        s = self.shape
        return max(s["semiaxes"]) if s["type"] == "ellipsoid" else s.get("half_length", 1.0)

def validate_particle(p, verbose=True):
    r6, i3 = p.res, p.inertia
    asym = np.linalg.norm(r6 - r6.T) / max(np.linalg.norm(r6), 1e-300)
    d = np.sqrt(np.abs(np.diag(r6))); d[d == 0] = 1.0
    ev_r = np.linalg.eigvalsh(symmetrize(r6) / np.outer(d, d))       # scale-free PD test
    ev_i = np.sort(np.linalg.eigvalsh(symmetrize(i3)))[::-1]
    rows = [("Res is 6x6 and Thermo is 6x3", r6.shape == (6, 6) and p.thermo.shape == (6, 3), ""),
            ("Res symmetric (reciprocal theorem): |R-R^T|/|R|", asym < 1e-8, f"{asym:.1e}"),
            ("Res positive definite (second law): min eig of scaled R", ev_r.min() > 0, f"{ev_r.min():.3g}"),
            ("Inertia symmetric", np.allclose(i3, i3.T, rtol=1e-12, atol=0), ""),
            ("Inertia positive definite", ev_i.min() > 0, np.array2string(ev_i, precision=4)),
            ("Principal moments obey I1 <= I2 + I3", ev_i[0] <= (ev_i[1] + ev_i[2])*(1 + 1e-9), ""),
            ("Mass > 0", p.mass > 0, f"{p.mass:.4g}")]
    if verbose:
        md = "| check | status | value |\n|---|---|---|\n" + "\n".join(
            f"| {c} | {'OK' if ok else '**CHECK**'} | {v} |" for c, ok, v in rows)
        display(Markdown(f"**{p.name}**\n\n" + md))
    return all(ok for _, ok, _ in rows)

def show_matrix(m, title=""):
    body = " \\\\ ".join(" & ".join("0" if x == 0 else f"{x:.4g}".replace("e", r"\mathrm{e}") for x in row) for row in m)
    display(Markdown(f"{title} $\\begin{{pmatrix}}{body}\\end{{pmatrix}}$"))

LEGACY_JS_CONFIG = dict(name="sim_symmetric.jsx defaults", M=1.0, g=9.81, Idiag=[0.02, 0.02, 0.008],
                        Xi=[2, 0, 0, 2, 0, 2], Gamma=[5, 0, 0, 5, 0, 5], Cr=[0.5, 0, 0, 0.5, 0, 0.5],
                        Lambda=[0, 0, 0, 0, 0, 0], ar=0.5, az=0.5,
                        x0=[0.5, 0.3, 0.2], v0=[0, 0, 0], w0=[0, 0.5, 2], tilt=30, tEnd=20, nSteps=1200)

def from_legacy_js(cfg):
    return Particle(name=cfg.get("name", "legacy JS configuration"), mass=cfg["M"],
                    inertia=np.diag(cfg["Idiag"]),
                    res=resistance_from_blocks(from_packed(cfg["Xi"]), Z3, from_packed(cfg["Cr"])),
                    thermo=thermo_from_blocks(from_packed(cfg["Gamma"]), from_packed(cfg["Lambda"])))

@dataclass
class Environment:
    g: float
    A: np.ndarray
    G0: np.ndarray
    g0_mode: str
    q_ref: np.ndarray
    field: Callable

def balanced_G0(p, g, q_ref=(1, 0, 0, 0)):
    R0 = quat_to_R(normalize_q(np.asarray(q_ref, float)))
    return np.array([0.0, 0.0, -p.mass * g * np.linalg.inv(R0 @ p.Gamma @ R0.T)[2, 2]])

def make_environment(p, alpha=(0.5, 0.5), g0_mode="Balanced", G0=None, q_ref=(1, 0, 0, 0),
                     custom_field=None, g=9.81):
    A = np.diag([alpha[0], alpha[0], alpha[1]]).astype(float) if np.ndim(alpha) == 1 else symmetrize(np.asarray(alpha, float))
    if g0_mode == "Balanced":   G0v = balanced_G0(p, g, q_ref)
    elif g0_mode == "LegacyJS": G0v = np.array([0.0, 0.0, -p.mass * g / p.Gamma[1, 1]])
    elif g0_mode == "Fixed":    G0v = np.asarray(G0, float)
    else: raise ValueError("g0_mode must be 'Balanced', 'Fixed' or 'LegacyJS'")
    fld = custom_field if custom_field is not None else (lambda x, A=A, G0v=G0v: A @ x + G0v)
    return Environment(g, A, G0v, g0_mode, np.asarray(q_ref, float), fld)

def alpha_from_trap_frequencies(p, fr, fz):
    G = p.Gamma
    return (p.mass*(2*np.pi*fr)**2 / np.mean([G[0, 0], G[1, 1]]), p.mass*(2*np.pi*fz)**2 / G[2, 2])

def required_temperature_gradient(env, T0):
    """|grad T| (K/m) at the trap centre."""
    return np.linalg.norm(env.G0) * T0

def conduction_consistency(env, x=(0, 0, 0), n=0.8):
    return {"div G (model)": np.trace(env.A),
            "div G required by source-free conduction": -(n + 1)*np.linalg.norm(env.field(np.asarray(x, float)))**2}

def make_rhs(p, env):
    Rc, Tc = p.at_com()
    M, I, Iinv, rH, g, fld = p.mass, p.inertia, np.linalg.inv(p.inertia), p.rH, env.g, env.field
    ez = np.array([0.0, 0.0, 1.0])

    def rhs(t, s):
        x, v, w = s[0:3], s[3:6], s[10:13]
        q = normalize_q(s[6:10])
        R = quat_to_R(q)
        Gb = R.T @ fld(x + R @ rH)
        U = np.concatenate([R.T @ v, w])
        FT, FD = -Tc @ Gb, -Rc @ U
        F = FT + FD
        return np.concatenate([v, R @ F[:3] / M - g*ez, q_dot(q, w),
                               Iinv @ (F[3:] - np.cross(w, I @ w)), [U @ FT, -(U @ FD)]])
    return rhs

@dataclass
class Trajectory:
    t: np.ndarray
    s: np.ndarray                     # (n, 15): 13 state variables + W_T, W_D
    particle: Particle
    env: Environment
    method: str
    dt: Optional[float] = None
    sol: Optional[Callable] = None

def initial_state(x0, v0, w0, tilt=0.0, axis=(1, 0, 0), q0=None):
    q = normalize_q(np.asarray(q0, float)) if q0 is not None else axis_angle_q(axis, tilt)
    return np.concatenate([np.asarray(x0, float), np.asarray(v0, float), q, np.asarray(w0, float)])

def rk4_simulate(p, env, s0, t_end, n_steps):
    rhs, dt = make_rhs(p, env), t_end / n_steps
    s = np.concatenate([s0, [0.0, 0.0]]); out = np.empty((n_steps + 1, 15)); out[0] = s
    for i in range(n_steps):
        k1 = rhs(0, s); k2 = rhs(0, s + 0.5*dt*k1); k3 = rhs(0, s + 0.5*dt*k2); k4 = rhs(0, s + dt*k3)
        s = s + dt*(k1 + 2*k2 + 2*k3 + k4)/6
        s[6:10] = normalize_q(s[6:10])
        out[i + 1] = s
    return Trajectory(dt*np.arange(n_steps + 1), out, p, env, "RK4", dt=dt)

def nd_simulate(p, env, s0, t_end, frames=601, rtol=1e-10, atol=1e-10, method="LSODA"):
    rhs, L, g = make_rhs(p, env), p.length_scale, env.g
    sc = np.concatenate([[L]*3, [np.sqrt(g*L)]*3, [1.0]*4, [np.sqrt(g/L)]*3, [p.mass*g*L]*2])
    sol = solve_ivp(lambda t, y: rhs(t, sc*y)/sc, (0, t_end), np.concatenate([s0, [0, 0]])/sc,
                    method=method, rtol=rtol, atol=atol, dense_output=True)
    if not sol.success: raise RuntimeError(sol.message)
    dense = lambda t: sol.sol(t) * (sc if np.ndim(t) == 0 else sc[:, None])
    ts = np.linspace(0, t_end, frames)
    S = dense(ts).T.copy()
    S[:, 6:10] /= np.linalg.norm(S[:, 6:10], axis=1, keepdims=True)
    return Trajectory(ts, S, p, env, f"solve_ivp/{method} ({sol.nfev} rhs calls)", sol=dense)

def trajectory_diagnostics(traj):
    p, env = traj.particle, traj.env
    Rc, Tc = p.at_com()
    d = {k: [] for k in ["x", "v", "w", "speed", "tilt", "qnorm_err", "Ekin", "Egrav", "FT_lab", "tauT", "axis", "xH", "PT", "PD"]}
    for s in traj.s:
        x, v, qn, w = s[0:3], s[3:6], s[6:10], s[10:13]
        R = quat_to_R(normalize_q(qn))
        Gb = R.T @ env.field(x + R @ p.rH)
        U = np.concatenate([R.T @ v, w]); FT = -Tc @ Gb
        for k, val in [("x", x), ("v", v), ("w", w), ("speed", np.linalg.norm(v)),
                       ("tilt", np.degrees(np.arccos(np.clip(R[2, 2], -1, 1)))), ("qnorm_err", np.sqrt(qn @ qn) - 1),
                       ("Ekin", 0.5*p.mass*v @ v + 0.5*w @ p.inertia @ w), ("Egrav", p.mass*env.g*x[2]),
                       ("FT_lab", R @ FT[:3]), ("tauT", FT[3:]), ("axis", R[:, 2]), ("xH", x + R @ p.rH),
                       ("PT", U @ FT), ("PD", U @ Rc @ U)]:
            d[k].append(val)
    d = {k: np.array(v) for k, v in d.items()}
    d["t"], d["WT"], d["WD"] = traj.t, traj.s[:, 13], traj.s[:, 14]
    d["dEmech"] = d["Ekin"] + d["Egrav"] - d["Ekin"][0] - d["Egrav"][0]
    d["E_residual"] = d["dEmech"] - (d["WT"] - d["WD"])
    scale = max(np.abs(d["dEmech"]).max(), np.abs(d["WT"]).max(), np.abs(d["WD"]).max(), 1e-300)
    d["E_residual_rel"] = d["E_residual"] / scale
    return d

def diagnostic_plots(traj, title=None):
    d, t = trajectory_diagnostics(traj), traj.t
    fig, axs = plt.subplots(4, 3, figsize=(13, 12))
    def ts(ax, ys, labels, ttl, ylab):
        for y, lab in zip(ys, labels): ax.plot(t, y, label=lab)
        ax.set_title(ttl); ax.set_xlabel("t (s)"); ax.set_ylabel(ylab)
        if any(labels): ax.legend(fontsize=7)
    ts(axs[0, 0], d["x"].T, "xyz", "COM position (lab)", "m")
    ts(axs[0, 1], [d["speed"]], [None], "Speed |v|", "m/s")
    ts(axs[0, 2], d["w"].T, ["ω1", "ω2", "ω3"], "Angular velocity (body frame)", "rad/s")
    ts(axs[1, 0], [d["tilt"]], [None], "Tilt of body axis e3 from vertical", "deg")
    ts(axs[1, 1], [d["dEmech"], d["WT"], -d["WD"]], ["Δ(Ekin+Egrav)", "work by thermophoresis", "−dissipated"], "Energy budget", "J")
    ts(axs[1, 2], [d["E_residual_rel"]], [None], "Energy-balance residual / energy scale", "")
    ts(axs[2, 0], [d["qnorm_err"]], [None], "Quaternion norm error |q|−1", "")
    ts(axs[2, 1], d["FT_lab"].T, ["Fx", "Fy", "Fz"], "Thermophoretic force (lab)", "N")
    ts(axs[2, 2], d["tauT"].T, ["τ1", "τ2", "τ3"], "Thermophoretic torque about COM (body)", "N m")
    axs[3, 0].plot(d["x"][:, 2], d["v"][:, 2]); axs[3, 0].set(title="Axial phase portrait", xlabel="z (m)", ylabel="vz (m/s)")
    axs[3, 1].plot(d["x"][:, 0], d["x"][:, 1]); axs[3, 1].set(title="Radial trajectory (top view)", xlabel="x (m)", ylabel="y (m)")
    axs[3, 1].set_aspect("equal", "datalim")
    ax = axs[3, 2]; floor = 1e-12*max(np.abs(d["PD"]).max(), np.abs(d["PT"]).max(), 1e-300)
    ax.semilogy(t, np.maximum(d["PD"], floor), label="dissipated power")
    ax.semilogy(t, np.maximum(np.abs(d["PT"]), floor), label="|thermophoretic power|")
    ax.set(title="Power", xlabel="t (s)", ylabel="W"); ax.legend(fontsize=7)
    fig.suptitle(title or f"{traj.particle.name}  —  {traj.method}", y=1.0)
    fig.tight_layout(); plt.show()

def linearize(p, env, x_star, q_star, h=1e-30):
    rhs, qs, xs = make_rhs(p, env), normalize_q(np.asarray(q_star, float)), np.asarray(x_star, float)
    def f12(y):
        q = q_mul(qs, np.concatenate([[1.0], y[6:9]/2]))
        d = rhs(0, np.concatenate([xs + y[0:3], y[3:6], q, y[9:12], [0, 0]]))
        return np.concatenate([d[0:3], d[3:6], y[9:12], d[10:13]])
    J = np.empty((12, 12))
    for j in range(12):
        y = np.zeros(12, complex); y[j] = 1j*h
        J[:, j] = f12(y).imag / h
    return J

def sorted_eigs(J):
    ev = np.linalg.eigvals(J)
    ev = np.where(np.abs(ev) < 1e-9*max(1, np.abs(ev).max()), 0, ev)
    return ev[np.lexsort((np.abs(ev.imag), -ev.real))]

def mode_table(J, title="Linear modes"):
    rows = []
    for lam in sorted_eigs(J):
        kind = "stable" if lam.real < 0 else ("neutral" if lam.real == 0 else "**UNSTABLE**")
        Q = f"{abs(lam)/(-2*lam.real):.3g}" if (lam.real < 0 and lam.imag != 0) else "–"
        rows.append(f"| {lam.real:.6g} {'+' if lam.imag >= 0 else '−'} {abs(lam.imag):.6g} i | {-lam.real:.4g} | {abs(lam.imag)/(2*np.pi):.4g} | {Q} | {kind} |")
    display(Markdown(f"**{title}**\n\n| eigenvalue (1/s) | decay rate (1/s) | frequency (Hz) | Q | type |\n|---|---|---|---|---|\n" + "\n".join(rows)))

def rk4_stability(J, dt):
    z = np.linalg.eigvals(J)*dt
    amp = np.abs(1 + z + z**2/2 + z**3/6 + z**4/24).max()
    return {"dt": dt, "max|λ|dt": np.abs(z).max(), "max RK4 amplification": amp,
            "verdict": "stable" if amp <= 1 + 1e-12 else "UNSTABLE: reduce dt or use nd_simulate"}

def equilibrium_prediction(p, env, q=(1, 0, 0, 0)):
    qn = normalize_q(np.asarray(q, float)); R = quat_to_R(qn)
    w = -p.mass*env.g*np.linalg.inv(R @ p.Gamma @ R.T) @ np.array([0, 0, 1.0])
    xH = np.linalg.pinv(env.A) @ (w - env.G0)
    xs = xH - R @ p.rH
    acc = make_rhs(p, env)(0, np.concatenate([xs, [0, 0, 0], qn, [0, 0, 0], [0, 0]]))
    return {"xCOM": xs, "xH": xH, "q": qn, "linear_acc_residual": acc[3:6], "angular_acc_residual": acc[10:13]}

def body_mesh(p, n_th=18, n_ph=28):
    """Return (X, Y, Z, facecolors) in body coordinates."""
    s = p.shape
    if s["type"] == "ellipsoid":
        a1, a2, a3 = s["semiaxes"]; eps = s.get("epsilon", 0.0)
        th, ph = np.meshgrid(np.linspace(0, np.pi, n_th), np.linspace(0, 2*np.pi, n_ph))
        X, Y, Z = a1*np.sin(th)*np.cos(ph), a2*np.sin(th)*np.sin(ph), a3*np.cos(th)
        shade = 0.3*np.ones_like(Z) if eps == 0 else np.clip((eps*Z + abs(eps)*a3)/(2*abs(eps)*a3), 0, 1)
        return X, Y, Z, plt.cm.Blues(0.25 + 0.7*shade)
    hl, r = s["half_length"], s["radius"]
    zz, ph = np.meshgrid(np.linspace(-hl, hl, 8), np.linspace(0, 2*np.pi, n_ph))
    X, Y, Z = r*np.cos(ph), r*np.sin(ph), zz
    fc = np.where((Z > 0.8*hl)[..., None], np.array(colors.to_rgba("tomato")), np.array(colors.to_rgba("#5a9fe6")))
    return X, Y, Z, fc

def scene_data(traj):
    xs = traj.s[:, 0:3]; L = traj.particle.length_scale
    lo, hi = xs.min(0), xs.max(0); ctr = (lo + hi)/2; half = max((hi - lo).max()/2 + 1.5*L, 2*L)
    return {"xs": xs, "Rs": np.array([quat_to_R(normalize_q(s[6:10])) for s in traj.s]),
            "mesh": body_mesh(traj.particle), "rH": traj.particle.rH, "L": L, "t": traj.t,
            "lims": np.stack([ctr - half, ctr + half], 1)}

def draw_frame(ax, sd, i, elev=25, azim=38):
    ax.cla()
    x, R, L = sd["xs"][i], sd["Rs"][i], sd["L"]
    X, Y, Z, fc = sd["mesh"]
    P = np.einsum("ij,jkl->ikl", R, np.stack([X, Y, Z])) + x[:, None, None]
    ax.plot_surface(P[0], P[1], P[2], facecolors=fc, linewidth=0, antialiased=False, shade=False, alpha=0.85)
    ax.plot(*sd["xs"][:i + 1].T, color="#4fc3f7", lw=1.2)
    for k, c in enumerate(["darkred", "darkgreen", "blue"]):
        ax.quiver(*x, *(1.8*L*R[:, k]), color=c, lw=1.5, arrow_length_ratio=0.15)
    ax.scatter(*x, color="k", s=12); ax.scatter(*(x + R @ sd["rH"]), color="red", s=12)
    ax.plot([x[0], x[0]], [x[1], x[1]], [x[2], sd["lims"][2, 0]], color="gray", ls="--", lw=0.8)
    ax.set_xlim(*sd["lims"][0]); ax.set_ylim(*sd["lims"][1]); ax.set_zlim(*sd["lims"][2])
    ax.set_box_aspect((1, 1, 1)); ax.view_init(elev, azim)
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_zlabel("z (m)")
    ax.set_title(f"t = {sd['t'][i]:.4g} s", fontsize=9)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis): a.set_tick_params(labelsize=6)

def _viewer_figure(traj):
    d = trajectory_diagnostics(traj)
    fig = plt.figure(figsize=(11, 5.6))
    ax3 = fig.add_axes([0.0, 0.03, 0.52, 0.92], projection="3d")
    panels = [("Position (lab)", d["x"].T, "m"), ("Body angular velocity", d["w"].T, "rad/s"),
              ("Tilt", [d["tilt"]], "deg"), ("Speed", [d["speed"]], "m/s")]
    cursors = []
    for k, (ttl, ys, yl) in enumerate(panels):
        ax = fig.add_axes([0.62, 0.78 - 0.235*k, 0.36, 0.17])
        for y in ys: ax.plot(traj.t, y, lw=1)
        ax.set_title(ttl, fontsize=8, pad=2); ax.set_ylabel(yl, fontsize=7); ax.tick_params(labelsize=6)
        cursors.append(ax.axvline(0, color="red", lw=1))
    return fig, ax3, cursors

def animate_trajectory(traj, n_frames=60, elev=25, azim=38, fps=15):
    sd = scene_data(traj)
    idx = np.unique(np.linspace(0, len(traj.t) - 1, n_frames).round().astype(int))
    fig, ax3, cursors = _viewer_figure(traj)
    def update(k):
        i = idx[k]; draw_frame(ax3, sd, i, elev, azim)
        for c in cursors: c.set_xdata([traj.t[i]] * 2)
    anim = animation.FuncAnimation(fig, update, frames=len(idx), interval=1000/fps)
    plt.close(fig)
    return anim

def show_animation(traj, dpi=55, **kw):
    """Embed an HTML/JS player (JPEG frames keep the notebook small)."""
    with plt.rc_context({"animation.frame_format": "jpeg", "figure.dpi": dpi, "savefig.dpi": dpi}):
        display(HTML(animate_trajectory(traj, **kw).to_jshtml(default_mode="once")))

def save_animation(traj, path, n_frames=120, fps=20, **kw):
    animate_trajectory(traj, n_frames=n_frames, fps=fps, **kw).save(path, writer=animation.PillowWriter(fps=fps))

def interactive_viewer(traj):
    try:
        import ipywidgets as W
    except ImportError:
        print("ipywidgets is not installed (pip install ipywidgets); use show_animation(traj) instead.")
        return
    sd = scene_data(traj); n = len(traj.t)
    play = W.Play(min=0, max=n - 1, step=max(1, n//200), interval=40)
    frame = W.IntSlider(min=0, max=n - 1, description="frame"); W.jslink((play, "value"), (frame, "value"))
    elev = W.IntSlider(25, -80, 80, description="Elev"); azim = W.IntSlider(38, -180, 180, description="Azim")
    out = W.Output()
    def redraw(*_):
        with out:
            out.clear_output(wait=True)
            fig, ax3, cursors = _viewer_figure(traj)
            draw_frame(ax3, sd, frame.value, elev.value, azim.value)
            for c in cursors: c.set_xdata([traj.t[frame.value]] * 2)
            plt.show()
    for wdg in (frame, elev, azim): wdg.observe(redraw, "value")
    redraw(); display(W.VBox([W.HBox([play, frame]), W.HBox([elev, azim]), out]))

KB, M_AIR, TORR = 1.380649e-23, 4.809e-26, 133.322

def gas_properties(T, p_torr):
    p = p_torr*TORR
    mu = 1.716e-5*(T/273.15)**1.5*(273.15 + 110.4)/(T + 110.4)
    kap = 0.00119 + 8.94685e-5*T - 2.18293e-8*T**2
    rho = p*M_AIR/(KB*T)
    return dict(T=T, p_torr=p_torr, p=p, mu=mu, kappa=kap, rho=rho, vbar=np.sqrt(8*KB*T/(np.pi*M_AIR)),
                lambda_mu=mu/p*np.sqrt(np.pi*KB*T/(2*M_AIR)), lambda_kappa=4*kap/(5*p)*np.sqrt(M_AIR*T/(2*KB)),
                kappa_tr=15/4*KB/M_AIR*mu, nu=mu/rho)

GAS_DEFAULT = gas_properties(260.0, 1.0)

def cunningham(kn):          return 1 + kn*(1.257 + 0.4*np.exp(-1.1/kn))

def rot_slip_length(gas):    return 4*gas["mu"]/(gas["rho"]*gas["vbar"])

def sphere_drag(a, gas):
    return 6*np.pi*gas["mu"]*a/cunningham(gas["lambda_mu"]/a), 8*np.pi*gas["mu"]*a**3/(1 + 3*rot_slip_length(gas)/a)

def f_T_2017(kn):            return 9*kn**3/((1 + 4.4844*kn**2)*(1 + kn))

def sphere_gamma(a, gas, model="APL2017", kp=2.2):
    if model == "APL2017":
        return f_T_2017(gas["lambda_kappa"]/a)*a**2*gas["kappa"]*gas["T"]/np.sqrt(2*KB*gas["T"]/M_AIR)
    if model == "Waldmann":
        return 32/15*a**2*gas["kappa_tr"]*gas["T"]/gas["vbar"]
    if model == "Talbot":
        kn, cs, ct, cm, kr = gas["lambda_mu"]/a, 1.17, 2.18, 1.14, gas["kappa"]/kp
        return 12*np.pi*gas["mu"]*gas["nu"]*a*cs*(kr + ct*kn)/((1 + 3*cm*kn)*(1 + 2*kr + 2*ct*kn))
    raise ValueError(f"unknown thermophoresis model {model!r}")

def ellipsoid_integrals(b):
    b2 = np.asarray(b, float)**2
    chi = 2*elliprf(*b2)
    al = np.array([2/3*elliprd(b2[(i+1) % 3], b2[(i+2) % 3], b2[i]) for i in range(3)])
    return chi, al

def ellipsoid_stokes(axes, mu):
    axes = np.asarray(axes, float); aeq = np.prod(axes)**(1/3); b = axes/aeq
    chi, al = ellipsoid_integrals(b)
    K = 16*np.pi*mu*aeq/(chi + b**2*al)
    Om = np.array([16*np.pi*mu*aeq**3*(b[j]**2 + b[k]**2)/(3*(b[j]**2*al[j] + b[k]**2*al[k]))
                   for j, k in [(1, 2), (2, 0), (0, 1)]])
    return K, Om

def prolate_coefficients(aspect):
    e = np.sqrt(1 - 1/aspect**2); L = np.log((1 + e)/(1 - e))
    return dict(XA=8/3*e**3/(-2*e + (1 + e**2)*L), YA=16/3*e**3/(2*e + (3*e**2 - 1)*L),
                XC=4/3*e**3*(1 - e**2)/(2*e - (1 - e**2)*L), YC=4/3*e**3*(2 - e**2)/(-2*e + (1 + e**2)*L))

def graded_ellipsoid_particle(semiaxes, epsilon=0.0, epsilon_a3=None, density=917.0, gas=None,
                              thermo_model="APL2017", kp=2.2, name=None):
    gas = GAS_DEFAULT if gas is None else gas
    ax = np.asarray(semiaxes, float); a1, a2, a3 = ax
    eps = epsilon_a3/a3 if epsilon_a3 is not None else float(epsilon)
    if abs(eps)*a3 >= 1: raise ValueError(f"density must stay positive: need |eps| a3 < 1 (got {eps*a3})")
    aeq = np.prod(ax)**(1/3)
    if np.allclose(ax, a1, rtol=1e-14):
        K0, Om0 = np.full(3, 6*np.pi*gas["mu"]*a1), np.full(3, 8*np.pi*gas["mu"]*a1**3)
    else:
        K0, Om0 = ellipsoid_stokes(ax, gas["mu"])
    kn = gas["lambda_mu"]/aeq
    K, Om = K0/cunningham(kn), Om0/(1 + 3*rot_slip_length(gas)/aeq)
    vT = sphere_gamma(aeq, gas, thermo_model, kp)/sphere_drag(aeq, gas)[0]
    res = resistance_from_blocks(np.diag(K), Z3, np.diag(Om))
    M = 4/3*np.pi*a1*a2*a3*density
    zc = eps*a3**2/5
    inertia = M/5*np.diag([a2**2 + a3**2, a3**2 + a1**2, a1**2 + a2**2]) - M*zc**2*np.diag([1, 1, 0])
    kind = "Graded ellipsoid" if eps != 0 else ("Sphere" if np.allclose(ax, a1) else "Ellipsoid")
    auto = f"{kind} {'/'.join(f'{1e6*a:.3g}' for a in ax)} µm" + (f", ε·a3 = {eps*a3:.3g}" if eps else "")
    return Particle(name or auto, M, inertia, res, vT*res[:, :3], rH=[0, 0, -zc],
                    shape={"type": "ellipsoid", "semiaxes": ax, "epsilon": eps},
                    meta=dict(kind=kind, aeq=aeq, Kn=kn, vT=vT, zc=zc, epsilon=eps, density=density,
                              thermo_model=thermo_model, gas=gas, length_scale=ax.max()))

def sphere_particle(radius, **kw):     return graded_ellipsoid_particle([radius]*3, epsilon=0.0, **kw)

def ellipsoid_particle(semiaxes, **kw): return graded_ellipsoid_particle(semiaxes, epsilon=0.0, **kw)

def righting_stiffness(p, g=9.81):     return p.mass*g*abs(p.meta["zc"])*p.thermo[0, 0]/p.thermo[2, 2]

PARTICLE_LIBRARY = {
    "Sphere":          dict(builder=sphere_particle, parameters=["radius", "density", "gas", "thermo_model", "kp"],
                            description="Homogeneous sphere of radius a (m)."),
    "Ellipsoid":       dict(builder=ellipsoid_particle, parameters=["semiaxes", "density", "gas", "thermo_model", "kp"],
                            description="Homogeneous triaxial ellipsoid, semi-axes (a1, a2, a3) in m along body axes."),
    "GradedEllipsoid": dict(builder=graded_ellipsoid_particle,
                            parameters=["semiaxes", "epsilon (1/m) or epsilon_a3", "density", "gas", "thermo_model", "kp"],
                            description="Ellipsoid with density rho (1 + eps z); COM displaced from the force centre."),
}

def register_particle_type(name, builder, parameters, description):
    PARTICLE_LIBRARY[name] = dict(builder=builder, parameters=parameters, description=description)

def build_particle(kind, **params):
    return PARTICLE_LIBRARY[kind]["builder"](**params)

def particle_summary(p):
    Rc, Tc = p.at_com(); m = p.meta
    rows = [("mass (kg)", f"{p.mass:.5e}"), ("principal moments about COM (kg m²)", np.array2string(np.sort(np.linalg.eigvalsh(p.inertia)), precision=4)),
            ("diag K (N s/m)", np.array2string(np.diag(p.res)[:3], precision=5)), ("diag Ω (N m s)", np.array2string(np.diag(p.res)[3:], precision=5)),
            ("diag Γ (N m)", np.array2string(np.diag(p.thermo[:3]), precision=5)), ("r_H: COM → force centre (m)", np.array2string(p.rH, precision=4)),
            ("Knudsen number (a_eq)", f"{m.get('Kn', float('nan')):.4g}"), ("thermophoretic velocity scale v_T (m²/s)", f"{m.get('vT', float('nan')):.4g}"),
            ("translational relaxation M/K (s)", f"{p.mass/p.res[0, 0]:.4g}"), ("rotational relaxation I/Ω (s)", f"{p.inertia[0, 0]/p.res[3, 3]:.4g}")]
    display(Markdown(f"**{p.name}**\n\n| quantity | value |\n|---|---|\n" + "\n".join(f"| {a} | {b} |" for a, b in rows)))
    show_matrix(Rc, r"$\mathcal R_{\rm COM}=$"); show_matrix(Tc, r"$\mathcal T_{\rm COM}=$")

def particle_diagram(p, ax=None, tilt=0.0):
    s = p.shape; a1, _, a3 = s["semiaxes"]; eps = s.get("epsilon", 0.0); L = max(a1, a3); u = 1e6
    ax = ax or plt.figure(figsize=(3.4, 4)).gca()
    c, sn = np.cos(np.radians(tilt)), np.sin(np.radians(tilt))
    rot = lambda X, Z: (c*X + sn*Z, -sn*X + c*Z)        # body x-z plane -> lab x-z plane, tilt about y
    X, Z = np.meshgrid(np.linspace(-a1, a1, 200), np.linspace(-a3, a3, 200))
    inside = (X/a1)**2 + (Z/a3)**2 <= 1
    dens = np.where(inside, 1 + eps*Z, np.nan)
    Xl, Zl = rot(X, Z)
    ax.pcolormesh(u*Xl, u*Zl, dens, cmap="Blues", vmin=1 - 1.6*max(abs(eps)*a3, 0.3), vmax=1 + 1.2*max(abs(eps)*a3, 0.3), shading="auto")
    t = np.linspace(0, 2*np.pi, 200); ex, ez = rot(a1*np.cos(t), a3*np.sin(t)); ax.plot(u*ex, u*ez, color="navy", lw=1)
    H, C = np.array(rot(0, 0)), np.array(rot(0, -p.rH[2]))
    ax.plot(*(u*H), "o", color="red", ms=5, label="force centre H"); ax.plot(*(u*C), "ko", ms=5, label="COM")
    if np.linalg.norm(p.rH) > 0:
        ax.annotate("", xy=u*H, xytext=u*C, arrowprops=dict(arrowstyle="->", color="k"))
    e_lo, e_hi = np.array(rot(0, -1.3*a3)), np.array(rot(0, 1.3*a3))
    ax.plot(*(u*np.stack([e_lo, e_hi], 1)), color="blue", ls=":", lw=1)
    ax.text(*(u*e_hi*1.04), r"$\hat e_3$", color="blue", ha="center")
    ax.annotate("", xy=u*(H + [0, 1.3*L]), xytext=u*H, arrowprops=dict(arrowstyle="-|>", color="green", lw=2))
    ax.annotate("", xy=u*(C - [0, 1.3*L]), xytext=u*C, arrowprops=dict(arrowstyle="-|>", color="dimgray", lw=2))
    ax.text(*(u*(H + [-0.1*L, 1.2*L])), r"$\mathbf{F}_T$", color="green", ha="right")
    ax.text(*(u*(C + [-0.1*L, -1.25*L])), r"$M\mathbf{g}$", color="dimgray", ha="right")
    ax.set_aspect("equal"); ax.set_xlabel("x (µm)"); ax.set_ylabel("z (µm)")
    ax.set_xlim(-1.6*u*L, 1.6*u*L); ax.set_ylim(-1.7*u*L, 1.8*u*L); ax.legend(fontsize=6, loc="lower right")
    ax.set_title(f"{p.name}\n$a_1$={u*a1:.3g} µm, $a_3$={u*a3:.3g} µm, $|z_c|$={u*abs(p.meta.get('zc', 0)):.3g} µm", fontsize=8)
    return ax
