"""Toy inverse problem with exact ground truth: does a regularized reconstruction invent localized
'perturbers' (false detections), and does a noise-free (Asimov) search predict when?

Model (all linear-Gaussian, so the Bayesian evidence is exact):
  data  y = B W(p, a) s + n         s: true source image (N x N), B: Gaussian PSF blur,
                                    W(p,a): local deflection of strength a centred at p ("perturber"), n: white noise
  prior s ~ N(0, (lam * D^T D + eps I)^-1)   D: finite-difference gradient -> smooth sources preferred
  evidence  ln Z = ln N(y | 0, C),  C = sigma^2 I + M P^-1 M^T,  M = B W
Detection score at position p:  d(p) = max_a max_lam ln Z(p,a)  -  max_lam ln Z(no perturber)
"""
import numpy as np

N = 16                      # image is N x N pixels
SIGMA_PSF = 1.0             # PSF width (pixels)
LAMS = (0.5, 3.0, 15.0)   # regularization strengths searched (like re-fitting 'reg')
AMPS = (0.6, 1.2, 2.4)      # perturber strengths searched (deflection scale, pixels)
EPS = 1e-3

yy, xx = np.mgrid[0:N, 0:N].astype(float)

def blur_matrix(sig=SIGMA_PSF):
    r = np.arange(N)
    g = np.exp(-0.5 * ((r[:, None] - r[None, :]) / sig) ** 2); g /= g.sum(1, keepdims=True)
    return np.kron(g, g)                         # separable 2D Gaussian blur, (N^2 x N^2)

def grad_prior():
    I = np.eye(N); Dx = np.eye(N, k=1) - I; Dx[-1] = 0
    D = np.vstack([np.kron(I, Dx), np.kron(Dx, I)])
    return D.T @ D

B = blur_matrix(); DTD = grad_prior(); EYE = np.eye(N * N)

def warp_matrix(p, a, core=1.5):
    """Bilinear-interpolation operator for a radial 'point-mass-like' deflection centred at p:
    observed pixel at r shows the source at r - alpha(r), alpha = a * (r-p)/|r-p|^2 * (1 - exp(-|r-p|^2 / 2 core^2))."""
    if a == 0: return EYE
    dy, dx = yy - p[0], xx - p[1]; r2 = dy ** 2 + dx ** 2 + 1e-9
    f = a * (1 - np.exp(-r2 / (2 * core ** 2))) / r2
    sy, sx = np.clip(yy - f * dy, 0, N - 1.001), np.clip(xx - f * dx, 0, N - 1.001)
    y0, x0 = np.floor(sy).astype(int), np.floor(sx).astype(int); ty, tx = sy - y0, sx - x0
    W = np.zeros((N * N, N * N)); rows = np.arange(N * N)
    for oy, ox, w in ((0, 0, (1 - ty) * (1 - tx)), (1, 0, ty * (1 - tx)), (0, 1, (1 - ty) * tx), (1, 1, ty * tx)):
        np.add.at(W, (rows, ((y0 + oy) * N + (x0 + ox)).ravel()), w.ravel())
    return W

def log_evidence(y, M, sigma, lam):
    P = lam * DTD + EPS * EYE
    C = sigma ** 2 * EYE + M @ np.linalg.solve(P, M.T)
    L = np.linalg.cholesky(C); z = np.linalg.solve(L, y)
    return -0.5 * z @ z - np.log(np.diag(L)).sum() - 0.5 * len(y) * np.log(2 * np.pi)

def best_ev(y, M, sigma): return max(log_evidence(y, M, sigma, l) for l in LAMS)

def search(y, sigma, positions, Wcache):
    """Detection map: d(p) = best evidence with a perturber at p minus best smooth evidence."""
    z0 = best_ev(y, B, sigma)
    return np.array([max(best_ev(y, B @ Wcache[(i, a)], sigma) for a in AMPS) - z0 for i in range(len(positions))])

def random_source(rng, n_smooth=4, n_sharp=0, sharp_w=0.6):
    """Smooth blobs plus optional compact (sharp) features that a smoothness prior represents poorly."""
    s = np.zeros((N, N))
    for _ in range(n_smooth):
        cy, cx = rng.uniform(5, N - 5, 2); w = rng.uniform(2.0, 4.0)
        s += rng.uniform(0.5, 1.0) * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * w ** 2))
    for _ in range(n_sharp):
        cy, cx = rng.uniform(5, N - 5, 2)
        s += rng.uniform(0.8, 1.5) * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * sharp_w ** 2))
    return s.ravel()
