"""Calibrated 2D echo matched filter for gravitationally lensed FRBs (CHIME/FRB Catalog 2 dynamic spectra).
Lens model: X(nu, t+d) = X(nu, t)/R for the whole main pulse (achromatic, same morphology, same DM).
Statistic for each delay d: echo amplitude a(d) = <Y_W, Y_{W+d}>/<Y_W,Y_W> (Y whitened per sub-band), significance
S = a/sigma_a, and a shape test chi2 of the residual Y_{W+d} - a Y_W. Detection = max_d S with chi2/dof < cmax.
"""
import numpy as np, h5py
B = '/home/claude/pbh/data_zhou/Z/FRB_data/canfar_downloads/'

def load_sub(path, nsub=64, D=None):
    if D is None:
        with h5py.File(path, 'r') as f:
            D = np.array(f['data'], float); G = np.array(f['good_freq']) & np.array(f['calibration/good_freq']); FL = np.array(f['flag'])
        D[~FL] = np.nan; D[~G] = np.nan
    nf, nt = D.shape; k = nf // nsub
    X = np.nanmean(D[:nsub * k].reshape(nsub, k, nt), 1)
    good = np.isfinite(X).all(1); X = X[good]
    return X

def whiten(X, on):
    off = ~on
    mu = X[:, off].mean(1, keepdims=True); sd = X[:, off].std(1, keepdims=True); sd[sd == 0] = 1
    return (X - mu) / sd

def main_window(Y):
    ts = Y.sum(0); n = len(ts); noise = ts[: n // 3].std() + 1e-9
    pk = int(np.argmax(ts)); lo = pk; hi = pk
    while lo > 0 and ts[lo - 1] > 3 * noise: lo -= 1
    while hi < n - 1 and ts[hi + 1] > 3 * noise: hi += 1
    return lo, hi, ts[lo:hi + 1].sum() / (noise * np.sqrt(hi - lo + 1))

def echo_scan(X, dmin=2, dmax=70, cmax=2.0):
    nt = X.shape[1]
    ts0 = np.nansum(X, 0); pk = int(np.argmax(ts0)); on = np.zeros(nt, bool); on[max(pk - 8, 0):min(pk + 75, nt)] = True
    Y = whiten(X, on); lo, hi, snr = main_window(Y); w = hi - lo + 1
    best = dict(S=-np.inf, d=0, a=0, chi=np.inf, snr=snr, w=w)
    for d in range(max(dmin, w), dmax + 1):
        if hi + d >= nt: break
        A = Y[:, lo:hi + 1]; Bk = Y[:, lo + d:hi + 1 + d]
        # region between end of main pulse and start of echo must not be counted twice: require echo window outside main
        aa = np.sum(A * A); ab = np.sum(A * Bk); a = ab / aa
        sig_a = np.sqrt(1.0 / aa)                         # whitened unit noise
        r = Bk - a * A; dof = r.size - 1
        chi = np.sum(r * r) / (1 + a * a) / dof
        S = a / sig_a
        if 0 < a < 1 and chi < cmax and S > best['S']:
            best = dict(S=S, d=d, a=a, chi=chi, snr=snr, w=w)
    return best
