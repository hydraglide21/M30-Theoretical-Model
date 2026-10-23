"""
Piste E corrigée — advection-diffusion linéaire 1D (U constant, Dirichlet) <-> Hatano-Nelson.
Convention : dψ_j/dt = c_left ψ_{j-1} + c_diag ψ_j + c_right ψ_{j+1}   (U>0 : écoulement vers +x)
Symétrisation : L_sym = S L S^-1, S = diag(e^{-α j}), α = ½ ln(c_left/c_right), exige c_left*c_right > 0.
Vecteurs propres droits de L : v = S^-1 u = e^{+α j} u  -> localisés à l'AVAL (j grand).
"""
import numpy as np
from scipy.linalg import expm

def operator(N, U, nu, dx, scheme):
    if scheme == "upwind":
        cl, cr, cd = U/dx + nu/dx**2, nu/dx**2, -(U/dx + 2*nu/dx**2)
    else:  # central
        cl, cr, cd = U/(2*dx) + nu/dx**2, -U/(2*dx) + nu/dx**2, -2*nu/dx**2
    L = np.diag(cd*np.ones(N)) + np.diag(cr*np.ones(N-1), 1) + np.diag(cl*np.ones(N-1), -1)
    return L, cl, cr, cd

def audit(N=80, U=2.0, nu=0.05, dx=1.0, scheme="upwind"):
    L, cl, cr, cd = operator(N, U, nu, dx, scheme)
    Pe = U*dx/nu
    print(f"\n--- {scheme} N={N} U={U} nu={nu} dx={dx}  Pe_cellule={Pe:.2f} ---")
    if cl*cr <= 0:
        print(f"c_left*c_right = {cl*cr:.3f} <= 0 : NON symétrisable par jauge réelle (spectre complexe).")
        ev = np.linalg.eigvals(L); print(f"max|Im(E)| = {np.abs(ev.imag).max():.3f}")
        return
    alpha = 0.5*np.log(cl/cr); j = np.arange(1, N+1)
    S = np.exp(-alpha*j)
    Ls = (S[:, None]*L)/S[None, :]
    print(f"alpha={alpha:.4f}  xi=1/alpha={1/alpha:.2f} mailles")
    print(f"asymétrie de L_sym : {np.abs(Ls-Ls.T).max():.2e}  (symétrisation exacte)")
    E, Uv = np.linalg.eigh((Ls+Ls.T)/2)                 # spectre réel exact
    ev = np.linalg.eigvals(L)
    print(f"spectre exact (eigh) : [{E.min():.3f}, {E.max():.3f}], réel par construction")
    print(f"eig dense sur L : max|Im| = {np.abs(ev.imag).max():.2e}  <- artefact de non-normalité")
    V = Uv/S[:, None]; V /= np.linalg.norm(V, axis=0)    # vecteurs propres droits de L
    d = np.abs(V)**2
    print(f"masse moyenne 10 premiers sites = {100*d[:10].sum(0).mean():.2f} %,  10 derniers (aval) = {100*d[-10:].sum(0).mean():.2f} %")
    print(f"conditionnement des vecteurs propres ~ e^(2 alpha (N-1)) : log10 = {2*alpha*(N-1)/np.log(10):.0f}")
    gains = [np.linalg.norm(expm(L*t), 2)**2 for t in np.linspace(0, 15, 60)]
    print(f"gain transitoire max ||e^(Lt)||^2 = {max(gains):.4f}  (L dissipatif : pas d'amplification)")
    fro = np.linalg.norm(L@L.T - L.T@L)
    L2 = operator(2*N, U, nu, dx, scheme)[0]
    print(f"||[L,L^T]||_F = {fro:.4f} (N={N}) ; {np.linalg.norm(L2@L2.T-L2.T@L2):.4f} (N={2*N}) : effet de bord, indépendant de N")

if __name__ == "__main__":
    for nu in [0.05, 0.5, 0.8]:
        audit(nu=nu, scheme="upwind")
    audit(nu=0.05, scheme="central")     # Pe=40 > 2 : signes mixtes
    audit(nu=1.5, scheme="central")      # Pe=1.33 < 2 : symétrisable aussi
