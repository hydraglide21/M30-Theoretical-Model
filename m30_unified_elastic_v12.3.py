"""
===============================================================================
M30 UNIFIED ELASTIC MODEL — VERSION V12.3 FINAL (KEATING FLEXION, ISO A=1)
===============================================================================
Validation multi-audit (Mistral, Claude, Meta) du socle dynamique et
thermodynamique.

Caractéristiques théoriques :
-----------------------------
1. Matrice dynamique D(q) : Keating-type angular flexion (12 triplets à 90°)
   garantissant l'invariance par rotation rigide et l'annulation stricte à Gamma.
2. Paramètres du Nombre d'Or (phi = (1+sqrt(5))/2) :
   - k2 = 1/phi ≈ 0.618034  (Ressorts diagonaux NNN)
   - kth = 1/(8*phi^2) ≈ 0.047746  (Flexion angulaire 90°)
   - kc = k2 + 8*kth = 1.000000  (Ressorts centraux NN)
3. Constantes élastiques macroscopiques exactes :
   - C11 = kc + 2*k2 = sqrt(5) ≈ 2.236068
   - C12 = k2 = 1/phi ≈ 0.618034
   - C44 = k2 + 4*kth = phi/2 ≈ 0.809017
4. Isotropie de Zener parfaite :
   - A = 2*C44 / (C11 - C12) = 1.000000
5. Stabilité spectrale (PSD) :
   - D(q) >= 0 sur toute la Zone de Brillouin 3D (pas de ramollissement).
6. Barrière de compression thermodynamique :
   - W(0) = 0, sigma(0) = 0 exacts à l'origine.
===============================================================================
"""

import numpy as np
import itertools

# --- CONSTANTES FONDAMENTALES ET GEOMETRIQUES ---
PHI = (1.0 + np.sqrt(5.0)) / 2.0  # Nombre d'or (~1.6180339887)

# Paramètres microscopiques issus du Nombre d'Or
K2 = 1.0 / PHI                  # ~0.6180339887 (NNN - ressorts diagonaux)
KTH = 1.0 / (8.0 * PHI**2)      # ~0.0477457514 (Flexion angulaire 90°)
KC = K2 + 8.0 * KTH             # = 1.0000000000 (NN - ressorts centraux, condition A=1)

# Voisins NN (6) et NNN (12)
NN = [np.array(v, float) for v in [(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]]
NNN = [np.array(v, float) for v in itertools.product([-1, 0, 1], repeat=3) if sum(abs(x) for x in v) == 2]

# 12 triplets orthogonaux à 90° par site pour la flexion angulaire de Keating
TRI = [(a, b) for a, b in itertools.combinations(NN, 2) if abs(a @ b) < 1e-9]
assert len(TRI) == 12, "Erreur : il doit y avoir exactement 12 triplets à 90°."


def dynamical_matrix(q, kc=KC, k2=K2, kth=KTH):
    """
    Calcule la matrice dynamique D(q) 3x3 avec phase de Bloch exacte,
    interactions centrales (NN, NNN) et flexion angulaire invariant par rotation.
    """
    M = np.zeros((3, 3), dtype=complex)
    
    # 1. Contributions des ressorts centraux (NN + NNN)
    for R, k in [(r, kc) for r in NN] + [(r, k2) for r in NNN]:
        Rh = R / np.linalg.norm(R)
        M += k * (1.0 - np.cos(q @ R)) * np.outer(Rh, Rh)
        
    # 2. Contribution de flexion angulaire de Keating (Triplets a, b à 90°)
    # g = -[a * e^(i q.b) + b * e^(i q.a) - (a + b)]
    for a, b in TRI:
        g = -(a * np.exp(1j * (q @ b)) + b * np.exp(1j * (q @ a)) - (a + b))
        M += kth * np.outer(g, np.conj(g))
        
    return M


def potential_compression(delta, kc=KC, cb=0.25):
    """
    Potentiel de compression W(delta) et contrainte sigma(delta)
    avec mur de repousse phenomenologique en 1/(1-delta).
    """
    if delta < 0 or delta >= 1.0:
        raise ValueError("La deformation delta doit etre dans [0, 1[.")
    
    # Densite d'energie W(delta)
    W = 0.5 * kc * delta**2 - cb * (delta + np.log(1.0 - delta))
    
    # Contrainte sigma(delta) = dW/ddelta
    sigma = kc * delta + cb * (delta / (1.0 - delta))
    
    return W, sigma


def run_full_audit():
    print("=================================================================")
    print("    AUDIT NUMERIQUE — MODEL M30 UNIFIE V12.3 (KEATING / ISO A=1)")
    print("=================================================================\n")
    
    # --- 1. VERIFICATION DES CONSTANTES ELASTIQUES EN GRANDES LONGUEURS D'ONDE ---
    e = 1e-4
    Dx = dynamical_matrix(np.array([e, 0.0, 0.0])).real / e**2
    Dxy = dynamical_matrix(np.array([e, e, 0.0]) / np.sqrt(2.0)).real / e**2
    
    C11 = Dx[0, 0]
    C44 = Dx[1, 1]
    C12 = 2.0 * Dxy[0, 1] - C44
    A_zener = (2.0 * C44) / (C11 - C12)
    
    print("--- 1. CONSTANTES ELASTIQUES MACROSCOPIQUES ---")
    print(f"  C11       = {C11:.6f}  (Attendu : sqrt(5)  = {np.sqrt(5):.6f})")
    print(f"  C12       = {C12:.6f}  (Attendu : 1/phi    = {1/PHI:.6f})")
    print(f"  C44       = {C44:.6f}  (Attendu : phi/2    = {PHI/2:.6f})")
    print(f"  A (Zener) = {A_zener:.6f}  (Attendu : ISO A=1  = 1.000000)")
    assert abs(A_zener - 1.0) < 1e-7, "Erreur : A doit etre egal a 1.0 !"
    print("  [OK] Isotropie de Zener parfaite confirmee.\n")

    # --- 2. VERIFICATION DU SPECTRE A GAMMA ET DANS LA ZONE DE BRILLOUIN 3D ---
    print("--- 2. STABILITE SPECTRALE ET ZONE DE BRILLOUIN ---")
    
    # Centre de zone Gamma
    D_gamma = dynamical_matrix(np.array([0.0, 0.0, 0.0]))
    eigs_gamma = np.sort(np.linalg.eigvalsh(D_gamma))
    print(f"  Valeurs propres a Gamma (q=0) : {eigs_gamma[0]:.2e}, {eigs_gamma[1]:.2e}, {eigs_gamma[2]:.2e}")
    assert np.allclose(eigs_gamma, 0.0, atol=1e-12), "Erreur : D(Gamma) doit s'annuler !"
    print("  [OK] Somme acoustique stricte respectee à Gamma.")
    
    # Hermiticite
    D_test = dynamical_matrix(np.array([0.3, 0.5, 0.7]))
    hermiticity_error = np.max(np.abs(D_test - np.conj(D_test.T)))
    print(f"  Ecart max Hermiticite (D - D†) : {hermiticity_error:.2e}")
    assert hermiticity_error < 1e-14, "Erreur : D(q) doit etre hermitique !"
    
    # Balayage de la Zone de Brillouin 3D (Grille 15^3)
    grid = np.linspace(-np.pi, np.pi, 15)
    min_eig_bz = 1e9
    for qx in grid:
        for qy in grid:
            for qz in grid:
                q = np.array([qx, qy, qz])
                if np.sum(np.abs(q)) < 1e-8:
                    continue
                w = np.linalg.eigvalsh(dynamical_matrix(q))[0]
                if w < min_eig_bz:
                    min_eig_bz = w
                    
    print(f"  Valeur propre minimale sur BZ 15^3 : {min_eig_bz:.6f} (> 0)")
    assert min_eig_bz >= -1e-12, "Erreur : Instabilite spectrale detectee !"
    print("  [OK] Positivite spectrale (PSD) et stabilite 3D confirmees.\n")

    # --- 3. BARRIERE DE COMPRESSION THERMODYNAMIQUE ---
    print("--- 3. BARRIERE DE COMPRESSION THERMODYNAMIQUE ---")
    deltas = [0.0, 0.20, 0.50, 0.95]
    for d in deltas:
        W, sigma = potential_compression(d)
        print(f"  delta = {d:.2f} | W(delta) = {W:10.6f} | sigma(delta) = {sigma:10.6f}")
        
    W0, sig0 = potential_compression(0.0)
    assert abs(W0) < 1e-12 and abs(sig0) < 1e-12, "Erreur : State de repos non nul !"
    print("  [OK] Repos W(0)=0 et sigma(0)=0 exacts. Mur de repousse d'accord avec Meta.\n")

    print("=================================================================")
    print("    VERDICT FINAL : SCRIPT V12.3 ENT IEREMENT VALIDE ET CONFORME")
    print("=================================================================")


if __name__ == "__main__":
    run_full_audit()