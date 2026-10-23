 # Model $M_{30}$ — Unification Géométrique, Spectrale & Cosmologique

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status: V10.1 Audited](https://img.shields.io/badge/Status-V10.1_Audited-brightgreen.svg)]()
[![GitHub Pages](https://img.shields.io/badge/Interactive_App-GitHub_Pages-blue.svg)](https://hydraglide21.github.io/M30-Theoretical-Model/)

Bienvenue sur le dépôt officiel du cadre théorique **$M_{30}$**. Ce projet explore l'unification des modèles de réseaux discrets, de l'élasticité quasicristalline, des effets topologiques non-hermitiens, de l'hydrodynamique péristaltique et de la cosmologie quantique à rebond (LQC).

---

## 📌 Aperçu de la Version 10.1 (Audite & Assainie)

La version **V10.1** représente la version de référence auditée par exécution indépendante des scripts. L'ensemble des résultats chiffrés est adossé à des exécutions reproductibles.

### Les 6 Pistes du Modèle $M_{30}$

| Piste | Domaine | Observable Clé | Statut d'Audit |
| :--- | :--- | :--- | :--- |
| **A** | Élasticité Quasicristalline | $C_{11}=\sqrt{5}$, $C_{44}=\varphi$, Zener $A=2.0$ | Archivé (non réexécuté) |
| **B** | Topologie 1D NH | $W=1$, Skin $\ge 98.45\%$, $V_{0,c}=2$ | **Validé par exécution** |
| **C** | Corner Skin 2D | $80.5\%$ quadrant $5\times 5$, Jauge $S_{2D}$ trivial | **Validé par exécution** |
| **D** | Lubrification | $\overline{Q}=0.2088$ (+25.3%), $C=-1.3827$ | **Validé** (Réserve $kH_0$) |
| **E** | Advection-Diffusion | $\alpha=\frac{1}{2}\ln(1+Pe)$, $G_{\max}\approx 1.005$ | **Validé par exécution** |
| **F** | Rebond LQC | $R_b=40.59$, Transmissions $\approx 1$, $n_s=-0.0358$ | **Validé avec réserves** |

> **Note sur l'Unification :** La jauge de similitude $S = \operatorname{diag}(e^{-\alpha x})$ assure l'unité conceptuelle entre les pistes B, C et E. Le lien avec la cosmologie LQC (Piste F) est qualifié comme une **analogie structurelle**.

---

## 🔮 Perspective : Version 11 — Propagation des Ondes Gravitationnelles

Dans le cadre du modèle $M_{30}$, **au point critique ($\rho = \rho_c$), rien ne s'effondre.** Le rebond non singulier ($a_{\text{min}} > 0, H = 0$) permet d'étudier la propagation continue des perturbations tensorielles à travers la haute densité.

La **Version 11** est dédiée à la dérivation des équations de propagation des **Ondes Gravitationnelles (OG) primordiales** :

$$\ddot{h}_{k} + 3H \dot{h}_{k} + \left( \frac{k^2}{a^2} + m_{\text{eff}}^2 \right) h_{k} = 0$$

- **Absence de divergence** : Le potentiel géométrique $U(\tau) = a''/a$ reste strictement fini.
- **Transmission spectrale** : Analyse de la signature stochastique des OG et des modes B du CMB.

---

## 📁 Structure du Dépôt

- `index.html` : Interface web interactive (compatible GitHub Pages).
- `Rapport-Final-M30-V10.1.pdf` : Rapport complet de l'audit indépendant.
- `audit_m30_v72.py` : Script Python de simulation et de vérification spectrale (Pistes B & C).
- `M30_PisteE_audit_multiIA.py` : Script d'advection-diffusion (Piste E).

---

## 🚀 Consulter le Site Interactif

L'application web dynamique permettant d'explorer les 6 pistes et la synthèse du modèle est accessible directement via GitHub Pages :  
👉 **[Accéder à l'application $M_{30}$](https://hydraglide21.github.io/M30-Theoretical-Model/)** *(Remplacez par votre lien réel)*

---
*Projet open-source sous licence MIT.*

