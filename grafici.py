#!/usr/bin/env python3
"""
Genera le figure del README a partire dai csv prodotti dal solutore.

Uso:
    ./cavita N=64 gamma=1 tol=1e-9 verbose=0   &&  mv campo.csv campo_cds.csv
    ./cavita N=64 gamma=0 tol=1e-9 verbose=0   &&  mv campo.csv campo_uds.csv
    ./cavita N=64 --spazzate
    python3 grafici.py

N non e' scritto nel codice: si ricava dal numero di righe del campo.
Cosi' le figure restano valide se cambi griglia.
"""

import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")          # niente finestre, si salva e basta
import matplotlib.pyplot as plt

USCITA = "figures"


def leggi_campo(nome):
    """Legge campo.csv e restituisce le matrici (N,N) di u, v, psi."""
    u, v, psi = [], [], []
    for r in csv.DictReader(open(nome)):
        u.append(float(r["u"]))
        v.append(float(r["v"]))
        psi.append(float(r["psi"]))
    N = int(round(len(u) ** 0.5))
    if N * N != len(u):
        raise ValueError(nome + ": il numero di celle non e' un quadrato")
    forma = (N, N)                       # riga = i (direzione y), colonna = j
    return (np.array(u).reshape(forma),
            np.array(v).reshape(forma),
            np.array(psi).reshape(forma), N)


def leggi_parametri(nome="parametri.csv"):
    return next(csv.DictReader(open(nome)))


def centro(campo):
    """Valore sulla mediana del dominio: media delle due file centrali."""
    N = campo.shape[0]
    return 0.5 * (campo[:, N // 2 - 1] + campo[:, N // 2])


# ---------------------------------------------------------------- spazzate

def figura_spazzate(nome="spazzate.csv"):
    righe = [r for r in csv.DictReader(open(nome)) if r["fase"] == "nswP"]
    griglie = sorted({int(r["N"]) for r in righe})

    fig, assi = plt.subplots(1, len(griglie), figsize=(5.2 * len(griglie), 4.0))
    if len(griglie) == 1:
        assi = [assi]

    for ax, N in zip(assi, griglie):
        d = sorted((int(r["nswP"]), int(r["iter"]), float(r["tempo"]))
                   for r in righe if int(r["N"]) == N)
        x = [a for a, _, _ in d]
        it = [b for _, b, _ in d]
        t = [c for _, _, c in d]

        ax.plot(x, it, "o-", color="tab:blue", label="iterazioni esterne")
        ax.set_xlabel("spazzate SIP sulla pressione (nswP)")
        ax.set_ylabel("iterazioni esterne", color="tab:blue")
        ax.tick_params(axis="y", labelcolor="tab:blue")
        ax.set_title("N = %d" % N)
        ax.grid(alpha=0.3)

        ax2 = ax.twinx()
        ax2.plot(x, t, "s--", color="tab:red", label="tempo")
        ax2.set_ylabel("tempo di calcolo [s]", color="tab:red")
        ax2.tick_params(axis="y", labelcolor="tab:red")

        # il minimo del tempo, che e' il punto della figura
        imin = t.index(min(t))
        ax2.annotate("minimo: nswP = %d" % x[imin],
                     xy=(x[imin], t[imin]), xytext=(14, 42),
                     textcoords="offset points", fontsize=9, color="tab:red",
                     arrowprops=dict(arrowstyle="->", color="tab:red", lw=0.8))

    fig.tight_layout()
    fig.savefig(os.path.join(USCITA, "spazzate.png"), dpi=150)
    plt.close(fig)


# ------------------------------------------------------------------ schemi

def figura_profili(f_cds="campo_cds.csv", f_uds="campo_uds.csv"):
    uc, vc, _, N = leggi_campo(f_cds)
    uu, vu, _, _ = leggi_campo(f_uds)
    U = float(leggi_parametri()["U"])

    y = (np.arange(N) + 0.5) / N
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.4, 4.2))

    a1.plot(centro(uc) / U, y, "-", label="CDS (gamma = 1)")
    a1.plot(centro(uu) / U, y, "--", label="UDS (gamma = 0)")
    a1.set_xlabel("u / U"); a1.set_ylabel("y / L")
    a1.set_title("mediana verticale")

    a2.plot(y, centro(vc.T) / U, "-", label="CDS (gamma = 1)")
    a2.plot(y, centro(vu.T) / U, "--", label="UDS (gamma = 0)")
    a2.set_xlabel("x / L"); a2.set_ylabel("v / U")
    a2.set_title("mediana orizzontale")

    for a in (a1, a2):
        a.grid(alpha=0.3)
        a.legend()

    fig.suptitle("Diffusione numerica dell'upwind, N = %d" % N)
    fig.tight_layout()
    fig.savefig(os.path.join(USCITA, "schemi_profili.png"), dpi=150)
    plt.close(fig)


def figura_correnti(f_cds="campo_cds.csv", f_uds="campo_uds.csv"):
    _, _, pc, N = leggi_campo(f_cds)
    _, _, pu, _ = leggi_campo(f_uds)
    U = float(leggi_parametri()["U"])
    L = float(leggi_parametri()["L"])

    x = (np.arange(N) + 0.5) / N
    X, Y = np.meshgrid(x, x)

    # livelli fitti vicino allo zero: i vortici secondari sono debolissimi
    forti = np.linspace(-0.12, -0.005, 12)
    deboli = np.array([1e-6, 5e-6, 2e-5, 1e-4, 5e-4, 1.5e-3])
    livelli = np.concatenate([forti, [0.0], deboli])

    fig, assi = plt.subplots(1, 2, figsize=(10.0, 4.9))
    for ax, psi, tit in ((assi[0], pc / (U * L), "CDS (gamma = 1)"),
                         (assi[1], pu / (U * L), "UDS (gamma = 0)")):
        ax.contour(X, Y, psi, levels=livelli, colors="k", linewidths=0.7)
        ax.set_aspect("equal")
        ax.set_title("%s      psi_min = %.5f" % (tit, psi.min()))
        ax.set_xlabel("x / L"); ax.set_ylabel("y / L")

    fig.suptitle("Linee di corrente, N = %d, Re = %s"
                 % (N, leggi_parametri()["Re"]))
    fig.tight_layout()
    fig.savefig(os.path.join(USCITA, "linee_corrente.png"), dpi=150)
    plt.close(fig)


# -------------------------------------------------------------------- main

if __name__ == "__main__":
    os.makedirs(USCITA, exist_ok=True)

    if os.path.exists("spazzate.csv"):
        figura_spazzate()
        print("scritto", USCITA + "/spazzate.png")

    if os.path.exists("campo_cds.csv") and os.path.exists("campo_uds.csv"):
        figura_profili()
        print("scritto", USCITA + "/schemi_profili.png")
        figura_correnti()
        print("scritto", USCITA + "/linee_corrente.png")
