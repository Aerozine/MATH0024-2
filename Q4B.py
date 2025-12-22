from matplotlib import pyplot as plt
import numpy as np


def von_neumann_analysis(c, h, dt, n_points=500):

    nu = c * dt / h
    xi_max = np.pi / h
    xi = np.linspace(-xi_max, xi_max, n_points)
    gamma = 1 - nu * (1 - np.exp(-1j * xi * h))
    e_a = np.abs(gamma)
    omega = -np.angle(gamma) / dt
    e_d = np.zeros_like(xi)
    # prevent NaN by dividing by 0
    mask = xi != 0
    e_d[mask] = omega[mask] / (xi[mask] * c)
    e_d[~mask] = 1

    return xi, e_a, e_d


c = 1.0
t_target = 1.0

h1 = 1 / 10
dt1 = 0.75 * h1

h2 = 1 / 4
dt2 = 0.75 * h2
xi1, ea1, ed1 = von_neumann_analysis(c, h1, dt1)
xi2, ea2, ed2 = von_neumann_analysis(c, h2, dt2)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

axes[0, 0].plot(xi1, ea1, "b-", linewidth=2,label=r'$e_a$')
axes[0, 0].set_xlabel(r"Wavenumber $\xi$", fontsize=12)
axes[0, 0].set_ylabel(r"Dissipation error $e_a = |\gamma(\xi)|$", fontsize=12)
axes[0, 0].set_title(f"Case 1")
axes[0, 0].grid(True, alpha=0.3)
axes[0, 0].set_ylim([0, 1.1])
axes[0, 0].axhline(y=1,color='gray',linestyle='--',linewidth=1,alpha=0.7,label='ideal case')
axes[0, 0].set_xticks([-np.pi/h1, -2*np.pi, 0, 2*np.pi, np.pi/h1])
axes[0, 0].set_xticklabels([r'$-\frac{\pi}{h}$',r'$-2\pi$',r'$0$',r'$2\pi$',r'$\frac{\pi}{h}$'])
axes[0, 0].legend()

axes[0, 1].plot(xi1, ed1, "r-", linewidth=2,label=r'$e_d$')
axes[0, 1].set_xlabel(r"Wavenumber $\xi$", fontsize=12)
axes[0, 1].set_ylabel(r"Dispersion error $e_d = \frac{\omega}{\xi c}$", fontsize=12)
axes[0, 1].set_title(f"Case 1")
axes[0, 1].grid(True, alpha=0.3)
axes[0, 1].set_ylim([0, 1.5])
axes[0, 1].axhline(y=1,color='gray',linestyle='--',linewidth=1,alpha=0.7,label='ideal case')
axes[0, 1].set_xticks([-np.pi/h1, -2*np.pi, 0, 2*np.pi, np.pi/h1])
axes[0, 1].set_xticklabels([r'$-\frac{\pi}{h}$',r'$-2\pi$',r'$0$',r'$2\pi$',r'$\frac{\pi}{h}$'])
axes[0, 1].legend()

axes[1, 0].plot(xi2, ea2, "b-", linewidth=2,label=r'$e_a$')
axes[1, 0].set_xlabel(r"Wavenumber $\xi$", fontsize=12)
axes[1, 0].set_ylabel(r"Dissipation error $e_a = |\gamma(\xi)|$", fontsize=12)
axes[1, 0].set_title(f"Case 2")
axes[1, 0].grid(True, alpha=0.3)
axes[1, 0].set_ylim([0, 1.1])
axes[1, 0].axhline(y=1,color='gray',linestyle='--',linewidth=1,alpha=0.7,label='ideal case')
axes[1, 0].set_xticks([-np.pi/h2, -2*np.pi, 0, 2*np.pi, np.pi/h2])
axes[1, 0].set_xticklabels([r'$-\frac{\pi}{h}$',r'$-2\pi$',r'$0$',r'$2\pi$',r'$\frac{\pi}{h}$'])
axes[1, 0].legend()

axes[1, 1].plot(xi2, ed2, "r-", linewidth=2,label=r'$e_b$')
axes[1, 1].set_xlabel(r"Wavenumber $\xi$", fontsize=12)
axes[1, 1].set_ylabel(r"Dispersion error$e_d = \frac{\omega}{\xi c}$", fontsize=12)
axes[1, 1].set_title(f"Case 2")
axes[1, 1].grid(True, alpha=0.3)
axes[1, 1].set_ylim([0, 1.5])
axes[1, 1].axhline(y=1,color='gray',linestyle='--',linewidth=1,alpha=0.7,label='ideal case')
axes[1, 1].set_xticks([-np.pi/h2, -2*np.pi, 0, 2*np.pi, np.pi/h2])
axes[1, 1].set_xticklabels([r'$-\frac{\pi}{h}$',r'$-2\pi$',r'$0$',r'$2\pi$',r'$\frac{\pi}{h}$'])
axes[1, 1].legend()

plt.tight_layout()
plt.savefig("graph/Q4B.svg", dpi=300, bbox_inches="tight")
plt.show()
