from matplotlib import pyplot as plt
import numpy as np


def initial_condition(x):
    u = np.zeros_like(x)
    mask = (x >= 0) & (x <= 2)
    u[mask] = np.sin(2 * np.pi * x[mask])
    return u


def exact_solution(x, t, c=1):
    x_shifted = x - c * t
    return initial_condition(x_shifted)


def upwind_method(c, h, dt, t_final, x_min=-1, x_max=5):
    nx = int((x_max - x_min) / h) + 1
    x = np.linspace(x_min, x_max, nx)
    u = initial_condition(x)
    nu = c * dt / h

    n_steps = int(t_final / dt)
    actual_t = n_steps * dt

    for n in range(n_steps):
        u_new = np.copy(u)
        for i in range(1, nx):
            u_new[i] = u[i] - nu * (u[i] - u[i - 1])
        u_new[0] = u[0] - nu * (u[0] - u[-1])
        u = u_new

    return x, u, actual_t, nu


c = 1.0
t_target = 1.0

# Case 1: h = 1/10, dt = 0.75h
h1 = 1 / 10
dt1 = 0.75 * h1
x1, u1_numerical, t1, nu1 = upwind_method(c, h1, dt1, t_target)
u1_exact = exact_solution(x1, t1, c)

x1_fine = np.linspace(-1, 5, 1000)
u1_exact_fine = exact_solution(x1_fine, t1, c)
# Case 2: h = 1/4, dt = 0.75h
h2 = 1 / 4
dt2 = 0.75 * h2
x2, u2_numerical, t2, nu2 = upwind_method(c, h2, dt2, t_target)
u2_exact = exact_solution(x2, t2, c)

x2_fine = np.linspace(-1, 5, 1000)
u2_exact_fine = exact_solution(x2_fine, t2, c)
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Case 1 plot
axes[0].plot(x1_fine, u1_exact_fine, "b-", label="Exact solution", linewidth=2)
axes[0].plot(
    x1,
    u1_numerical,
    "r--",
    label="Upwind method",
    linewidth=2,
    marker="o",
    markevery=5,
    markersize=4,
)
axes[0].set_xlabel("x", fontsize=12)
axes[0].set_ylabel("u(x,t)", fontsize=12)
axes[0].set_title(
    rf"Case 1: h = 1/10, $\Delta t$ = 0.75h\nt ≈ {t1:.4f}, $\nu$ = {nu1:.4f}", fontsize=11
)
axes[0].legend(fontsize=10)
axes[0].grid(True, alpha=0.3)
axes[0].set_xlim([-0.5, 3.5])

# Case 2 plot
axes[1].plot(x2_fine, u2_exact_fine, "b-", label="Exact solution", linewidth=2)
axes[1].plot(
    x2,
    u2_numerical,
    "r--",
    label="Upwind method",
    linewidth=2,
    marker="s",
    markevery=2,
    markersize=6,
)
axes[1].set_xlabel("x", fontsize=12)
axes[1].set_ylabel("u(x,t)", fontsize=12)
axes[1].set_title(
    f"Case 2: h = 1/4, $\Delta t$ = 0.75h\nt ≈ {t2:.4f}, $\nu$ = {nu2:.4f}", fontsize=11
)
axes[1].legend(fontsize=10)
axes[1].grid(True, alpha=0.3)
axes[1].set_xlim([-0.5, 3.5])

plt.tight_layout()
plt.savefig("graph/Q4A.svg", dpi=300, bbox_inches="tight")
