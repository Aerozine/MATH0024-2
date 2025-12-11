import numpy as np
import matplotlib.pyplot as plt

lambda_val = -5
t_max = 1.5
u0 = 1
t = np.linspace(0, t_max, 100)
u_exact = u0 * np.exp(lambda_val * t)


def forward_euler(delta_t, lambda_val, u0, t_max):
    n_steps = int(t_max / delta_t)
    u = np.zeros(n_steps)
    u[0] = u0
    for n in range(1, n_steps):
        u[n] = (1 + delta_t * lambda_val) * u[n - 1]
    return u

# stable plot
fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(9, 6))
for dt in [0.01, 0.15, 0.2, 0.3]:
    u_stable = forward_euler(dt, lambda_val, u0, t_max)
    ttime = np.linspace(0, t_max, len(u_stable))
    ax1.plot(ttime, u_stable, label=rf"Stable ($\Delta t$ = {dt})")
    ax2.plot(
        ttime,
        np.abs(u_stable - u0 * np.exp(lambda_val * ttime)),
        label=rf"Stable ($\Delta t$ = {dt})",
    )

ax1.plot(t, u_exact, label="Exact Solution", color="black", linestyle="--", linewidth=2)
ax1.grid(True)
ax1.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
ax2.grid(True)
plt.tight_layout()
plt.savefig("graph/Q2S.svg")
# unstable plot 
t_max = 4
fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(9, 6))
for dt in [0.4, 0.5, 0.6]:
    u_stable = forward_euler(dt, lambda_val, u0, t_max)
    ttime = np.linspace(0, t_max, len(u_stable))
    ax1.plot(ttime, u_stable, label=rf"Unstable ($\Delta t$ = {dt})")
    ax2.plot(
        ttime,
        np.abs(u_stable - u0 * np.exp(lambda_val * ttime)),
        label=rf"Stable ($\Delta t$ = {dt})",
    )

ax1.plot(t, u_exact, label="Exact Solution", color="black", linestyle="--", linewidth=2)

ax1.grid(True)
ax1.legend(
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    borderaxespad=0,
)
ax2.grid(True)
ax1.set_xlabel("Time t")
ax1.set_ylabel("y")
ax2.set_xlabel("Time t")
ax2.set_ylabel("absolute error")
plt.tight_layout()
plt.savefig("graph/Q2U.svg")
