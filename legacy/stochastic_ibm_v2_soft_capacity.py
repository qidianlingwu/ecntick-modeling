import numpy as np
import matplotlib.pyplot as plt

# ===== 参数（平滑生长版本） =====
mu_max = 0.8            # 最大生长速率
K_soft = 6000           # 软性承载（生长速率减半时的菌量）
r_plasmid = 0.7         # 质粒复制速率
d_plasmid = 0.1
PCN_threshold = 30
f_survive = 0.15
lysis_delay = 0.5
N_A = 1000
MAX_PCN = 500

n_init = 100
pcns_init = np.random.poisson(15, n_init).astype(np.int64)
dt = 0.02
T_end = 100
steps = int(T_end / dt)

times = np.linspace(0, T_end, steps+1)
pop_h = np.zeros(steps+1)
pcn_h = np.zeros(steps+1)
lysis_events = []
nanobody = 0.0
nb_h = np.zeros(steps+1)

pcn = pcns_init.copy()
n_cells = n_init
lysis_timer = np.zeros(n_cells)

for i, t in enumerate(times):
    pop_h[i] = n_cells
    pcn_h[i] = np.mean(pcn) if n_cells > 0 else 0
    nb_h[i] = nanobody
    if n_cells == 0:
        break

    # 1. 分裂：软性承载
    crowd = 1.0 / (1.0 + n_cells / K_soft)
    div_prob = mu_max * dt * crowd / (1 + pcn / PCN_threshold)
    div_prob = np.clip(div_prob, 0, 1)
    will_divide = np.random.binomial(1, div_prob)
    n_div = will_divide.sum()
    if n_div > 0:
        dividing_idx = np.where(will_divide)[0]
        new_pcn = np.zeros(n_div, dtype=np.int64)
        for j, idx in enumerate(dividing_idx):
            mother = pcn[idx]
            d1 = np.random.binomial(mother, 0.5) if mother > 0 else 0
            d2 = mother - d1
            pcn[idx] = d1
            new_pcn[j] = d2
        pcn = np.concatenate([pcn, new_pcn])
        n_cells = len(pcn)
        lysis_timer = np.zeros(n_cells)

    # 2. 质粒复制与丢失
    if n_cells > 0:
        rep = np.random.binomial(pcn.astype(np.int64), r_plasmid * dt)
        pcn = pcn + rep
        loss = np.random.binomial(pcn.astype(np.int64), d_plasmid * dt)
        pcn = pcn - loss
        pcn = np.clip(pcn, 1, MAX_PCN)

    # 3. 裂解触发
    in_risk = pcn >= PCN_threshold
    lysis_timer[in_risk] += dt
    lysis_timer[~in_risk] = 0
    to_lyse = lysis_timer >= lysis_delay
    if np.any(to_lyse):
        n_lysing = np.sum(to_lyse)
        lysis_events.append((t, n_lysing))
        survive = np.random.rand(n_lysing) < f_survive
        surviving_idx = np.where(to_lyse)[0][survive]
        keep = ~to_lyse
        keep[surviving_idx] = True
        pcn = pcn[keep]
        lysis_timer = lysis_timer[keep]
        n_cells = len(pcn)
        nanobody += N_A * (n_lysing - survive.sum())

# 绘图
fig, axes = plt.subplots(4,1,figsize=(12,10), sharex=True)
axes[0].plot(times, pop_h, color='blue')
axes[0].set_ylabel('Total Viable Cells')
axes[0].set_title('Soft Carrying Capacity – No Hard Ceiling')
axes[0].grid(alpha=0.3)

axes[1].plot(times, pcn_h, color='orange')
axes[1].axhline(PCN_threshold, color='red', ls='--', label='Threshold')
axes[1].set_ylabel('Mean PCN'); axes[1].legend(); axes[1].grid(alpha=0.3)

if lysis_events:
    t_ev, n_ev = zip(*lysis_events)
    axes[2].scatter(t_ev, n_ev, s=20, color='red', alpha=0.7)
axes[2].set_ylabel('Lysis Events'); axes[2].grid(alpha=0.3)

axes[3].plot(times, nb_h, color='green')
axes[3].set_ylabel('Nanobody (a.u.)'); axes[3].set_xlabel('Time (h)')
axes[3].grid(alpha=0.3)

plt.tight_layout()
plt.show()