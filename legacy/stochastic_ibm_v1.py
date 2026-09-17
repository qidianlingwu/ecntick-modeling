import numpy as np
import matplotlib.pyplot as plt

# ===== 参数（经过结构平衡） =====
mu_max = 0.5            # 最大生长速率 (/h)
K_bact = 1000
           # 环境承载容量（保证计算效率）
r_plasmid = 1.2        # 质粒复制速率
d_plasmid = 0.1         # 质粒丢失速率
PCN_threshold = 30      # 裂解阈值（pBBR1稳态~40，低于稳态）
f_survive = 0.15        # 裂解后存活率
lysis_delay = 0.3       # 裂解蛋白积累延迟 (h)
N_A = 1000              # 每裂解释放纳米抗体
MAX_PCN = 500           # 拷贝数硬上限

n_init = 300            # 初始细胞数
pcns_init = np.random.poisson(15, n_init).astype(np.int64)  # 初始PCN约15
dt = 0.02               # 模拟步长 (h)
T_end = 80              # 模拟时长 (h)
steps = int(T_end / dt)

# 记录数组
times = np.linspace(0, T_end, steps+1)
pop_history = np.zeros(steps+1)
mean_pcn_history = np.zeros(steps+1)
lysis_events = []
nanobody = 0.0
nanobody_history = np.zeros(steps+1)

# 初始种群
pcn = pcns_init.copy()
n_cells = n_init
lysis_timer = np.zeros(n_cells)

# 主模拟循环
for i, t in enumerate(times):
    # 保存当前状态
    pop_history[i] = n_cells
    mean_pcn_history[i] = np.mean(pcn) if n_cells > 0 else 0
    nanobody_history[i] = nanobody

    if n_cells == 0:
        break

    # --- 1. 细胞分裂 (只有存活细胞且未达容量时才分裂) ---
    if n_cells < K_bact:
        # 分裂概率受拷贝数抑制（高拷贝负担大）
        div_prob = mu_max * dt / (1 + pcn / PCN_threshold)
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
            lysis_timer = np.zeros(n_cells)  # 新生细胞无积累蛋白
    else:
        # 超过容量限制，不再分裂，避免计算爆炸
        pass

    # --- 2. 质粒复制与丢失 ---
    if n_cells > 0:
        rep = np.random.binomial(pcn.astype(np.int64), r_plasmid * dt)
        pcn = pcn + rep
        loss = np.random.binomial(pcn.astype(np.int64), d_plasmid * dt)
        pcn = pcn - loss
        pcn = np.clip(pcn, 1, MAX_PCN)

    # --- 3. 裂解触发 (拷贝数高于阈值并持续一段时间) ---
    in_risk = pcn >= PCN_threshold
    lysis_timer[in_risk] += dt
    lysis_timer[~in_risk] = 0
    to_lyse = lysis_timer >= lysis_delay

    if np.any(to_lyse):
        n_lysing = np.sum(to_lyse)
        lysis_events.append((t, n_lysing))
        # 每个裂解细胞独立存活
        survive = np.random.rand(n_lysing) < f_survive
        surviving_idx = np.where(to_lyse)[0][survive]
        keep = ~to_lyse
        keep[surviving_idx] = True
        pcn = pcn[keep]
        lysis_timer = lysis_timer[keep]
        n_cells = len(pcn)
        nanobody += N_A * (n_lysing - survive.sum())

    # (不再需要额外的容量裁剪，分裂时已控制)

# ===== 绘图 =====
fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)
axes[0].plot(times[:i+1], pop_history[:i+1], color='blue')
axes[0].set_ylabel('Total Viable Cells')
axes[0].set_title('Robust Stochastic Oscillations (Gillespie)')
axes[0].grid(alpha=0.3)

axes[1].plot(times[:i+1], mean_pcn_history[:i+1], color='orange')
axes[1].axhline(PCN_threshold, color='red', ls='--', label='Lysis Threshold')
axes[1].set_ylabel('Mean Plasmid Copy Number')
axes[1].legend(); axes[1].grid(alpha=0.3)

if lysis_events:
    t_ev, n_ev = zip(*lysis_events)
    axes[2].scatter(t_ev, n_ev, s=20, color='red', alpha=0.7)
axes[2].set_ylabel('Lysis Events')
axes[2].grid(alpha=0.3)

axes[3].plot(times[:i+1], nanobody_history[:i+1], color='green')
axes[3].set_ylabel('Released Nanobody (a.u.)')
axes[3].set_xlabel('Time (h)')
axes[3].grid(alpha=0.3)

plt.tight_layout()
plt.show()