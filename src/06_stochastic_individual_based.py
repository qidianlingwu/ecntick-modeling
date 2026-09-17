"""
随机个体模型（stochastic individual-based model, IBM）

【方法说明 —— 重要】
本脚本**不是 Gillespie 算法（SSA）**。原文件名为 Gillespie.py，属命名不当，已更正。
实际实现的是：固定时间步长 dt 的随机个体模型，每个细胞的分裂、质粒复制/丢失、
裂解存活等事件用二项/泊松分布逐步采样。它介于确定性 ODE 与精确 SSA 之间，
接近 tau-leaping 的思路。

与 Gillespie SSA 的区别：
  - SSA：由反应倾向（propensity）之和抽取指数分布的等待时间，事件逐个发生，时间步不等距；
  - 本脚本：时间步固定为 dt，每个 dt 内对所有细胞并行采样事件。

【建模对象说明】
本脚本建模的是**早期 QS（群体感应）依赖的设计**（含 AHL 信号 S 及其阈值 S_threshold），
该方案在项目后期被"无 QS 依赖的弛豫振荡"方案取代。保留此脚本的目的是：
用随机模拟检验群体层面的振荡在有限细胞数下是否会被噪声抹掉。
对应最终设计（弛豫振荡）的随机版本尚未实现，见 LIMITATIONS.md。

输出：gillespie_output.png（文件名沿用历史，未改）
"""
import numpy as np
import matplotlib.pyplot as plt

# ===== 参数 =====
mu_max = 0.5            # 最大生长速率 (/h)
r_plasmid = 1.2         # 质粒复制速率 (/h)
d_plasmid = 0.1         # 质粒丢失速率 (/h)
PCN_threshold = 30      # 裂解所需的最低 PCN
f_survive = 0.15        # 裂解后存活率
lysis_delay = 0.3       # 裂解蛋白积累延迟 (h)
N_A = 1000              # 每裂解释放纳米抗体 (a.u.)
MAX_PCN = 500           # 拷贝数硬上限

# 群体感应参数
S_prod = 0.08           # 每细胞信号分子生产率 (/h)
S_deg = 1.2             # 信号分子降解速率 (/h)
S_threshold = 60        # 裂解能力激活阈值

n_init = 200            # 初始细胞数（模拟低密度接种）
pcns_init = np.random.poisson(15, n_init).astype(np.int64)
dt = 0.01
T_end = 80
steps = int(T_end / dt)

times = np.linspace(0, T_end, steps + 1)
pop_history = np.zeros(steps + 1)
mean_pcn_history = np.zeros(steps + 1)
lysis_events = []
nanobody = 0.0
nanobody_history = np.zeros(steps + 1)
signal_history = np.zeros(steps + 1)

pcn = pcns_init.copy()
n_cells = n_init
lysis_timer = np.zeros(n_cells)
S = 0.0  # 信号分子初始浓度

for i, t in enumerate(times):
    pop_history[i] = n_cells
    mean_pcn_history[i] = np.mean(pcn) if n_cells > 0 else 0
    nanobody_history[i] = nanobody
    signal_history[i] = S

    if n_cells == 0:
        break

    # --- 1. 细胞分裂 (仅受PCN负担抑制，无容量限制) ---
    div_prob = mu_max * dt / (1 + pcn / PCN_threshold)
    div_prob = np.clip(div_prob, 0, 1)
    will_divide = np.random.binomial(1, div_prob)
    n_div = will_divide.sum()
    if n_div > 0:
        dividing_idx = np.where(will_divide)[0]
        new_pcn = np.zeros(n_div, dtype=np.int64)
        new_timers = np.zeros(n_div)
        for j, idx in enumerate(dividing_idx):
            mother = pcn[idx]
            d1 = np.random.binomial(mother, 0.5) if mother > 0 else 0
            d2 = mother - d1
            pcn[idx] = d1
            new_pcn[j] = d2
            lysis_timer[idx] = 0.0  # 分裂稀释裂解蛋白
        pcn = np.concatenate([pcn, new_pcn])
        lysis_timer = np.concatenate([lysis_timer, new_timers])
        n_cells = len(pcn)

    # --- 2. 质粒复制与丢失 ---
    if n_cells > 0:
        rep = np.random.binomial(pcn.astype(np.int64), r_plasmid * dt)
        pcn = pcn + rep
        loss = np.random.binomial(pcn.astype(np.int64), d_plasmid * dt)
        pcn = pcn - loss
        pcn = np.clip(pcn, 1, MAX_PCN)

    # --- 3. 群体感应信号 ---
    # 信号分子由所有细胞持续产生，按一级动力学降解
    S += S_prod * n_cells * dt
    S -= S_deg * S * dt
    if S < 0:
        S = 0.0

    # --- 4. 裂解触发 (需同时满足: PCN >= 阈值 AND S >= 阈值) ---
    if n_cells > 0:
        # 只有信号浓度够高时，高PCN细胞才开始积累裂解蛋白
        lysis_competent = (pcn >= PCN_threshold) & (S >= S_threshold)
        lysis_timer[lysis_competent] += dt
        lysis_timer[~lysis_competent] = 0
        to_lyse = lysis_timer >= lysis_delay

        if np.any(to_lyse):
            n_lysing = np.sum(to_lyse)
            lysis_events.append((t, n_lysing))
            survive = np.random.rand(n_lysing) < f_survive
            lysing_idx = np.where(to_lyse)[0]
            surviving_idx = lysing_idx[survive]
            n_dead = n_lysing - survive.sum()

            lysis_timer[surviving_idx] = 0.0
            pcn[surviving_idx] = np.random.poisson(5, len(surviving_idx)).astype(np.int64)
            keep = ~to_lyse
            keep[surviving_idx] = True
            pcn = pcn[keep]
            lysis_timer = lysis_timer[keep]
            n_cells = len(pcn)
            nanobody += N_A * n_dead

# ===== 绘图 =====
fig, axes = plt.subplots(5, 1, figsize=(12, 12), sharex=True)

axes[0].plot(times[:i + 1], pop_history[:i + 1], color='blue')
axes[0].set_ylabel('Total Viable Cells')
axes[0].set_title('Quorum-Sensing Model: Plasmid-Driven Oscillations')
axes[0].grid(alpha=0.3)

axes[1].plot(times[:i + 1], signal_history[:i + 1], color='purple')
axes[1].axhline(S_threshold, color='red', ls='--', label='QS Threshold')
axes[1].set_ylabel('Signal (AHL)')
axes[1].legend()
axes[1].grid(alpha=0.3)

axes[2].plot(times[:i + 1], mean_pcn_history[:i + 1], color='orange')
axes[2].axhline(PCN_threshold, color='red', ls='--', label='PCN Threshold')
axes[2].set_ylabel('Mean PCN')
axes[2].legend()
axes[2].grid(alpha=0.3)

if lysis_events:
    t_ev, n_ev = zip(*lysis_events)
    axes[3].scatter(t_ev, n_ev, s=20, color='red', alpha=0.7)
axes[3].set_ylabel('Lysis Events')
axes[3].grid(alpha=0.3)

axes[4].plot(times[:i + 1], nanobody_history[:i + 1], color='green')
axes[4].set_ylabel('Nanobody (a.u.)')
axes[4].set_xlabel('Time (h)')
axes[4].grid(alpha=0.3)

plt.tight_layout()
plt.savefig('gillespie_output.png', dpi=150, bbox_inches='tight')
plt.show()

print(f"Final population: {pop_history[i]:.0f}")
print(f"Population range: {pop_history[:i+1].min():.0f} – {pop_history[:i+1].max():.0f}")
print(f"Lysis events: {len(lysis_events)}")
print(f"Final nanobody: {nanobody_history[i]:.1f}")
print(f"Signal range: {signal_history[:i+1].min():.1f} – {signal_history[:i+1].max():.1f}")
if lysis_events:
    t_ev, n_ev = zip(*lysis_events)
    print(f"Lysis timing: first at t={t_ev[0]:.1f}h, last at t={t_ev[-1]:.1f}h")
