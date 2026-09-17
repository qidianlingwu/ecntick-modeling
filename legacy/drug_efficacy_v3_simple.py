import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# ==================== 参数 ====================
r_T = 0.08
K_T = 3.0
k_kill = 0.05
nb_EC50 = 3.6e6
n_kill = 2.0

mu = 0.5
K_bact = 1e9
lysis_thr = 8.0
survive = 0.25
alpha_E = 2.0
delta_E = 0.2
burst = 1.0e6
nb_deg = 0.015

scale = 1e9
K_bact_norm = K_bact / scale

# ==================== 脉冲模型 ====================
class PulsedModel:
    def model(self, t, y):
        B, E, NB, T, NB_total = y
        growth = mu * (1 - B / K_bact_norm)
        dE = alpha_E - (delta_E + growth) * E
        dNB = -nb_deg * NB
        dNB_total = 0.0
        dB = growth * B
        kill = k_kill * NB**n_kill / (NB**n_kill + nb_EC50**n_kill)
        dT = r_T * T * (1 - T / K_T) - kill * T
        return [dB, dE, dNB, dT, dNB_total]

    def lysis_event(self, t, y):
        return y[1] - lysis_thr

    def apply_lysis(self, y):
        B, E, NB, T, NB_total = y
        released = burst * (1 - survive) * B
        return [survive * B, 0.0, NB + released, T, NB_total + released]

def simulate_pulsed(t_span, y0):
    pm = PulsedModel()

    # 事件包装函数
    def event_func(t, y):
        return pm.lysis_event(t, y)
    event_func.terminal = True
    event_func.direction = 1

    t, y = t_span[0], np.array(y0)
    ts, ys = [t], [y.copy()]
    events = 0
    while t < t_span[1] and events < 100:
        sol = solve_ivp(pm.model, [t, t_span[1]], y, events=event_func, max_step=0.5)
        if len(sol.t) > 1:
            ts.extend(sol.t[1:])
            ys.extend(sol.y[:, 1:].T)
        t = sol.t[-1]
        y = sol.y[:, -1]
        if sol.t_events[0].size > 0:
            t = sol.t_events[0][0]
            y = pm.apply_lysis(sol.y_events[0][0])
            ts.append(t)
            ys.append(y)
            events += 1
        else:
            break
    return np.array(ts), np.array(ys)

# ==================== 持续给药模型 ====================
class ContinuousModel:
    def __init__(self, R):
        self.R = R
    def model(self, t, y):
        NB, T = y
        dNB = self.R - nb_deg * NB
        kill = k_kill * NB**n_kill / (NB**n_kill + nb_EC50**n_kill)
        dT = r_T * T * (1 - T / K_T) - kill * T
        return [dNB, dT]

# ==================== 运行并绘图 ====================
def run_and_plot():
    y0_p = [0.08, 0.0, 0.0, 1.0, 0.0]
    ts_p, ys_p = simulate_pulsed([0, 200], y0_p)
    total_dose = ys_p[-1, 4]
    R_inf = total_dose / 200.0

    sol_c = solve_ivp(ContinuousModel(R_inf).model, [0, 200], [0.0, 1.0], method='LSODA', max_step=0.5)
    ts_c, ys_c = sol_c.t, sol_c.y.T

    T_p_final = max(ys_p[-1, 3], 0.0)
    T_c_final = max(ys_c[-1, 1], 0.0)
    sup_p = 1 - T_p_final / 1.0
    sup_c = 1 - T_c_final / 1.0

    # 绘图
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    axes[0,0].plot(ts_c, ys_c[:,1], 'g-', lw=2.5, label='持续给药')
    axes[0,0].plot(ts_p, ys_p[:,3], 'r-', lw=2, label='脉冲释放')
    axes[0,0].axhline(1.0, color='gray', ls='--', alpha=0.5)
    axes[0,0].set_ylabel('肿瘤大小')
    axes[0,0].set_title('肿瘤生长抑制')
    axes[0,0].legend()
    axes[0,0].grid(alpha=0.3)

    axes[0,1].plot(ts_c, ys_c[:,0], 'g-', lw=2.5, label='持续')
    axes[0,1].plot(ts_p, ys_p[:,2], 'r-', lw=1.5, label='脉冲')
    axes[0,1].axhline(nb_EC50, color='purple', ls=':', label='EC50')
    axes[0,1].set_ylabel('药物浓度')
    axes[0,1].set_title('药物暴露曲线')
    axes[0,1].legend()
    axes[0,1].grid(alpha=0.3)

    kill_p = k_kill * ys_p[:,2]**n_kill / (ys_p[:,2]**n_kill + nb_EC50**n_kill)
    kill_c = k_kill * ys_c[:,0]**n_kill / (ys_c[:,0]**n_kill + nb_EC50**n_kill)
    axes[1,0].plot(ts_c, kill_c, 'g-', lw=2.5, label='持续杀伤')
    axes[1,0].plot(ts_p, kill_p, 'r-', lw=1.5, label='脉冲杀伤')
    axes[1,0].set_ylabel('瞬时杀伤速率')
    axes[1,0].set_title('药效输出')
    axes[1,0].legend()
    axes[1,0].grid(alpha=0.3)

    # 关键指标对比：同一指标内归一化到两者中的较大值，不做任何人为缩放
    # 注意 sup_p/sup_c 是相对“初始肿瘤大小”的回缩率，不是相对无治疗对照的 TGI，
    # 因此指标名使用“肿瘤回缩率”，避免与 TGI 混淆。
    metrics = ['肿瘤回缩率', '峰值浓度', '平均浓度']
    val_sup_p, val_sup_c = sup_p, sup_c
    val_peak_p = np.max(ys_p[:, 2])
    val_peak_c = np.max(ys_c[:, 0])
    # 事件驱动模型的时间网格是非等间距的，必须用时间加权均值，不能用 np.mean
    val_avg_p = np.trapezoid(ys_p[:, 2], ts_p) / (ts_p[-1] - ts_p[0])
    val_avg_c = np.trapezoid(ys_c[:, 0], ts_c) / (ts_c[-1] - ts_c[0])

    vals_p, vals_c = [], []
    for a, b in zip([val_sup_p, val_peak_p, val_avg_p],
                    [val_sup_c, val_peak_c, val_avg_c]):
        m = max(abs(a), abs(b), 1e-12)
        vals_p.append(a / m)
        vals_c.append(b / m)

    x_pos = np.arange(len(metrics))
    width = 0.35
    axes[1,1].bar(x_pos - width/2, vals_p, width, label='脉冲', color='#FF4B4B')
    axes[1,1].bar(x_pos + width/2, vals_c, width, label='持续', color='#4CAF50')
    axes[1,1].set_xticks(x_pos)
    axes[1,1].set_xticklabels(metrics)
    axes[1,1].set_ylabel('相对值 (归一化)')
    axes[1,1].set_title('关键指标对比')
    axes[1,1].legend()
    axes[1,1].grid(axis='y', alpha=0.3)
    axes[1,1].set_ylim(0, 1.0)

    plt.tight_layout()
    plt.show()

    print(f"脉冲抑制率: {sup_p:.3f}, 持续抑制率: {sup_c:.3f}")

# 运行
run_and_plot()