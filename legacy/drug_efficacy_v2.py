import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# ==================== 参数 ====================
r_T, K_T = 0.08, 3.0
A_0, k_rec, k_ex, K_ex, n_ex = 0.8, 0.03, 0.04, 1.5e6, 4
k_kill, nb_EC50, n_kill, A_threshold = 0.15, 2.0e6, 2.0, 0.3
mu, K_bact = 0.5, 1e9
lysis_thr, survive, alpha_E, delta_E, burst, nb_deg = 8.0, 0.25, 2.0, 0.2, 1.5e6, 0.015
scale, K_bact_norm = 1e9, 1e9 / 1e9

# ==================== 脉冲模型 (状态: B, E, NB, T, NB_total, A) ====================
class PulsedModel:
    def model(self, t, y):
        B, E, NB, T, NB_total, A = y
        growth = mu * (1 - B / K_bact_norm)
        dE = alpha_E - (delta_E + growth) * E
        dNB, dNB_total, dB = -nb_deg * NB, 0.0, growth * B
        drug_kill = k_kill * NB**n_kill / (NB**n_kill + nb_EC50**n_kill)
        kill = drug_kill * (A / (A + A_threshold)) * T
        dT = r_T * T * (1 - T / K_T) - kill
        exhaustion = k_ex * NB**n_ex / (NB**n_ex + K_ex**n_ex) * A
        dA = -exhaustion + k_rec * (1 - A)
        return [dB, dE, dNB, dT, dNB_total, dA]
    def lysis_event(self, t, y): return y[1] - lysis_thr
    lysis_event.terminal, lysis_event.direction = True, 1
    def apply_lysis(self, y):
        B, E, NB, T, NB_total, A = y
        released = burst * (1 - survive) * B
        return [survive * B, 0.0, NB + released, T, NB_total + released, A]

def simulate_pulsed(t_span, y0):
    t, y, ts, ys, events = t_span[0], np.array(y0), [t_span[0]], [np.array(y0)], 0
    while t < t_span[1] and events < 100:
        sol = solve_ivp(PulsedModel().model, [t, t_span[1]], y,
                        events=PulsedModel().lysis_event, max_step=0.5)
        if len(sol.t) > 1:
            ts.extend(sol.t[1:])
            ys.extend(sol.y[:, 1:].T)
        t, y = sol.t[-1], sol.y[:, -1]
        if sol.t_events[0].size > 0:
            t = sol.t_events[0][0]
            y = PulsedModel().apply_lysis(sol.y_events[0][0])
            ts.append(t)
            ys.append(y)
            events += 1
        else:
            break
    return np.array(ts), np.array(ys)

# ==================== 持续给药模型 (状态: NB, T, A) ====================
class ContinuousModel:
    def __init__(self, R): self.R = R
    def model(self, t, y):
        NB, T, A = y
        drug_kill = k_kill * NB**n_kill / (NB**n_kill + nb_EC50**n_kill)
        kill = drug_kill * (A / (A + A_threshold)) * T
        dT = r_T * T * (1 - T / K_T) - kill
        exhaustion = k_ex * NB**n_ex / (NB**n_ex + K_ex**n_ex) * A
        dA = -exhaustion + k_rec * (1 - A)
        return [self.R - nb_deg * NB, dT, dA]

# ==================== 运行 ====================
def run():
    print("=== Pulsed vs Continuous (T-cell exhaustion model) ===")
    # 1. 脉冲组模拟
    y0_p = [0.08, 0.0, 0.0, 1.0, 0.0, A_0]
    ts_p, ys_p = simulate_pulsed([0, 200], y0_p)
    total_dose = ys_p[-1, 4]
    R_inf = total_dose / 200.0
    # 2. 持续组模拟
    sol_c = solve_ivp(ContinuousModel(R_inf).model, [0, 200], [0.0, 1.0, A_0],
                      method='LSODA', max_step=0.5)
    ts_c, ys_c = sol_c.t, sol_c.y
    # 3. 提取数据 (从统一索引)
    T_p, A_p = ys_p[:, 3], ys_p[:, 5]   # 脉冲组肿瘤在索引3，T细胞活性在索引5
    T_c, A_c = ys_c[1], ys_c[2]         # 持续组肿瘤在索引1，T细胞活性在索引2
    NB_p, NB_c = ys_p[:, 2], ys_c[0]    # 药物浓度
    sup_p = 1 - T_p[-1]; sup_c = 1 - T_c[-1]
    print(f"总剂量: {total_dose:.2e}")
    print(f"脉冲  抑制率: {sup_p:.1%}  | 终末T活性: {A_p[-1]:.3f}")
    print(f"持续  抑制率: {sup_c:.1%}  | 终末T活性: {A_c[-1]:.3f}")
    print(f"优势: +{(sup_p - sup_c):.1%}")

    # 4. 绘图 (所有数据来自同一来源)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes[0,0].plot(ts_c, T_c, 'g-', lw=2.5, label='持续给药')
    axes[0,0].plot(ts_p, T_p, 'r-', lw=2, label='脉冲释放')
    axes[0,0].axhline(1.0, color='gray', ls='--')
    axes[0,0].set_ylabel('肿瘤大小'); axes[0,0].set_title('肿瘤生长抑制')
    axes[0,0].legend(); axes[0,0].grid(alpha=0.3)
    
    axes[0,1].plot(ts_c, A_c, 'g-', lw=2.5, label='持续')
    axes[0,1].plot(ts_p, A_p, 'r-', lw=2, label='脉冲')
    axes[0,1].set_ylabel('T细胞活性'); axes[0,1].set_title('T细胞活性动态')
    axes[0,1].legend(); axes[0,1].grid(alpha=0.3)
    
    axes[1,0].plot(ts_c, NB_c, 'g-', lw=2.5, label='持续')
    axes[1,0].plot(ts_p, NB_p, 'r-', lw=1.5, label='脉冲')
    axes[1,0].axhline(nb_EC50, color='purple', ls=':', label='EC50')
    axes[1,0].set_ylabel('药物浓度'); axes[1,0].set_title('药物暴露曲线')
    axes[1,0].legend(); axes[1,0].grid(alpha=0.3)
    
    metrics = ['抑制率', '终末T活性', '平均药物']
    val_p = [sup_p*100, A_p[-1], np.mean(NB_p)]
    val_c = [sup_c*100, A_c[-1], np.mean(NB_c)]
    max_v = [max(a,b) for a,b in zip(val_p, val_c)]
    norm_p = [v/m for v,m in zip(val_p, max_v)]
    norm_c = [v/m for v,m in zip(val_c, max_v)]
    x = np.arange(len(metrics))
    axes[1,1].bar(x-0.2, norm_p, 0.4, label='脉冲', color='#FF4B4B')
    axes[1,1].bar(x+0.2, norm_c, 0.4, label='持续', color='#4CAF50')
    axes[1,1].set_xticks(x); axes[1,1].set_xticklabels(metrics)
    axes[1,1].set_ylabel('相对值'); axes[1,1].set_title('综合对比')
    axes[1,1].legend(); axes[1,1].grid(axis='y', alpha=0.3)
    axes[1,1].set_ylim(0, 1.2)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run()