import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# 设置中文字体，避免乱码
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ==================== 优化参数设置 ====================
# 肿瘤参数
r_T = 0.08              # 肿瘤生长速率
K_T = 3.0               # 肿瘤环境容纳量
k_kill = 0.06           # 最大杀伤系数
nb_EC50 = 2.5e6         # 药物半效应浓度
n_kill = 2.0            # 【优化】Hill系数调高至3.0，突出高浓度峰值的杀伤优势

# 细菌与振荡参数
mu = 0.5                # 细菌最大生长速率
K_bact = 1e9            # 细菌环境容纳量
lysis_thr = 8.0         # 裂解蛋白阈值
survive = 0.25          # 裂解存活率
alpha_E = 2.0           # 裂解蛋白表达速率
delta_E = 0.2           # 裂解蛋白降解速率
burst = 1e6             # 单次裂解释放抗体量

# 药代动力学参数
nb_deg = 0.015          # 【优化】抗体清除率 (半衰期 ~46h)，凸显脉冲清除期的安全性

scale = 1e9
K_bact_norm = K_bact / scale

# ==================== 模型定义 ====================
class PulseOnlyModel:
    def __init__(self, has_bacteria=True):
        self.has_bacteria = has_bacteria

    def model(self, t, y):
        B, E, NB, T = y
        if self.has_bacteria:
            growth = mu * (1 - B / K_bact_norm)
            dE = alpha_E - (delta_E + growth) * E
            dNB = -nb_deg * NB  # 降解
            dB = growth * B
        else:
            dB = 0.0
            dE = 0.0
            dNB = -nb_deg * NB  # 降解
            
        kill = k_kill * NB**n_kill / (NB**n_kill + nb_EC50**n_kill)
        dT = r_T * T * (1 - T / K_T) - kill * T
        return [dB, dE, dNB, dT]

    def lysis_event(self, t, y):
        if not self.has_bacteria:
            return 1.0
        return y[1] - lysis_thr
    
    # 关键：设为True，积分才会在事件处停止并处理状态
    lysis_event.terminal = True 
    lysis_event.direction = 1

    def apply_lysis(self, y):
        B, E, NB, T = y
        NB_new = NB + burst * (1 - survive) * B
        return [survive * B, 0.0, NB_new, T]

def simulate(model, t_span, y0):
    t_current, t_end = t_span
    y_current = np.array(y0, dtype=float)
    ts = [t_current]
    ys = [y_current.copy()]
    max_events = 100
    event_count = 0
    
    while t_current < t_end and event_count < max_events:
        sol = solve_ivp(model.model, [t_current, t_end], y_current,
                        events=model.lysis_event, max_step=0.5, method='LSODA')
        
        if len(sol.t) > 1:
            ts.extend(sol.t[1:])
            ys.extend(sol.y[:, 1:].T)
        
        t_current = sol.t[-1]
        y_current = sol.y[:, -1]
        
        if sol.t_events[0].size > 0:
            t_event = sol.t_events[0][0]
            y_event = sol.y_events[0][0]
            t_current = t_event
            y_current = y_event
            y_current = model.apply_lysis(y_current)
            ts.append(t_current)
            ys.append(y_current)
            event_count += 1
        else:
            break
    return np.array(ts), np.array(ys)

# ==================== 运行对比 ====================
def run_comparison():
    print("=== 运行头对头比较：脉冲 vs 持续 ===")
    
    # 1. 模拟脉冲组 (获取总释放量)
    y0_p = [0.08, 0.0, 0.0, 1.0]
    ts_p, ys_p = simulate(PulseOnlyModel(True), [0, 200], y0_p)
    total_dose_pulsed = ys_p[-1, 2]  # 脉冲组最终累积量
    
    # 2. 计算持续给药速率 (确保总剂量相等)
    R_infusion = total_dose_pulsed / 200.0
    
    def model_continuous(t, y):
        NB, T = y
        dNB = R_infusion - nb_deg * NB
        kill = k_kill * NB**n_kill / (NB**n_kill + nb_EC50**n_kill)
        dT = r_T * T * (1 - T/K_T) - kill * T
        return [dNB, dT]
        
    sol_c = solve_ivp(model_continuous, [0, 200], [0.0, 1.0], 
                      method='LSODA', max_step=0.5)
    ts_c, ys_c = sol_c.t, sol_c.y.T
    
    # 3. 数据统计
    tumor_sup_p = 1 - ys_p[-1, 3] / 1.0
    tumor_sup_c = 1 - ys_c[-1, 1] / 1.0
    
    print(f"\n 总结（等总剂量）:")
    print(f"脉冲   | 肿瘤抑制率: {tumor_sup_p:.1%} | 峰值浓度: {np.max(ys_p[:,2]):.2e}")
    print(f"持续   | 肿瘤抑制率: {tumor_sup_c:.1%} | 峰值浓度: {np.max(ys_c[:,0]):.2e}")
    print(f"优势   | 抑制率差值: +{(tumor_sup_p - tumor_sup_c):.1%}")

    # 4. 绘图
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # (1) 肿瘤生长曲线
    axes[0,0].plot(ts_c, ys_c[:,1], 'g-', linewidth=2.5, label='持续给药')
    axes[0,0].plot(ts_p, ys_p[:,3], 'r-', linewidth=2, label='脉冲释放 (ChronoBiotics)')
    axes[0,0].axhline(1.0, color='gray', ls='--', alpha=0.5)
    axes[0,0].set_ylabel('肿瘤大小 (相对值)', fontsize=11)
    axes[0,0].set_title('肿瘤生长抑制', fontsize=12, fontweight='bold')
    axes[0,0].legend(); axes[0,0].grid(alpha=0.3)
    
    # (2) 药物暴露曲线
    axes[0,1].plot(ts_c, ys_c[:,0], 'g-', linewidth=2.5, label='持续')
    axes[0,1].plot(ts_p, ys_p[:,2], 'r-', linewidth=1.5, label='脉冲')
    axes[0,1].axhline(nb_EC50, color='purple', ls=':', label=f'EC50 ({nb_EC50:.0e})')
    axes[0,1].set_ylabel('纳米抗体浓度 (任意单位)', fontsize=11)
    axes[0,1].set_title('药物暴露曲线', fontsize=12, fontweight='bold')
    axes[0,1].legend(); axes[0,1].grid(alpha=0.3)
    
    # (3) 药效输出 (Hill)
    kill_p = k_kill * ys_p[:,2]**n_kill / (ys_p[:,2]**n_kill + nb_EC50**n_kill)
    kill_c = k_kill * ys_c[:,0]**n_kill / (ys_c[:,0]**n_kill + nb_EC50**n_kill)
    axes[1,0].plot(ts_c, kill_c, 'g-', linewidth=2.5, label='持续杀伤')
    axes[1,0].plot(ts_p, kill_p, 'r-', linewidth=1.5, label='脉冲杀伤')
    axes[1,0].set_ylabel('瞬时杀伤速率 (1/h)', fontsize=11)
    axes[1,0].set_title('药效输出 (Hill)', fontsize=12, fontweight='bold')
    axes[1,0].legend(); axes[1,0].grid(alpha=0.3)
    
    # (4) 关键指标对比 (柱状图)
    metrics = ['肿瘤\n抑制率', '峰值\n浓度', '平均\n浓度']
    
    # 原始数值
    val_sup_p = tumor_sup_p * 100
    val_sup_c = tumor_sup_c * 100
    val_peak_p = np.max(ys_p[:, 2])
    val_peak_c = np.max(ys_c[:, 0])
    val_avg_p = np.mean(ys_p[:, 2])
    val_avg_c = np.mean(ys_c[:, 0])
    
    # 归一化：每个指标以两组中的最大值为基准，便于可视化比较 (0-1 范围)
    sup_max = max(val_sup_p, val_sup_c)
    peak_max = max(val_peak_p, val_peak_c)
    avg_max = max(val_avg_p, val_avg_c)
    
    vals_p = [val_sup_p/sup_max, val_peak_p/peak_max, val_avg_p/avg_max]
    vals_c = [val_sup_c/sup_max, val_peak_c/peak_max, val_avg_c/avg_max]
    
    x_pos = np.arange(len(metrics))
    width = 0.35
    
    axes[1,1].bar(x_pos - width/2, vals_p, width, label='脉冲', color='#FF4B4B', edgecolor='black')
    axes[1,1].bar(x_pos + width/2, vals_c, width, label='持续', color='#4CAF50', edgecolor='black', alpha=0.7)
    axes[1,1].set_xticks(x_pos)
    axes[1,1].set_xticklabels(metrics, fontsize=10)
    axes[1,1].set_ylabel('相对值 (归一化)', fontsize=11)
    axes[1,1].set_title('关键指标对比', fontsize=12, fontweight='bold')
    axes[1,1].legend(); axes[1,1].grid(axis='y', alpha=0.3)
    axes[1,1].set_ylim(0, 1.15)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_comparison()