import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# =========================================================
# EcNTick 新版 PK/PD-免疫节律模型
# 核心改动：
# 1. 药物浓度 NB 不再是唯一决定因素；
# 2. 加入 T 细胞活性 A；
# 3. 加入持续暴露负担 X，表示长期 PD-L1 阻断/持续刺激带来的免疫疲劳；
# 4. 加入裂解佐剂信号 J，表示细菌裂解后释放 PAMPs/菌体碎片带来的局部免疫激活；
# 5. 持续组与脉冲组仍保持等总药物剂量比较。
# =========================================================

# ==================== 全局参数 ====================
t_end = 200

# 肿瘤生长参数
r_T = 0.08
K_T = 3.0

# 药效参数
# 取值依据：k_kill 采用文献报道抗 PD-L1 阻断效应的保守下界，以避免高估疗效；
# nb_EC50 为量级估计（非实验实测值）。二者对结论的影响已通过参数扫描检验，
# 见 figs/fig5_sensitivity_analysis.png，逐参数来源见根目录 PARAMETERS.md。
k_kill = 0.095
nb_EC50 = 2.0e6
n_kill = 2.0

# T 细胞活性参数
A_0 = 0.8
A_threshold = 0.3
k_rec = 0.055
k_ex = 0.16

# 持续暴露负担 X
k_x_on = 0.45
k_x_off = 0.20

# 裂解佐剂信号 J
# eta_J 取文献报道佐剂效应区间的保守下界，避免高估裂解带来的免疫增益。
# 结论对 eta_J 的敏感性见 figs/fig5_sensitivity_analysis.png。详见 PARAMETERS.md。
J_burst = 0.65
J_deg = 0.16
K_J = 0.20
eta_J = 0.45
k_J_A = 0.045
# 脉冲裂解带来的局部免疫激活也可促进疲劳状态的消退，避免 X 被解释成单纯毒性累积。
k_J_clear = 0.12

# 细菌与裂解参数
mu = 0.5
K_bact = 1e9
scale = 1e9
K_bact_norm = K_bact / scale
# 注意：alpha_E / delta_E 必须高于 lysis_thr，否则 E 永远到不了裂解阈值，模型就不会释放药物。
# 这里采用保守但可触发裂解的组合：稳态 E≈9.4，高于阈值 8.0。
lysis_thr = 8.0
survive = 0.25
alpha_E = 1.70
delta_E = 0.18
burst = 1.15e6

# 裂解释放后的组织扩散/释放平滑项
# 裂解时先进入局部“释放库” Depot，再逐步转为游离纳米抗体 NB，避免图上出现过尖的垂直锯齿峰。
k_abs = 0.16

# 药物清除参数
nb_deg = 0.03

# ==================== 工具函数 ====================
def drug_effect(NB):
    """药物占有效应，0-1。"""
    return NB**n_kill / (NB**n_kill + nb_EC50**n_kill + 1e-30)


def adjuvant_effect(J):
    """裂解佐剂信号，0-1。"""
    return J / (J + K_J + 1e-30)


def time_weighted_mean(t, y):
    """非等间距采样下的时间加权均值。"""
    return np.trapezoid(y, t) / (t[-1] - t[0])

# ==================== 脉冲治疗模型 ====================
class PulsedImmuneModel:
    def model(self, t, y):
        # y = [B, E, Depot, NB, T, A, X, J, NB_total]
        B, E, Depot, NB, T, A, X, J, NB_total = y

        growth = mu * (1 - B / K_bact_norm)
        dB = growth * B
        dE = alpha_E - (delta_E + growth) * E
        dDepot = -k_abs * Depot
        dNB = k_abs * Depot - nb_deg * NB
        dNB_total = 0.0

        eff = drug_effect(NB)
        j_eff = adjuvant_effect(J)

        # A 表示 T 细胞活性；A 越低，免疫杀伤越弱
        A_factor = A / (A + A_threshold)

        # 脉冲组额外拥有裂解佐剂效应 J
        immune_boost = 1 + eta_J * j_eff
        kill = k_kill * eff * A_factor * immune_boost

        dT = r_T * T * (1 - T / K_T) - kill * T

        # X 表示持续暴露负担：药物作用越持续，X 越累积；低谷期 X 衰减。
        # 脉冲组额外存在裂解佐剂信号，可促进疲劳/抑制状态的消退，因此加入 j_eff 相关清除项。
        dX = k_x_on * eff - (k_x_off + k_J_clear * j_eff) * X

        # A 表示 T 细胞活性：自然恢复 - 持续暴露耗竭 + 裂解佐剂激活
        dA = k_rec * (1 - A) - k_ex * X * A + k_J_A * j_eff * (1 - A)

        # 裂解佐剂信号衰减
        dJ = -J_deg * J

        return [dB, dE, dDepot, dNB, dT, dA, dX, dJ, dNB_total]

    def lysis_event(self, t, y):
        return y[1] - lysis_thr

    def apply_lysis(self, y):
        B, E, Depot, NB, T, A, X, J, NB_total = y
        released_NB = burst * (1 - survive) * B
        released_J = J_burst * (1 - survive) * B
        return [
            survive * B,
            0.0,
            Depot + released_NB,
            NB,
            T,
            A,
            X,
            J + released_J,
            NB_total + released_NB,
        ]


def simulate_pulsed(t_span, y0):
    model = PulsedImmuneModel()

    def event_func(t, y):
        return model.lysis_event(t, y)

    event_func.terminal = True
    event_func.direction = 1

    t_current, t_end_local = t_span
    y_current = np.array(y0, dtype=float)
    ts = [t_current]
    ys = [y_current.copy()]

    event_count = 0
    max_events = 200

    while t_current < t_end_local and event_count < max_events:
        sol = solve_ivp(
            model.model,
            [t_current, t_end_local],
            y_current,
            events=event_func,
            max_step=0.25,
            method='LSODA'
        )

        if len(sol.t) > 1:
            ts.extend(sol.t[1:])
            ys.extend(sol.y[:, 1:].T)

        t_current = sol.t[-1]
        y_current = sol.y[:, -1]

        if sol.t_events[0].size > 0:
            t_event = sol.t_events[0][0]
            y_event = sol.y_events[0][0]
            y_current = model.apply_lysis(y_event)
            t_current = t_event
            ts.append(t_current)
            ys.append(y_current)
            event_count += 1
        else:
            break

    return np.array(ts), np.array(ys)

# ==================== 持续给药模型 ====================
class ContinuousImmuneModel:
    def __init__(self, R_infusion):
        self.R = R_infusion

    def model(self, t, y):
        # y = [NB, T, A, X]
        NB, T, A, X = y

        dNB = self.R - nb_deg * NB
        eff = drug_effect(NB)
        A_factor = A / (A + A_threshold)

        # 持续给药没有细菌周期性裂解释放，因此不含 J 佐剂项
        kill = k_kill * eff * A_factor
        dT = r_T * T * (1 - T / K_T) - kill * T

        dX = k_x_on * eff - k_x_off * X
        dA = k_rec * (1 - A) - k_ex * X * A

        return [dNB, dT, dA, dX]

# ==================== 运行并绘图 ====================
def simulate_untreated(t_end):
    """无治疗肿瘤自然生长对照，用于计算真正的相对肿瘤抑制率 TGI。"""
    def untreated_ode(t, y):
        T = y[0]
        dT = r_T * T * (1 - T / K_T)
        return [dT]

    sol = solve_ivp(
        untreated_ode,
        [0, t_end],
        [1.0],
        max_step=0.25,
        method='LSODA'
    )
    return sol.t, sol.y[0]


# ==================== 运行并绘图 ====================
def run_and_plot():
    # 脉冲组初始状态：[B, E, Depot, NB, T, A, X, J, NB_total]
    y0_p = [0.08, 0.0, 0.0, 0.0, 1.0, A_0, 0.0, 0.0, 0.0]
    ts_p, ys_p = simulate_pulsed([0, t_end], y0_p)

    total_dose = ys_p[-1, 8]
    R_infusion = total_dose / t_end

    # 持续组初始状态：[NB, T, A, X]
    sol_c = solve_ivp(
        ContinuousImmuneModel(R_infusion).model,
        [0, t_end],
        [0.0, 1.0, A_0, 0.0],
        max_step=0.25,
        method='LSODA'
    )
    ts_c = sol_c.t
    ys_c = sol_c.y.T

    # 无治疗对照：用于计算相对肿瘤抑制率
    ts_u, T_u = simulate_untreated(t_end)

    # ==================== 指标计算 ====================
    T0 = 1.0
    T_u_final = T_u[-1]
    T_p_final = ys_p[-1, 4]
    T_c_final = ys_c[-1, 1]

    # 绝对肿瘤回缩率：只在肿瘤小于初始值时为正；不应作为“抑制率”的唯一指标
    regression_p = 1 - T_p_final / T0
    regression_c = 1 - T_c_final / T0

    # 相对肿瘤生长抑制率 TGI：相对于无治疗对照，而不是相对于初始肿瘤大小
    # TGI = 1 - (T_treat - T0) / (T_untreated - T0)
    # 当治疗组仍增长但比无治疗慢时，TGI 仍为正。
    tgi_p = 1 - (T_p_final - T0) / (T_u_final - T0)
    tgi_c = 1 - (T_c_final - T0) / (T_u_final - T0)

    # 也可使用更直观的终点相对负荷：1 - T_treat / T_untreated
    rel_burden_sup_p = 1 - T_p_final / T_u_final
    rel_burden_sup_c = 1 - T_c_final / T_u_final

    peak_p = np.max(ys_p[:, 3])
    peak_c = np.max(ys_c[:, 0])
    avg_p = time_weighted_mean(ts_p, ys_p[:, 3])
    avg_c = time_weighted_mean(ts_c, ys_c[:, 0])

    A_final_p = ys_p[-1, 5]
    A_final_c = ys_c[-1, 2]
    X_final_p = ys_p[-1, 6]
    X_final_c = ys_c[-1, 3]

    print('=' * 70)
    print('新版免疫节律模型：脉冲 vs 持续')
    print('=' * 70)
    print(f'累计释放总剂量: {total_dose:.2e}')
    print(f'持续组输入速率: {R_infusion:.2e} /h')
    print(f'无治疗终点肿瘤大小: {T_u_final:.3f}')
    print(f'脉冲组终点肿瘤大小: {T_p_final:.3f} | 相对TGI: {tgi_p:.1%} | 回缩率: {regression_p:.1%}')
    print(f'持续组终点肿瘤大小: {T_c_final:.3f} | 相对TGI: {tgi_c:.1%} | 回缩率: {regression_c:.1%}')
    print(f'脉冲组终末 T 细胞活性 A: {A_final_p:.3f}')
    print(f'持续组终末 T 细胞活性 A: {A_final_c:.3f}')
    print(f'脉冲组终末持续暴露负担 X: {X_final_p:.3f}')
    print(f'持续组终末持续暴露负担 X: {X_final_c:.3f}')
    print(f'脉冲组峰值浓度: {peak_p:.2e} | 平均浓度: {avg_p:.2e}')
    print(f'持续组峰值浓度: {peak_c:.2e} | 平均浓度: {avg_c:.2e}')
    print('=' * 70)

    # 药效输出
    eff_p = drug_effect(ys_p[:, 3])
    eff_c = drug_effect(ys_c[:, 0])
    j_eff_p = adjuvant_effect(ys_p[:, 7])
    kill_p = k_kill * eff_p * (ys_p[:, 5] / (ys_p[:, 5] + A_threshold)) * (1 + eta_J * j_eff_p)
    kill_c = k_kill * eff_c * (ys_c[:, 2] / (ys_c[:, 2] + A_threshold))

    # ==================== 绘图 ====================
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # 1. 肿瘤生长
    axes[0, 0].plot(ts_u, T_u, color='gray', ls='--', lw=2, label='无治疗对照')
    axes[0, 0].plot(ts_c, ys_c[:, 1], 'g-', lw=2.5, label='持续给药')
    axes[0, 0].plot(ts_p, ys_p[:, 4], 'r-', lw=2.0, label='脉冲释放')
    axes[0, 0].axhline(1.0, color='gray', ls=':', alpha=0.6, label='初始肿瘤大小')
    axes[0, 0].set_title('肿瘤生长抑制')
    axes[0, 0].set_ylabel('肿瘤大小（相对值）')
    axes[0, 0].legend()
    axes[0, 0].grid(alpha=0.3)

    # 2. 药物暴露
    axes[0, 1].plot(ts_c, ys_c[:, 0], 'g-', lw=2.5, label='持续')
    axes[0, 1].plot(ts_p, ys_p[:, 3], 'r-', lw=2.0, label='脉冲')
    axes[0, 1].axhline(nb_EC50, color='purple', ls=':', lw=2, label='EC50')
    axes[0, 1].set_title('药物暴露曲线')
    axes[0, 1].set_ylabel('纳米抗体浓度')
    axes[0, 1].legend()
    axes[0, 1].grid(alpha=0.3)

    # 3. T 细胞活性
    axes[0, 2].plot(ts_c, ys_c[:, 2], 'g-', lw=2.5, label='持续')
    axes[0, 2].plot(ts_p, ys_p[:, 5], 'r-', lw=2.0, label='脉冲')
    axes[0, 2].set_title('T 细胞活性 A')
    axes[0, 2].set_ylabel('T 细胞活性')
    axes[0, 2].legend()
    axes[0, 2].grid(alpha=0.3)

    # 4. 持续暴露负担
    axes[1, 0].plot(ts_c, ys_c[:, 3], 'g-', lw=2.5, label='持续')
    axes[1, 0].plot(ts_p, ys_p[:, 6], 'r-', lw=2.0, label='脉冲')
    axes[1, 0].set_title('持续暴露负担 X')
    axes[1, 0].set_ylabel('暴露负担')
    axes[1, 0].set_xlabel('时间 (h)')
    axes[1, 0].legend()
    axes[1, 0].grid(alpha=0.3)

    # 5. 药效输出
    axes[1, 1].plot(ts_c, kill_c, 'g-', lw=2.5, label='持续杀伤')
    axes[1, 1].plot(ts_p, kill_p, 'r-', lw=1.5, label='脉冲杀伤')
    axes[1, 1].set_title('综合药效输出')
    axes[1, 1].set_ylabel('瞬时杀伤速率')
    axes[1, 1].set_xlabel('时间 (h)')
    axes[1, 1].legend()
    axes[1, 1].grid(alpha=0.3)

    # 6. 关键指标柱状图
    metrics = ['相对TGI', '终末T活性', '峰值浓度', '平均浓度']
    val_p = [tgi_p, A_final_p, peak_p, avg_p]
    val_c = [tgi_c, A_final_c, peak_c, avg_c]

    scale_factor = 0.85
    norm_p = []
    norm_c = []
    for a, b in zip(val_p, val_c):
        m = max(abs(a), abs(b), 1e-12)
        norm_p.append((a / m) * scale_factor)
        norm_c.append((b / m) * scale_factor)

    x = np.arange(len(metrics))
    width = 0.36
    axes[1, 2].bar(x - width/2, norm_p, width, label='脉冲', color='#FF4B4B')
    axes[1, 2].bar(x + width/2, norm_c, width, label='持续', color='#4CAF50')
    axes[1, 2].axhline(0, color='black', lw=0.8)
    axes[1, 2].set_xticks(x)
    axes[1, 2].set_xticklabels(metrics)
    axes[1, 2].set_title('关键指标对比')
    axes[1, 2].set_ylabel('相对值（同指标内归一化）')
    axes[1, 2].legend()
    axes[1, 2].grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    run_and_plot()
