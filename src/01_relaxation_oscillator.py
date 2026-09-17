import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# 设置中文字体（避免乱码）
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS']  # Windows 常用黑体
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# =========================================================
# LEVEL 1:
# 生长-稀释 弛豫振荡器
# =========================================================

# ---------------- 参数 ----------------
mu_max = 0.8          # 最大细菌生长速率 (/h)
K = 1.0               # 环境容纳量（归一化）

alpha_E = 2.5         # 裂解蛋白合成速率
delta_E = 0.2        # 固有降解速率

E_threshold = 8.0     # 裂解阈值
survival = 0.2        # 裂解后存活比例

t_end = 60            # 模拟时长 (h)

# ---------------- 常微分方程组 ----------------
def oscillator_ode(t, y):

    B, E = y

    # Logistic 生长
    growth_rate = mu_max * (1 - B / K)

    # 细菌数量变化
    dBdt = growth_rate * B

    # 裂解蛋白动态
    #
    # 关键：growth_rate * E 表示生长稀释耦合
    #
    dEdt = alpha_E - (delta_E + growth_rate) * E

    return [dBdt, dEdt]

# ---------------- 事件：触发裂解 ----------------
def lysis_event(t, y):

    return y[1] - E_threshold

lysis_event.terminal = True
lysis_event.direction = 1

# ---------------- 裂解重置 ----------------
def apply_lysis(y):

    B, E = y

    # 大部分细菌裂解死亡
    B_new = survival * B

    # 胞内裂解蛋白急剧下降
    E_new = 0.0

    return [B_new, E_new]

# ---------------- 仿真循环 ----------------
def simulate():

    t_current = 0.0

    # 初始状态
    y_current = [0.05, 0.0]

    ts = [t_current]
    ys = [y_current]

    max_events = 100
    event_count = 0

    while t_current < t_end and event_count < max_events:

        sol = solve_ivp(
            oscillator_ode,
            [t_current, t_end],
            y_current,
            events=lysis_event,
            max_step=0.1,
            method='LSODA'
        )

        # 保存轨迹
        if len(sol.t) > 1:
            ts.extend(sol.t[1:])
            ys.extend(sol.y[:,1:].T)

        # 更新当前状态
        t_current = sol.t[-1]
        y_current = sol.y[:,-1]

        # 如果发生裂解
        if len(sol.t_events[0]) > 0:

            t_event = sol.t_events[0][0]
            y_event = sol.y_events[0][0]

            # 应用裂解重置
            y_current = apply_lysis(y_event)

            # 保存间断点
            ts.append(t_event)
            ys.append(y_current)

            t_current = t_event

            event_count += 1

        else:
            break

    return np.array(ts), np.array(ys)

# ---------------- 运行仿真 ----------------
ts, ys = simulate()

B = ys[:,0]
E = ys[:,1]

# =========================================================
# 可视化（中文标签）
# =========================================================

fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)

# ---- 细菌数量 ----
axes[0].plot(ts, B, linewidth=2.5)

axes[0].set_ylabel("归一化细菌浓度")
axes[0].set_title("生长-稀释弛豫振荡器")
axes[0].grid(alpha=0.3)

# ---- 裂解蛋白 ----
axes[1].plot(ts, E, linewidth=2.5)

axes[1].axhline(
    E_threshold,
    linestyle='--',
    alpha=0.7,
    label='裂解阈值'
)

axes[1].set_ylabel("裂解蛋白 E")
axes[1].set_xlabel("时间 (小时)")
axes[1].legend()
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.show()