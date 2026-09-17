import numpy as np
import matplotlib.pyplot as plt

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS']  # Windows 常用黑体
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# ==================== 参数 ====================
beta_A = 0.1             # Asd / phiX174E 降解速率 (min⁻¹), 半衰期 ~6.9 min
beta_E = 0.1
gamma_DAP_eff = 0.12     # DAP 有效降解速率 (min⁻¹)
lysis_threshold = 8.0    # phiX174E 裂解阈值 (a.u.)
DAP_crit = 0.2           # DAP 临界浓度 (a.u.)
k_death = 1.5            # DAP 耗尽后细菌死亡率 (min⁻¹)

A_0 = 18.82              # 初始 Asd
E_0 = 18.82              # 初始 phiX174E
DAP_0 = 20.0             # 初始 DAP
B_0 = 1.0                # 初始细菌数量 (归一化)

# 时间轴
t = np.linspace(0, 60, 500)

# 解析衰减
A = A_0 * np.exp(-beta_A * t)
E = E_0 * np.exp(-beta_E * t)
DAP = DAP_0 * np.exp(-gamma_DAP_eff * t)

# 关键时间点
t_E_safe = -np.log(lysis_threshold / E_0) / beta_E
t_DAP_crit = -np.log(DAP_crit / DAP_0) / gamma_DAP_eff
t_A_10pct = -np.log(0.1) / beta_A

# 细菌数量：DAP耗尽后开始指数死亡
B = np.ones_like(t)
death_mask = t >= t_DAP_crit
B[death_mask] = B_0 * np.exp(-k_death * (t[death_mask] - t_DAP_crit))
B[t > t_DAP_crit + 5] = 0  # 5分钟后视为完全清除

# ==================== 鲁棒性扫描数据 ====================
half_lives = np.linspace(5, 60, 100)
beta_vals = np.log(2) / half_lives
E_tumor_ss = 20.0
t_lysis_loss = -np.log(lysis_threshold / E_tumor_ss) / beta_vals
t_asd_depletion = -np.log(0.1) / beta_vals
safety_window = t_asd_depletion - t_lysis_loss

# 寻找"安全窗口归零"的临界半衰期。
# 注意：不能写 critical_idx = np.argmax(safety_window < 0) —— 当没有任何元素 < 0 时
# argmax 会返回 0，从而被误解读为"临界值等于第一个扫描点"。
_closed = np.where(safety_window < 0)[0]
if _closed.size > 0:
    critical_hl = half_lives[_closed[0]]
    critical_msg = f"{critical_hl:.1f} 分钟（安全窗口归零）"
else:
    critical_msg = (f"在 {half_lives[0]:.0f}–{half_lives[-1]:.0f} 分钟的扫描范围内"
                    f"安全窗口始终为正，未出现归零点")
design_idx = np.argmin(np.abs(half_lives - 7))

# ==================== 绘图 ====================
fig, axes = plt.subplots(1, 2, figsize=(18, 6))

# ---- 左图：衰减曲线 + 细菌数量 (双Y轴) ----
ax1 = axes[0]
ax1.plot(t, A, 'b-', linewidth=2.5, label='Asd (酶)')
ax1.plot(t, E, 'r-', linewidth=2.5, label='phiX174E (裂解蛋白)')
ax1.plot(t, DAP, 'g-', linewidth=2.5, label='DAP (细胞壁成分)')

ax1.axhline(lysis_threshold, color='red', linestyle='--', alpha=0.6)
ax1.text(62, lysis_threshold, f'{lysis_threshold}', color='red', fontsize=9, va='center')
ax1.axhline(DAP_crit, color='green', linestyle='--', alpha=0.6)
ax1.text(62, DAP_crit, f'{DAP_crit}', color='green', fontsize=9, va='center')

# 关键事件标注
ax1.axvline(t_E_safe, color='red', linestyle=':', alpha=0.7)
ax1.annotate(f'裂解无法触发\n(t={t_E_safe:.1f} 分钟)', xy=(t_E_safe, lysis_threshold),
            xytext=(t_E_safe+8, lysis_threshold+3), fontsize=9, color='red', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='red'))

ax1.axvline(t_DAP_crit, color='green', linestyle=':', alpha=0.7)
ax1.annotate(f'细胞壁缺陷\n(t={t_DAP_crit:.1f} 分钟)', xy=(t_DAP_crit, DAP_crit),
            xytext=(t_DAP_crit+8, DAP_crit+2), fontsize=9, color='green', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='green'))

# 安全窗口
ax1.axvspan(t_E_safe, t_DAP_crit, alpha=0.1, color='green')
ax1.text((t_E_safe+t_DAP_crit)/2, 17, f'安全窗口\n{t_DAP_crit-t_E_safe:.1f} 分钟',
        ha='center', fontsize=11, color='darkgreen', fontweight='bold',
        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

ax1.set_xlabel('进入正常组织后的时间 (分钟)', fontsize=12)
ax1.set_ylabel('蛋白质 / 代谢物浓度 (任意单位)', fontsize=12)
ax1.set_xlim(0, 65)
ax1.set_ylim(0, 22)

# 右侧Y轴：细菌数量
ax1b = ax1.twinx()
ax1b.plot(t, B, 'k--', linewidth=2.5, label='存活细菌')
ax1b.set_ylabel('存活细菌 (相对值)', fontsize=12)
ax1b.set_ylim(0, 1.2)
ax1b.legend(loc='center right', fontsize=9)

# 合并图例
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax1b.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=9)
ax1.set_title('安全模块：降解与细菌死亡', fontsize=14, fontweight='bold')
ax1.grid(alpha=0.3)

# ---- 右图：鲁棒性扫描 ----
ax2 = axes[1]
ax2.plot(half_lives, t_lysis_loss, 'r-', linewidth=2, label='裂解失效 (E < 8.0)')
ax2.plot(half_lives, t_asd_depletion, 'b--', linewidth=2, label='Asd 耗尽 (< 10%)')
ax2.fill_between(half_lives, t_lysis_loss, t_asd_depletion, alpha=0.12, color='green', label='安全窗口')
ax2.axvline(7, color='gray', linestyle=':', alpha=0.7, label='设计值 (半衰期=7 分钟)')
ax2.axhline(0, color='red', linestyle='--', alpha=0.5, linewidth=1)

hl_design = half_lives[design_idx]
win_design = safety_window[design_idx]
ax2.annotate(f'安全窗口 = {win_design:.1f} 分钟', xy=(hl_design, (t_lysis_loss[design_idx]+t_asd_depletion[design_idx])/2),
            xytext=(hl_design+10, (t_lysis_loss[design_idx]+t_asd_depletion[design_idx])/2+8),
            fontsize=10, color='darkgreen', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='darkgreen'))

ax2.set_xlabel('蛋白质半衰期 (分钟)', fontsize=12)
ax2.set_ylabel('信号丢失后时间 (分钟)', fontsize=12)
ax2.set_title('安全窗口的鲁棒性', fontsize=14, fontweight='bold')
ax2.legend(loc='upper left', fontsize=9)
ax2.grid(alpha=0.3)
ax2.set_xlim(5, 60)
ax2.set_ylim(0, 100)

plt.tight_layout()
plt.show()

# ==================== 输出关键数值 ====================
print(f"设计半衰期 (7 分钟): 裂解失活 {t_lysis_loss[design_idx]:.1f} 分钟, "
      f"Asd耗尽 {t_asd_depletion[design_idx]:.1f} 分钟, "
      f"安全窗口 {safety_window[design_idx]:.1f} 分钟")
print(f"临界半衰期: {critical_msg}")
print(f"DAP 耗尽后细菌清除时间: ~5 分钟 (死亡率 k_death = {k_death} min⁻¹)")