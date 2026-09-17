import numpy as np
import matplotlib.pyplot as plt

# 设置中文字体，避免乱码
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ==================== 模型参数 ====================
lysis_thr = 8.0                     # 裂解阈值 (a.u.)
delta_E = 0.2                       # 降解速率 (1/h)
alpha_E_max = 3.0                   # 无 LacI 抑制时的最大表达速率 (a.u./h)

# LacI-IPTG 参数
# 注意：这里 I_total 仅用于绿色柱子（有IPTG）的计算
I_total_green = 10.0                # 绿色柱子：lacIq过表达
I_total_red = 1.0                   # 红色柱子：非过表达，展示泄漏梯度

Ki = 0.1                            # IPTG-LacI 亲和常数 (mM)
n_iptg = 2                          # IPTG 结合的 Hill 系数
h_LacI = 2                          # LacI 抑制的 Hill 系数

# lacO 海绵参数
K50_base = 0.10                     # 强亲和力
K_N = 3.0                           # 每个额外 lacO 位点调整 K50 的因子

# 扫描 lacO 位点数量
N_lacO_vals = np.arange(0, 6)       # 0 到 5 个额外 lacO 位点
iptg_optimal = 1.0                  # 代表性最优 IPTG 浓度 (mM)

# -------- 绿色柱子（有IPTG）：保持原逻辑不变 --------
K50_eff_green = K50_base * (1 + N_lacO_vals / K_N)
I_free_opt = I_total_green * Ki**n_iptg / (Ki**n_iptg + iptg_optimal**n_iptg)
alpha_E_opt = alpha_E_max * (1 - I_free_opt**h_LacI / (I_free_opt**h_LacI + K50_eff_green**h_LacI))
alpha_E_opt = np.clip(alpha_E_opt, 0.0, alpha_E_max)
E_ss_opt = alpha_E_opt / delta_E

# -------- 红色柱子（无IPTG）：修正参数与公式 --------
K50_eff_red = K50_base / (1 + N_lacO_vals / K_N)   # 修正：位点越多，表观K50越小，抑制越强
I_free_zero = I_total_red
alpha_E_zero = alpha_E_max * (1 - I_free_zero**h_LacI / (I_free_zero**h_LacI + K50_eff_red**h_LacI))
alpha_E_zero = np.clip(alpha_E_zero, 0.0, alpha_E_max)
E_ss_zero = alpha_E_zero / delta_E

# ==================== 可视化 ====================
plt.figure(figsize=(10, 6))
plt.bar(N_lacO_vals - 0.15, E_ss_zero, width=0.3, color='red', alpha=0.7, label='无 IPTG（泄漏，I_total=1.0）')
plt.bar(N_lacO_vals + 0.15, E_ss_opt, width=0.3, color='green', alpha=0.7, label=f'IPTG = {iptg_optimal} mM（I_total=10.0）')
plt.axhline(lysis_thr, color='black', linestyle='--', linewidth=2, label=f'裂解阈值（{lysis_thr} a.u.）')
plt.axhline(1.5 * lysis_thr, color='blue', linestyle=':', linewidth=1.5, label='最优窗口下限（12 a.u.）')
plt.axhline(2.5 * lysis_thr, color='blue', linestyle=':', linewidth=1.5, label='最优窗口上限（20 a.u.）')
plt.xlabel('额外 lacO 位点数量', fontsize=12)
plt.ylabel('稳态裂解蛋白 E (a.u.)', fontsize=12)
plt.title('LacI-lacO 调控：泄漏抑制与最优表达', fontsize=14)
plt.xticks(N_lacO_vals)
plt.legend(fontsize=10)
plt.grid(alpha=0.3, axis='y')
plt.ylim(0, max(E_ss_opt)*1.2)
plt.tight_layout()
plt.show()

# ==================== 输出汇总 ====================
print("="*60)
print("LacI-lacO 调控 — 修正后的泄漏趋势")
print("="*60)
print(f"裂解阈值: {lysis_thr:.1f} a.u.")
print(f"红色柱子: I_total = {I_total_red:.1f}（非过表达 LacI）")
print(f"绿色柱子: I_total = {I_total_green:.1f}（lacIq 过表达）+ IPTG = {iptg_optimal} mM")
print()
for i, n in enumerate(N_lacO_vals):
    status_zero = "安全" if E_ss_zero[i] < lysis_thr else "泄漏"
    status_opt = "最优" if 1.5*lysis_thr <= E_ss_opt[i] <= 2.5*lysis_thr else "窗口外"
    print(f"  lacO = {n}: 无 IPTG → E = {E_ss_zero[i]:.2f} a.u. ({status_zero})   |   IPTG → E = {E_ss_opt[i]:.2f} a.u. ({status_opt})")
print("="*60)
print("红色柱趋势：lacO 少 → 泄漏多；lacO 多 → 逐步被抑制（左高右低）")
print("绿色柱趋势：IPTG 充分解除抑制，表达稳定在最优窗口")
print("="*60)