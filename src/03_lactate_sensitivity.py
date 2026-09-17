import numpy as np
import matplotlib.pyplot as plt

# 设置中文字体，避免乱码
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ==================== 乳酸敏感性分析 ====================
# 扫描乳酸浓度，评估振荡可行性风险

L_values = np.linspace(2.0, 10.0, 50)
E_steady = []
can_lyse = []

for L in L_values:
    # 稳态 E（保守估计，忽略生长稀释）
    # ALPaGA 活性建模为 Hill 函数，EC50 = 5.0 mM
    alpha = 0.01 + (3.0 - 0.01) * (L**2 / (L**2 + 5.0**2))
    E_ss = alpha / 0.2  # delta_E = 0.2 /h
    E_steady.append(E_ss)
    can_lyse.append(E_ss > 8.0)  # 阈值 = 8.0 a.u.

E_steady = np.array(E_steady)
can_lyse = np.array(can_lyse)

# ==================== 可视化 ====================
plt.figure(figsize=(8, 5))
plt.plot(L_values, E_steady, 'b-', linewidth=2.5, label='稳态 E')
plt.axhline(8.0, color='red', linestyle='--', linewidth=2, label='裂解阈值 (8.0 a.u.)')
plt.axvline(5.0, color='gray', linestyle=':', linewidth=2, label='ALPaGA EC50 (5.0 mM)')
plt.fill_between(L_values, 0, E_steady, where=can_lyse, color='green', alpha=0.2, label='振荡区域')
plt.fill_between(L_values, 0, E_steady, where=~can_lyse, color='red', alpha=0.2, label='稳态陷阱')
plt.xlabel('局部乳酸浓度 (mM)', fontsize=12)
plt.ylabel('稳态裂解蛋白 E (a.u.)', fontsize=12)
plt.title('风险评估：乳酸阈值与振荡可行性', fontsize=14, fontweight='bold')
plt.legend(fontsize=10)
plt.grid(alpha=0.3)
plt.ylim(0, max(E_steady) * 1.1)
plt.tight_layout()
plt.show()

# ==================== 控制台输出 ====================
# 寻找振荡变得可行的临界乳酸浓度
critical_idx = np.argmax(can_lyse)  # 第一个 True
critical_L = L_values[critical_idx] if np.any(can_lyse) else None

# 计算在 ALPaGA EC50 (L=5.0 mM) 处的 E_ss
idx_EC50 = np.argmin(np.abs(L_values - 5.0))
E_at_EC50 = E_steady[idx_EC50]

print("=" * 60)
print("乳酸敏感性分析")
print("=" * 60)
print(f"ALPaGA EC50: 5.0 mM")
print(f"裂解阈值: 8.0 a.u.")
print(f"L = 5.0 mM 时的稳态 E: {E_at_EC50:.2f} a.u.")
if critical_L is not None:
    print(f"振荡所需的临界乳酸浓度: L >= {critical_L:.2f} mM")
    print(f"  -> 低于此值：稳态陷阱 (E < 8.0)")
    print(f"  -> 高于此值：振荡区域 (E >= 8.0)")
else:
    print("在任何扫描的乳酸浓度下均无法触发振荡。")
print("-" * 60)
print("解读：")
if critical_L is not None and critical_L > 5.0:
    print(f"  警告：ALPaGA 激活 (L>=5.0 mM) 与振荡可行性 (L>={critical_L:.2f} mM) 之间存在间隙。")
    print("  这就是“稳态陷阱”——代谢自调节必须确保乳酸维持在此临界值以上。")
elif critical_L is not None and critical_L <= 5.0:
    print("  一旦 ALPaGA 被激活，振荡即可发生。")
    print("  无明显稳态陷阱风险。")
print("=" * 60)