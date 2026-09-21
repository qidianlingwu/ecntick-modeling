# EcNTick 建模：乳酸驱动的工程菌脉冲释药系统

> 2026 年四川大学微生物应用设计大赛参赛项目的**机理建模（mechanistic modeling）**部分。
> 项目全称：**EcNTick — 用于结直肠癌免疫治疗的自主节律益生菌节拍器**

---

## 代码说明

用**事件驱动的常微分方程（ODE）+ 随机个体模型**，验证一个工程菌回路在动力学上是否可行：
> 乳酸感知 → 裂解蛋白积累到阈值 → 周期性裂解 → 脉冲式释放抗 PD-L1 纳米抗体

**这是一个纯机理模型，不包含任何真实组学或实验数据。** 所有参数为设计取值或量级估计，
逐条来源见 [`PARAMETERS.md`](PARAMETERS.md)。

---

## 个人贡献

| 项 | 说明 |
|---|---|
| 时间 | 2026-04-23 → 2026-05-23，约一个月 |
| 规模 | 35 个脚本版本，最终归档 6 个模块（`src/`）+ 6 个被取代的早期版本（`legacy/`）|
| 建模模块 | 弛豫振荡 / PK-PD 免疫节律 / 乳酸敏感性 / LacI-lacO 阻尼 / 安全窗鲁棒性 / 随机个体模型 |

### 技术要点

1. **事件驱动 ODE + 状态重置**
   用 `scipy.integrate.solve_ivp` 的 `events` 参数做终末事件检测
   （`terminal=True, direction=1`），触发裂解后用 `apply_lysis()` 重置状态，
   再以新初值继续积分，循环至多 200 次。关键是处理**间断点 + 分段重积分**。

2. **等总剂量下的公平对比**
   脉冲组与持续给药组按**相同累计剂量**比较，避免"脉冲组靠堆剂量取胜"：

   ```python
   total_dose = ys_p[-1, 8]          # 脉冲组累计释放总量
   R_infusion = total_dose / t_end   # 持续组按同样总量换算输入速率
   ```

3. **时间加权均值（非等间距采样）**
   事件驱动模型的时间网格是非等间距的，因此**不能用 `np.mean`**：

   ```python
   def time_weighted_mean(t, y):
       return np.trapezoid(y, t) / (t[-1] - t[0])
   ```

4. **相对 TGI 用无治疗对照，而不是初始肿瘤大小**
   专门实现了 `simulate_untreated()` 作为对照。
   绝对回缩率在肿瘤持续生长时会变成负数，不能当作"抑制率"使用 ——
   两个指标在代码中并列输出，标签明确区分。

5. **参数扫描与鲁棒性分析**
   乳酸浓度、IPTG × lacO 二维网格、蛋白半衰期、杀伤/表达速率等均做了扫描。

---

## 失败与修正

**第一版模型失败了。**

最初只构建了基础 PK/PD 模型（只考虑药物浓度对肿瘤的直接杀伤）。
结果发现：**只看平均药物暴露的话，持续给药也能产生明显抑瘤效果**——
这意味着模型无法支持"脉冲优于持续"这个项目核心假设。

于是把模型升级为 **PK/PD-免疫节律耦合模型**，补入三个变量：

| 变量 | 含义 |
|---|---|
| `A` | T 细胞活性（决定免疫杀伤能否真正发挥） |
| `X` | 持续暴露负担（长期受体占用导致的免疫疲劳） |
| `J` | 裂解佐剂信号（菌体裂解释放 PAMPs 带来的局部免疫激活） |

脉冲释放的优势因此不再是"峰值更高"，而是体现在三处：
**周期性高峰浓度 + 低谷期免疫恢复窗口 + 裂解佐剂效应**。

> 相关代码：`src/02_pkpd_immune_rhythm.py`；
> 中间版本：`legacy/pkpd_immune_rhythm_v1_no_depot.py`

---

## 这套模型**不能**证明什么

- 不能证明实际工程菌一定能稳定振荡 —— 需要 OD600 / CFU 实验验证
- 不能预测绝对药物浓度 —— 参数多为量级估计，单位是归一化的 a.u.
- 不能替代疗效实验 —— 全部为机理推演，无任何真实数据
- 不能作为安全性结论 —— 安全窗口是模型预测，需体内分布与毒性实验
- 随机模拟部分（`src/06`）建模的是**早期 QS 依赖设计**，非最终设计

完整边界说明见 [`LIMITATIONS.md`](LIMITATIONS.md)。

---

## 目录结构

```
ecntick-modeling/
├── README.md                 ← 本文件
├── PARAMETERS.md             ← 逐参数来源表（哪些是估计、哪些待实验标定）
├── LIMITATIONS.md            ← 模型边界与已知缺陷
├── HISTORY.md                ← 一个月、35 个版本的迭代说明
├── requirements.txt
├── src/                      ← 最终归档的 6 个模块
│   ├── 01_relaxation_oscillator.py       弛豫振荡：能否周期性裂解
│   ├── 02_pkpd_immune_rhythm.py          PK/PD-免疫节律：脉冲 vs 持续（核心）
│   ├── 03_lactate_sensitivity.py         乳酸阈值与"稳态陷阱"
│   ├── 04_laci_laco_damper.py            LacI-lacO 阻尼：泄露 vs 可激活
│   ├── 05_safety_window_robustness.py    安全窗口与半衰期鲁棒性
│   └── 06_stochastic_individual_based.py 随机个体模型（非 Gillespie，见文件头）
├── legacy/                   ← 被取代的早期版本，保留以供追溯
└── figs/                     ← 14 张输出图
```

---

## 运行说明

```bash
conda create -n ecntick python=3.11 -c conda-forge
conda activate ecntick
pip install -r requirements.txt

cd src
python 01_relaxation_oscillator.py
python 02_pkpd_immune_rhythm.py
# ...
```

各脚本运行后会弹出 matplotlib 图窗，并在终端打印关键指标。
`src/06_stochastic_individual_based.py` 会在当前目录写出 `gillespie_output.png`。

**依赖**：Python ≥ 3.9，numpy ≥ 2.0（使用了 `np.trapezoid`，旧版请改为 `np.trapz`）、scipy、matplotlib。

---

## 已知的命名遗留问题

1. 原文件名 `Gillespie.py` / `re-Gillespie.py` / `Gillespie(1).py` **命名不当**：
   实际实现的是固定步长随机个体模型，**不是 Gillespie SSA**（详见 `src/06` 文件头说明）。
2. `src/06` 输出文件仍叫 `gillespie_output.png`，未改名。
3. `src/02` 中原有的"目标是显示…"类注释已改写为取值依据说明（原注释见 `HISTORY.md`）。

---

## 致谢与协作说明

本项目建模工作由本人负责（模型设计、参数选择、结果判断与结论边界）；
调试过程中使用了 AI 辅助工具。所有方程、参数取值理由与每张图的含义均可解释。

---

## 引用

若需引用本项目：

> EcNTick: 用于结直肠癌免疫治疗的自主节律益生菌节拍器. 2026 年四川大学微生物应用设计大赛参赛项目（机理建模部分）. 2026.
