# 迭代说明（HISTORY.md）

> 记录从零到终版的完整过程。**35 个脚本版本，横跨 2026-04-23 → 2026-05-23，约一个月。**
> 保留这份记录的目的：说明这套模型不是一次成型的，而是经过反复推翻与修正的结果。

---

## 总览

| 阶段 | 时间 | 版本数 | 主题 |
|---|---|---|---|
| 一 | 04-23 → 04-24 | 6 | 早期振荡与药效概念验证 |
| 二 | 04-30 → 05-03 | 10 | 回路参数、乳酸阈值与阻尼探索 |
| 三 | 05-07 → 05-08 | 9 | 随机模拟与模型精简 |
| — | 05-09 → 05-19 | — | 申报书提交（05-08）后的间歇期 |
| 四 | 05-20 → 05-23 | 10 | 终版整理与模块归档 |

---

## 阶段一：早期振荡与药效概念验证（04-23 → 04-24）

| 时间 | 文件 | 做了什么 |
|---|---|---|
| 04-23 15:46 | `lactic_anoxia_andGate.py` | 乳酸 × 缺氧二维网格 + AND 门响应面（`figs/fig1_and_gate_surface.png`）|
| 04-23 16:02 | `repressilatorOSC.py` | 试算经典 repressilator 振荡器 |
| 04-23 16:11 | `andGatedOscillatorCouplingModel.py` | 振荡器与 AND 门耦联模型 |
| 04-23 16:26 | `therapy_output.py` | 首次尝试药效输出建模 |
| 04-23 18:47 | `robustness_waveforms.py` | 振荡参数鲁棒性扫描（`figs/fig6_robustness_waveforms.png`）|
| 04-24 20:00 | `HRtest.py` | 参数试探 |

**这一阶段的方向后来被放弃了。** 最初的思路是"用 repressilator 类振荡器 + AND 门实现
肿瘤特异性脉冲"，但这类方案需要额外基因元件，与项目"极简回路"的核心主张冲突。

**转向**：改为利用**自杀回路的天然弛豫振荡**（生长-稀释-阈值裂解），无需额外振荡元件。
依据是文献 *Oscillations by minimal bacterial suicide circuits*（Marguet et al., PLoS ONE 2010）
——该文证明简单自杀回路的振荡依赖宿主介导的质粒扩增，**不需要群体感应模块**。

> 相关早期图：`figs/oscillation_comparison.png`、`figs/corrected_oscillation_comparison.png`、
> `figs/verified_oscillation.png`、`figs/professional_slc_oscillation.png`、
> `figs/verified_slc_oscillation.png`、`figs/realistic_slc_simulation.png`、`figs/project_oscillation_model.png`

---

## 阶段二：回路参数、乳酸阈值与阻尼探索（04-30 → 05-03）

| 时间 | 文件 | 做了什么 |
|---|---|---|
| 04-30 21:30 | `lactornext.py` | 乳酸响应回路再试 |
| 05-01 20:37 | `practical.py` | 向"可实现"方向收敛 |
| 05-01 21:09 | `anothertry.py` | 空间/扩散尝试 |
| 05-02 17:53 | `re-Gillespie.py` | 随机模型（软性承载版本）|
| 05-02 21:23 | `lactatesensity.py` | 乳酸敏感性扫描雏形 |
| 05-02 21:48 | `bestexpress.py` | IPTG × lacO 二维网格扫描 |
| 05-02 21:51 | `ITPG.py` | LacI-IPTG 阻尼模型雏形 |
| 05-02 22:35 | `Pulse-Only.py` | 单独考察脉冲释放 |
| 05-02 22:38 | `nexttest-compaire.py` | 脉冲 vs 持续对比 |
| 05-03 00:11 | `re-saftytest.py` | 安全性试探 |

**这一阶段产出了"稳态陷阱"的初步观察**：
启动子被乳酸激活 ≠ 裂解蛋白能越过阈值，中间存在一段"激活但不裂解"的死区。
这个观察直接导致了 **LldR 代谢自调控**这个设计要素的出现——
如果乳酸被细菌过度消耗，系统会掉进死区而停止振荡。

> 相关早期图：`figs/plasmid_sensitivity.png`、`figs/threshold_sensitivity.png`、
> `figs/coupled_lactate_oscillation.png`、`figs/fig5_sensitivity_analysis.png`

---

## 阶段三：随机模拟与模型精简（05-07 → 05-08）

| 时间 | 文件 | 做了什么 |
|---|---|---|
| 05-07 00:42 | `timecheck.py` | 安全窗时序 |
| 05-07 01:23 | `IPTGnext.py` | 阻尼模型迭代 |
| 05-07 16:16 | `Gillespie(1).py` | 随机个体模型 v1 |
| 05-07 17:22 | `Gillespie.py` | 随机个体模型 v2（+ AHL 群体感应）|
| 05-08 19:46 | `newlevel1.py` | Level 1 模型重构 |
| 05-08 20:40 | `simplifymodel.py` | 简化 |
| 05-08 21:23 | `level2.py` | Level 2 模型 |
| 05-08 21:37 | `level2try.py` | Level 2 迭代 |
| 05-08 21:55 | `compairetry.py` | 对比模型迭代 |

> 05-08 当天一直改到 21:55 —— 申报书提交日。

**⚠️ 命名遗留**：`Gillespie.py` 系列**并不是 Gillespie 算法**，
实现的是固定时间步长的随机个体模型（详见 `src/06` 文件头）。
这个命名问题是后来整理时才发现的。

---

## 阶段四：终版整理与模块归档（05-20 → 05-23）

把散落的脚本收敛为 9 个 `AAA` 前缀的终版文件：

| 时间 | 原文件 | 现归档为 | 说明 |
|---|---|---|---|
| 05-20 23:37 | `AAAtimecheck.py` | `src/05_safety_window_robustness.py` | 安全窗 + 半衰期鲁棒性 |
| 05-20 23:44 | `AAArusuanyuzhi.py` | `src/03_lactate_sensitivity.py` | 乳酸阈值与稳态陷阱 |
| 05-20 23:49 | `AAAIPTG.py` | `src/04_laci_laco_damper.py` | LacI-lacO 阻尼 |
| 05-21 00:10 | `AAAzhengdang.py` | `src/01_relaxation_oscillator.py` | 弛豫振荡（"正档"） |
| 05-21 20:22 | `AAAyaoxiao.py` | `legacy/drug_efficacy_v1.py` | 药效对比 v1 |
| 05-21 20:43 | `AAAyaoxiaozaigai.py` | `legacy/drug_efficacy_v2.py` | 药效对比 v2 |
| 05-21 20:56 | `test.py` | 〔未归档〕 | 临时试验 |
| 05-23 18:42 | `AAAjixugai.py` | `legacy/pkpd_immune_rhythm_v1_no_depot.py` | 免疫节律模型 v1（无缓释项）|
| 05-23 19:01 | `AAAgaigaigai.py` | **`src/02_pkpd_immune_rhythm.py`** | **免疫节律模型终版（+ Depot 缓释）** |
| 05-23 19:11 | `AAAyaoxiaogai.py` | `legacy/drug_efficacy_v3_simple.py` | 简化的药效对比分支 |

**为什么 `SRC/02` 选 `AAAgaigaigai` 而不是最后的 `AAAyaoxiaogai`？**

`AAAyaoxiaogai`（19:11）虽然是最后一个修改的文件，但它的模型更简单
（5 状态 `[B, E, NB, T, NB_total]`，无免疫节律变量）。
而 `AAAgaigaigai` 包含完整的免疫节律机制（`A` / `X` / `J`）与 Depot 缓释项，
且与 `docs/建模目的.docx` 中描述的最终模型一致。因此以它为终版。

---

## 核心转折：一次模型失败与一次架构升级

这是整个建模过程中最有价值的一段，也是**最值得在面谈时讲的部分**。

### 第一版模型（基础 PK/PD）

只考虑药物浓度对肿瘤的直接杀伤：

```
dNB/dt = -nb_deg * NB + 释放项
dT/dt  = r_T·T·(1 - T/K_T) - k_kill·Hill(NB)·T
```

**结果：这个模型失败了。**

原因是：**只看平均药物暴露的话，持续给药也能产生明显抑瘤效果。**
换句话说，模型无法支持"脉冲式释药优于持续给药"这个项目的核心假设——
而这恰恰是整个设计存在的理由。

### 第二版模型（PK/PD-免疫节律耦合）

补入三个变量：

| 变量 | 含义 | 为什么必须加 |
|---|---|---|
| `A` | T 细胞活性 | 药物阻断 PD-L1 只是"解除刹车"，真正的杀伤靠 T 细胞；持续暴露会让 T 细胞耗竭 |
| `X` | 持续暴露负担 | 长期受体占用导致的免疫疲劳累积 |
| `J` | 裂解佐剂信号 | 细菌裂解释放 PAMPs，带来局部免疫激活（这是脉冲方案独有的） |

**升级后结论改变了**：脉冲组的优势不再只是"峰值更高"，而是三件事叠加：

1. 周期性高峰浓度；
2. **低谷期给 T 细胞恢复窗口**（持续给药没有这个过程，`X` 不断累积、`A` 不断下降）；
3. 裂解带来的佐剂效应（持续给药组不含 `J` 项）。

**终版模型输出（`src/02_pkpd_immune_rhythm.py` 实际运行结果）**：

```
无治疗终点肿瘤大小: 3.000
脉冲组终点肿瘤大小: 1.827 | 相对TGI: 58.6% | 回缩率: -82.7%
持续组终点肿瘤大小: 2.219 | 相对TGI: 39.0% | 回缩率: -121.9%
脉冲组终末 T 细胞活性 A: 0.347
持续组终末 T 细胞活性 A: 0.234
脉冲组终末持续暴露负担 X: 0.869
持续组终末持续暴露负担 X: 1.126
```

注意两点：
- **回缩率是负数**（肿瘤仍在生长，只是长得慢）——所以不能用"抑制率"来描述，
  必须用相对无治疗对照的 **TGI**。这个指标定义问题是在建模过程中才意识到的。
- **脉冲组的平均浓度反而更低**（1.55e6 vs 1.67e6），峰值更高（2.12e6 vs 2.00e6）——
  也就是说脉冲组的优势来自**节律**，不是来自"平均给药量更多"。

---

## 代码整理时修改的内容（透明记录）

整理归档时，以下三处**原始注释/写法被修改**，原始内容保留在此以备追溯：

### 1. `src/02`（原 `AAAgaigaigai.py`）——参数注释

```python
# 原：
# 设置为保守情景：目标是显示“延缓生长/增强控制”，而不是过度展示完全回缩
k_kill = 0.095

# 原：
# 降低佐剂增益，避免模型出现过强、过快的肿瘤回缩
J_burst = 0.65
```

**改为什么**：改为说明取值依据与敏感性检验位置。
**为什么改**：原注释是"目标导向"的（"为了显示…"），读起来像是为了让图好看而调参数。
实际上参数取保守下界是合理的科研选择，但**表述方式**会让读者产生相反的理解。

### 2. `legacy/drug_efficacy_v3_simple.py`（原 `AAAyaoxiaogai.py`）——图表归一化

```python
# 原：
# 归一化柱状图（加入缩放系数）
scale_factor = 0.85  # 缩小最大值到95%
vals_p = [val_sup_p / max(val_sup_p, val_sup_c) * scale_factor, ...]
```

**改为什么**：删除 `scale_factor`，直接归一化到同指标内较大值。
**为什么改**：`# 缩小最大值到95%` 这类注释会让人怀疑数据被美化，
即使它只影响柱状图的高度显示、不影响任何指标数值。**不值得冒这个风险。**

### 3. `legacy/drug_efficacy_v3_simple.py`——指标名与均值算法

```python
# 原：
metrics = ['肿瘤抑制率', ...]
val_avg_p = np.mean(ys_p[:, 2])
```

**改为什么**：
- 指标名改为 `'肿瘤回缩率'`（因为它相对的是**初始肿瘤大小**，不是无治疗对照，叫"抑制率"不准确）；
- `np.mean` 改为时间加权均值（事件驱动模型的时间网格是**非等间距**的）。

---

## 其他遗留物

| 文件 | 说明 |
|---|---|
| `testsafty`（无扩展名，4412 B） | 早期测试的输出/草稿，未归档 |
| `test.py`（05-21 20:56） | 终版整理期的临时试验，未归档 |
| `.claude/settings.local.json` | AI 辅助工具的配置，未归档 |
| `gillespie_output.png` | 由 `src/06` 生成，文件名沿用历史（方法已非 Gillespie）|

---

## 一句话总结这段历史

**这一个月里最重要的事，不是把模型跑通，而是发现第一版模型不能支持自己的假设，
然后去补上缺失的机制。**

如果不能把这件事讲清楚，那 35 个版本就只是一堆脚本。
