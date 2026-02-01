# 磁性体系计算工作流

以 PdCrO2 钒掺杂磁性计算为例，展示 VASP 非共线磁性计算的完整工作流程。

## 📖 目录结构

```
magnetic-systems/
├── README.md                    # 本文档 - 完整工作流说明
├── model/                       # 磁性模型与结构
│   ├── magnetic_structure.md    # 磁性结构理论与建模
│   ├── sqs_generation.md        # SQS 结构生成方法
│   └── examples/                # 示例结构文件
└── calculation/                 # 计算设置与分析
    ├── noncollinear_setup.md    # 非共线磁性计算设置
    ├── dftu_parameters.md       # DFT+U 参数选择
    ├── convergence_tips.md      # 收敛技巧
    └── analysis_scripts/        # 分析脚本
```

## 🎯 项目背景：PdCrO2 钒掺杂

**研究目标：** 研究钒掺杂对 PdCrO2 磁性三角晶格材料磁学性质的影响

**挑战：**
- 复杂的 6 层磁性结构（18 个子晶格）
- 非共线磁矩（3D 磁矩向量）
- 自旋轨道耦合（SOC）效应
- 强关联电子（需要 DFT+U 修正）

## 🔄 完整工作流程

### 步骤 1：结构准备

#### 1.1 获取基础结构

```bash
# 从 Materials Project 或实验数据获取 PdCrO2 晶体结构
# 或使用 VESTA/Pymatgen 构建
```

#### 1.2 创建掺杂结构（SQS 方法）

使用 ICET 生成特殊准随机结构（Special Quasirandom Structure）：

```python
# 参考：sqs_generation/icet.ipynb
from icet import ClusterSpace, StructureContainer
from ase.io import read

# 读取原始结构
prim = read('POSCAR')

# 定义掺杂配置
# Cr -> V 掺杂，浓度 30%, 60%, 90%
```

**输出文件：**
- `POSCAR_SQS_V30.vasp` - 30% V 掺杂
- `POSCAR_SQS_V60.vasp` - 60% V 掺杂
- `POSCAR_SQS_V90.vasp` - 90% V 掺杂

详细说明见：[model/sqs_generation.md](model/sqs_generation.md)

---

### 步骤 2：磁性结构设计

#### 2.1 确定磁性配置

基于 **Takatsu et al. (2014)** 实验结果，PdCrO2 具有：
- **6 层反铁磁结构**
- **120° 磁矩旋转** (三角晶格特征)
- **非共线磁序**

#### 2.2 计算初始磁矩

使用 Python 脚本根据磁性模型计算每个原子的初始磁矩：

```python
# 参考：magnetic_analysis/magmom.ipynb
# 输出：3D 磁矩向量 (mx, my, mz)

# 示例输出：
# Cr atom 1: mx=0.584587, my=0.351255, mz=-0.731354
# Cr atom 2: mx=-0.697972, my=-0.674024, mz=-0.241922
```

**关键参数：**
- `α_n`: 层内旋转角
- `φ_n`: 层间相位差
- `γ_n`: 倾斜角
- `ξ_n`: 磁矩大小

详细说明见：[model/magnetic_structure.md](model/magnetic_structure.md)

---

### 步骤 3：DFT 计算设置

#### 3.1 自洽场（SCF）计算

**INCAR 关键设置：**

```bash
# ========== 基本精度 ==========
ALGO = Normal
PREC = Normal
EDIFF = 1e-6          # 收敛判据
ENCUT = 400           # 截断能
NELM = 500            # 最大电子步数

# ========== 电荷密度 ==========
ISTART = 0            # 从头开始
ICHARG = 2            # 叠加原子电荷
LCHARG = True         # 输出 CHGCAR
LWAVE = True          # 输出 WAVECAR

# ========== 非共线磁性 + SOC ==========
LSORBIT = True        # 开启自旋轨道耦合
LNONCOLLINEAR = True  # 非共线磁性（自动开启）

# 初始磁矩（3D 向量）
MAGMOM = 108*0.0 \      # 108 个 Pd 原子（非磁性）
  0.584587  0.351255 -0.731354 \   # Cr1: 120° 旋转
 -0.697972 -0.674024 -0.241922 \   # Cr2
  0.250611  0.150583  0.956305 \   # Cr3
  ... (共 18 个 Cr 磁性原子)
  54*0.0                # 54 个 O 原子（非磁性）

GGA_COMPAT = .FALSE.  # 兼容性设置

# ========== DFT+U 修正 ==========
LDAU = True           # 启用 DFT+U
LDAUTYPE = 2          # Dudarev 形式
LDAUL = -1 2 2 -1     # Pd=off, Cr=d, V=d, O=off
LDAUU = 0 3 4 0       # U: Cr=3eV, V=4eV
LDAUJ = 0 0 0.9 0     # J 值（Dudarev 中不使用）
LMAXMIX = 6           # 混合参数

# ========== 收敛加速 ==========
BMIX = 3              # 混合参数
AMIN = 0.01           # 最小混合参数
NELMDL = -5           # 跳过前 5 步对角化
ISMEAR = 0            # Gaussian smearing
SIGMA = 0.05          # 展宽宽度

# ========== 并行化（Perlmutter） ==========
NPAR = 32             # 并行组数
KPAR = 2              # k-point 并行
NCORE = 8             # 每组核心数
```

**KPOINTS：**

```bash
Automatic
0
Gamma
4 4 2    # 根据体系调整
```

**提交作业：**

```bash
sbatch submit_vasp.slurm
```

详细说明见：[calculation/noncollinear_setup.md](calculation/noncollinear_setup.md)

---

#### 3.2 态密度（DOS）计算

SCF 收敛后计算态密度：

```bash
# INCAR 修改
ISTART = 1            # 读取 WAVECAR
ICHARG = 11           # 读取 CHGCAR
LORBIT = 11           # 输出投影 DOS
NSW = 0               # 静态计算
IBRION = -1

# 增加 k-point 密度
# KPOINTS: 8 8 4
```

---

#### 3.3 能带结构计算

```bash
# 生成高对称点路径
# 对于三角晶格：Γ-M-K-Γ-A-L-H-A

ISTART = 1
ICHARG = 11
LORBIT = 11
LWAVE = False
```

详细说明见：[calculation/dftu_parameters.md](calculation/dftu_parameters.md)

---

### 步骤 4：计算监控与收敛检查

#### 4.1 监控作业状态

```bash
# 查看队列
squeue -u $USER

# 查看 OSZICAR（收敛情况）
tail -f OSZICAR

# 检查磁矩演化
grep "magnetization (x)" OSZICAR
```

#### 4.2 收敛判据

- **能量收敛：** `ΔE < 1e-6 eV`
- **磁矩稳定：** 磁矩在最后几步变化 < 0.01 μB
- **力收敛：** （如果做结构优化）`F < 0.01 eV/Å`

**常见问题与解决：**

| 问题 | 原因 | 解决方法 |
|------|------|----------|
| 磁矩振荡 | 初始磁矩不合理 | 调整 MAGMOM，增大 BMIX |
| 不收敛 | 电荷密度混合问题 | 减小 BMIX，使用 ALGO = All |
| 磁矩归零 | 对称性破缺不够 | ISYM = 0，增大初始磁矩 |

详细说明见：[calculation/convergence_tips.md](calculation/convergence_tips.md)

---

### 步骤 5：结果分析

#### 5.1 提取磁性性质

```bash
# 总磁矩
grep "number of electron" OUTCAR | tail -1

# 每个原子磁矩
grep "magnetization (x)" OUTCAR -A [原子数+4]
```

**使用 Python 分析：**

```python
from pymatgen.io.vasp import Outcar

# 读取 OUTCAR
outcar = Outcar("OUTCAR")

# 获取磁矩
magnetization = outcar.magnetization
print(f"Total magnetic moment: {sum(magnetization):.3f} μB")

# 绘制磁矩分布
import matplotlib.pyplot as plt
import numpy as np

mag_magnitudes = [np.linalg.norm(m) for m in magnetization]
plt.bar(range(len(mag_magnitudes)), mag_magnitudes)
plt.xlabel('Atom Index')
plt.ylabel('Magnetic Moment (μB)')
plt.title('Magnetic Moment Distribution')
plt.show()
```

#### 5.2 VESTA 可视化

将磁性结构导出为 MCIF 格式在 VESTA 中可视化：

```python
# 生成 vesta.mcif 文件
# 可以在 VESTA 中查看 3D 磁矩箭头
```

#### 5.3 态密度分析

```python
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter

vasprun = Vasprun("vasprun.xml")
dos = vasprun.complete_dos

# 绘制总 DOS
plotter = DosPlotter()
plotter.add_dos("Total", dos)
plotter.show()

# 分自旋 DOS
dos_up = dos.densities[Spin.up]
dos_down = dos.densities[Spin.down]
```

---

## 📊 完整计算流程图

```
┌─────────────────────────────────────────────────────────────┐
│  1. 结构准备                                                  │
│     • 基础晶体结构（Materials Project / 实验）                │
│     • SQS 掺杂结构生成（ICET）                                │
└────────────────┬────────────────────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────────────────────┐
│  2. 磁性模型设计                                              │
│     • 确定磁性配置（6 层，18 子晶格）                         │
│     • 计算初始 3D 磁矩（基于 Takatsu Model 4）                │
│     • 生成 MAGMOM 参数                                        │
└────────────────┬────────────────────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────────────────────┐
│  3. SCF 计算（Non-collinear + SOC + DFT+U）                   │
│     • LSORBIT = True                                          │
│     • MAGMOM = 3D 向量                                        │
│     • LDAU: Cr(U=3), V(U=4)                                   │
│     • 收敛判据：1e-6 eV                                        │
└────────────────┬────────────────────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────────────────────┐
│  4. 后处理计算                                                │
│     • DOS 计算（LORBIT=11）                                   │
│     • Band 计算（高对称点路径）                               │
│     • 电荷密度分析（CHGCAR）                                  │
└────────────────┬────────────────────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────────────────────┐
│  5. 结果分析与可视化                                          │
│     • 磁矩提取与分析（pymatgen）                              │
│     • VESTA 磁性结构可视化（MCIF）                            │
│     • 电子结构分析（vaspvis）                                 │
│     • 磁性相图构建                                             │
└─────────────────────────────────────────────────────────────┘
```

---

## 🔧 关键技术要点

### 1. 非共线磁性计算

- **自动开启：** 设置 `LSORBIT = True` 后，VASP 自动启用非共线模式
- **磁矩格式：** `MAGMOM = mx1 my1 mz1 mx2 my2 mz2 ...`
- **对称性：** 通常需要 `ISYM = 0` 或 `ISYM = -1` 破坏对称性

### 2. DFT+U 参数选择

| 元素 | U (eV) | J (eV) | 参考 |
|------|--------|--------|------|
| Cr   | 3.0-4.0| 0.9    | 文献 |
| V    | 4.0    | 0.9    | 经验 |
| Mn   | 4.0    | 1.0    | 标准 |

**参数调试：**
- 先用文献值
- 对比实验能隙/磁矩
- 扫描 U 值优化

详细说明见：[calculation/dftu_parameters.md](calculation/dftu_parameters.md)

### 3. 收敛加速技巧

**问题：** 磁性体系收敛慢

**解决方案：**
```bash
ALGO = All            # 尝试更稳定算法
BMIX = 5              # 增大混合参数
AMIX = 0.2
NELMDL = -10          # 跳过更多初始步骤
MAXMIX = 80           # 增加混合历史
```

详细说明见：[calculation/convergence_tips.md](calculation/convergence_tips.md)

---

## 📚 参考项目文件

本工作流基于以下项目：

```
PdCrO2_Doping_Magnetic_calculation/
├── magnetic_analysis/
│   └── magmom.ipynb              # 磁矩计算脚本
├── sqs_generation/
│   └── icet.ipynb                # SQS 生成脚本
└── dft_calculations/
    ├── nodoping/scf_noncollinear/
    │   └── INCAR                 # 非共线 SCF 设置
    └── doping/PdCrO2_V30%_U_3/
        └── scf/INCAR             # 掺杂体系设置
```

---

## 🎓 学习路径

### 初学者
1. 阅读 [model/magnetic_structure.md](model/magnetic_structure.md) 理解磁性模型
2. 学习基本的 INCAR 参数设置
3. 运行简单的铁磁/反铁磁计算

### 进阶
1. 掌握非共线磁性计算
2. 学习 DFT+U 参数调优
3. 使用 SQS 方法研究掺杂体系

### 专家
1. 复杂磁性结构建模（18 子晶格）
2. 磁性相图构建
3. 高通量磁性材料筛选

---

## 🛠️ 实用脚本

### 快速提取磁矩

```bash
#!/bin/bash
# extract_magmom.sh

echo "Total magnetization:"
grep "number of electron" OUTCAR | tail -1

echo -e "\nPer-atom magnetization:"
python << EOF
from pymatgen.io.vasp import Outcar
import numpy as np

out = Outcar("OUTCAR")
for i, mag in enumerate(out.magnetization):
    mag_size = np.linalg.norm(mag)
    if mag_size > 0.1:  # 只显示磁性原子
        print(f"Atom {i:3d}: |M| = {mag_size:6.3f} μB, "
              f"({mag[0]:7.3f}, {mag[1]:7.3f}, {mag[2]:7.3f})")
EOF
```

### 批量检查收敛

```bash
#!/bin/bash
# check_convergence.sh

for dir in calc_*/; do
    echo "=== $dir ==="
    cd "$dir"
    tail -1 OSZICAR | awk '{print "Energy: " $3 " eV"}'
    grep "magnetization (x)" OSZICAR | tail -1
    cd ..
done
```

---

## 📖 扩展阅读

- [VASP Wiki: Non-collinear Calculations](https://www.vasp.at/wiki/index.php/Non-collinear_calculations)
- [VASP Wiki: DFT+U](https://www.vasp.at/wiki/index.php/LDAU)
- Takatsu et al., Phys. Rev. B 89, 104408 (2014)
- ICET Documentation: https://icet.materialsmodeling.org/

---

## 💡 提示与最佳实践

1. **始终备份初始结构和 INCAR**
2. **记录所有参数选择的理由**
3. **对比文献中的实验数据验证**
4. **多次独立计算确认结果可靠性**
5. **使用版本控制（Git）管理计算**

---

## 🙏 致谢

本工作流基于 PdCrO2 钒掺杂磁性计算项目总结而成。感谢 Takatsu et al. 提供的实验磁性结构数据。
