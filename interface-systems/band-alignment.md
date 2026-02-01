# 界面能带对齐（Band Alignment）

## 理论背景

异质界面的能带对齐决定了电荷传输特性，是理解器件性能的关键。主要关注：
- **VBO (Valence Band Offset)** - 价带偏移
- **CBO (Conduction Band Offset)** - 导带偏移

## 计算方法

### 方法一：宏观平均静电势法

这是最常用的方法，通过对齐界面两侧的平均静电势来确定能带偏移。

#### 步骤 1：计算块体材料

分别计算两种材料的块体结构：

```bash
# INCAR - 块体 A
SYSTEM = Bulk A
ISTART = 0
ICHARG = 2
LVHAR = .TRUE.      # 输出 LOCPOT（Hartree势）
LVTOT = .TRUE.      # 输出总势能
LORBIT = 11         # 计算投影态密度
```

记录：
- VBM (价带顶)
- CBM (导带底)
- 平均静电势 V̄_bulk

#### 步骤 2：计算异质界面

建立界面超胞并计算：

```bash
# INCAR - 界面
SYSTEM = Interface A/B
LVHAR = .TRUE.
LVTOT = .TRUE.
LCHARGE = .TRUE.    # 输出电荷密度
```

#### 步骤 3：提取静电势

使用脚本计算沿界面法向的宏观平均静电势：

```python
import numpy as np
from pymatgen.io.vasp import Locpot

# 读取 LOCPOT
locpot = Locpot.from_file("LOCPOT")
avg_pot = locpot.get_average_along_axis(2)  # z方向平均

# 绘制势能曲线
import matplotlib.pyplot as plt
plt.plot(avg_pot)
plt.xlabel('z (Å)')
plt.ylabel('Electrostatic Potential (eV)')
plt.axhline(y=V_bulk_A, color='r', linestyle='--', label='Bulk A')
plt.axhline(y=V_bulk_B, color='b', linestyle='--', label='Bulk B')
plt.legend()
plt.show()
```

#### 步骤 4：计算能带偏移

```
ΔE_v = (VBM_A - V̄_A_interface) - (VBM_B - V̄_B_interface)
ΔE_c = (CBM_A - V̄_A_interface) - (CBM_B - V̄_B_interface)
```

或者：
```
ΔE_c = ΔE_v + (E_g^A - E_g^B)
```

## 实战技巧

### 💡 界面模型构建

1. **真空层厚度**：至少 15 Å，确保两侧势能趋于平台
2. **界面层数**：每侧至少 3-5 个原子层
3. **晶格匹配**：应变 < 5%

```python
# 使用 OgreInterface 或手动构建
# 确保对称性以避免偶极矩
```

### 💡 宏观平均处理

为消除原子级涨落，需要做宏观平均：

```python
def macroscopic_average(potential, window=5):
    """
    对势能做宏观平均
    window: 平均窗口大小（单位：Å）
    """
    avg = np.convolve(potential, np.ones(window)/window, mode='same')
    return avg
```

### 💡 VBM/CBM 确定

从 DOS 和能带结构确定：

```bash
# 提取费米能级
E_f=$(grep "E-fermi" OUTCAR | awk '{print $3}')

# 或从 DOSCAR 分析
# VBM = 最高占据态
# CBM = 最低未占据态
```

## 常见陷阱

### 🐛 偶极矩问题

**问题：** 不对称界面产生宏观偶极矩，导致势能不收敛

**解决方法：**
```bash
LDIPOL = .TRUE.     # 修正偶极矩
IDIPOL = 3          # 偶极方向（1=x, 2=y, 3=z）
```

或构建对称界面（两端相同终止面）

### 🐛 真空层不足

**问题：** 势能未趋于平台，无法确定参考能级

**解决方法：**
- 增加真空层至 20+ Å
- 检查 LOCPOT 势能曲线

### 🐛 k-point 不足

**问题：** VBM/CBM 位置不准确

**解决方法：**
```bash
# 使用更密的 k-point 网格
# 特别是界面法向（通常是 z 方向）
KPOINTS
Automatic
0
Gamma
8 8 1    # xy 方向密集，z 方向稀疏
```

## 完整工作流示例

```bash
# 1. 计算块体 A
cd bulk_A
vasp_std > vasp.out
E_vbm_A=$(grep "VBM" analysis.log)
V_bulk_A=$(python get_avg_potential.py)

# 2. 计算块体 B
cd ../bulk_B
vasp_std > vasp.out
E_vbm_B=$(grep "VBM" analysis.log)
V_bulk_B=$(python get_avg_potential.py)

# 3. 计算界面
cd ../interface
vasp_std > vasp.out

# 4. 分析界面势能
python << EOF
from analyze_interface import get_band_offset
VBO, CBO = get_band_offset(
    E_vbm_A, V_bulk_A,
    E_vbm_B, V_bulk_B,
    'LOCPOT'
)
print(f"VBO = {VBO:.3f} eV")
print(f"CBO = {CBO:.3f} eV")
EOF
```

## 可视化脚本

```python
# plot_band_alignment.py
import numpy as np
import matplotlib.pyplot as plt

def plot_alignment(VBO, CBO, Eg_A, Eg_B):
    fig, ax = plt.subplots(figsize=(8, 6))

    # 材料 A
    ax.plot([0, 1], [0, 0], 'b-', linewidth=3, label='A VBM')
    ax.plot([0, 1], [Eg_A, Eg_A], 'b--', linewidth=3, label='A CBM')

    # 材料 B
    ax.plot([2, 3], [VBO, VBO], 'r-', linewidth=3, label='B VBM')
    ax.plot([2, 3], [VBO+Eg_B, VBO+Eg_B], 'r--', linewidth=3, label='B CBM')

    # 标注
    ax.text(0.5, -0.3, f'VBO = {VBO:.2f} eV', ha='center')
    ax.text(0.5, Eg_A+0.3, f'CBO = {CBO:.2f} eV', ha='center')

    ax.set_ylabel('Energy (eV)')
    ax.set_xlim(-0.5, 3.5)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('band_alignment.png', dpi=300)
    plt.show()

# 使用示例
plot_alignment(VBO=0.5, CBO=0.8, Eg_A=3.0, Eg_B=2.5)
```

## 参考文献

1. Van de Walle & Martin, *Phys. Rev. B* 35, 8154 (1987) - 经典方法
2. Schleife et al., *Appl. Phys. Lett.* 94, 012104 (2009) - 实际应用
3. [VASP Wiki: Band Alignment](https://www.vasp.at/wiki/index.php/Band_alignment)

## 相关脚本

- `get_avg_potential.py` - 提取平均静电势
- `analyze_interface.py` - 自动化分析
- `plot_band_alignment.py` - 可视化
