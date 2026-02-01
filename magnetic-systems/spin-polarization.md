# 自旋极化计算

## 理论背景

自旋极化计算用于处理具有磁性的材料体系。VASP 通过 ISPIN 标签控制是否考虑自旋自由度。

## 基本参数设置

### ISPIN 标签

```bash
ISPIN = 1    # 非自旋极化（默认）
ISPIN = 2    # 自旋极化计算
```

### MAGMOM 标签

定义每个原子的初始磁矩（单位：μB）

```bash
# 示例1：4个Fe原子（5 μB）+ 4个O原子（0 μB）
MAGMOM = 4*5.0 4*0.0

# 示例2：详细指定每个原子
MAGMOM = 5.0 5.0 -5.0 -5.0 0.0 0.0 0.0 0.0

# 示例3：大体系简化写法
MAGMOM = 100*0.0    # 100个非磁性原子
```

## 计算流程

### 1. 结构优化（带自旋极化）

```bash
# INCAR 关键设置
ISPIN = 2
MAGMOM = ...        # 根据体系设置
IBRION = 2          # 优化算法
NSW = 100           # 最大优化步数
EDIFFG = -0.01      # 收敛标准
```

### 2. 静态自洽计算

```bash
ISPIN = 2
MAGMOM = ...
IBRION = -1         # 静态计算
NSW = 0
LORBIT = 11         # 输出投影态密度
```

## 实战技巧

### 💡 磁矩初值设置

1. **查阅文献**：参考已发表的相似体系
2. **原子价电子**：通常设为未配对电子数
3. **尝试不同初值**：避免陷入局域极小

**常见元素初始磁矩：**
- Fe: 4-5 μB
- Co: 3-4 μB
- Ni: 2-3 μB
- Mn: 4-5 μB
- O, C, N: 0 μB（一般情况）

### 💡 收敛检查

查看 OSZICAR 文件中的磁矩变化：

```bash
grep "magnetization" OSZICAR
```

最终磁矩在 OUTCAR 中：

```bash
grep "number of electron" OUTCAR | tail -1
```

### 💡 反铁磁计算

对于反铁磁体系，需要给不同子晶格相反的初始磁矩：

```bash
# 例：NiO 反铁磁
# Ni-up, Ni-down, O, O
MAGMOM = 5.0 -5.0 0.0 0.0
```

## 常见陷阱

### 🐛 磁矩收敛到 0

**原因：**
- 初始磁矩设置不合理
- 体系本身非磁性
- 对称性破缺不够

**解决方法：**
```bash
# 增大初始磁矩
MAGMOM = 8*10.0  # 尝试更大的初值

# 破坏对称性
ISYM = 0         # 关闭对称性
```

### 🐛 不同磁态能量接近

**问题：** 铁磁与反铁磁能量差异很小

**解决方法：**
- 增加 k-point 采样
- 使用更精确的交换关联泛函（如 HSE06）
- 考虑 DFT+U 修正

```bash
LDAU = .TRUE.
LDAUTYPE = 2
LDAUL = 2 -1        # Fe的d轨道，O的p轨道
LDAUU = 4.0 0.0     # U值
LDAUJ = 0.0 0.0
```

## 输出分析

### 查看最终磁矩

```bash
# 总磁矩
grep "number of electron" OUTCAR | tail -1

# 每个原子的磁矩
grep "magnetization (x)" OUTCAR -A [原子数+4]
```

### 绘制磁矩演化

```python
import numpy as np
import matplotlib.pyplot as plt

# 从 OSZICAR 提取磁矩数据
data = []
with open('OSZICAR', 'r') as f:
    for line in f:
        if 'mag=' in line:
            mag = float(line.split('mag=')[1].split()[0])
            data.append(mag)

plt.plot(data)
plt.xlabel('SCF Step')
plt.ylabel('Total Magnetization (μB)')
plt.title('Magnetization Convergence')
plt.show()
```

## 参考资源

- [VASP Wiki: ISPIN](https://www.vasp.at/wiki/index.php/ISPIN)
- [VASP Wiki: MAGMOM](https://www.vasp.at/wiki/index.php/MAGMOM)
- [Magnetism in DFT](https://www.vasp.at/wiki/index.php/Magnetism)
