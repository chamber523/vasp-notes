# VASP 计算笔记

系统整理的 VASP 第一性原理计算实践笔记，记录磁性体系、界面体系和数据可视化的经验总结。

## 📚 内容概览

本笔记集涵盖以下主题：

### 1. 磁性体系计算
- 自旋极化设置（ISPIN, MAGMOM）
- 磁矩初始化与收敛技巧
- 反铁磁与铁磁计算
- 磁各向异性能计算

### 2. 界面体系与能带对齐
- 异质结构建模
- Band Alignment 计算方法
- 静电势对齐技术
- 能带偏移（VBO/CBO）计算
- 界面能量计算

### 3. 数据可视化（vaspvis）
- 能带结构绘制
- 态密度（DOS/PDOS）可视化
- 电荷密度图
- 能带投影分析
- 自定义绘图技巧

### 4. 常见问题与优化
- 计算不收敛问题排查
- 参数优化策略
- 性能调优经验
- 错误信息解读

## 📖 使用说明

每个主题包含：
- ✅ **理论背景** - 简要的物理/化学背景
- ⚙️ **参数设置** - 关键 INCAR/KPOINTS 参数
- 💡 **实战技巧** - 基于实践的经验总结
- 🐛 **常见陷阱** - 需要注意的问题
- 📝 **示例代码** - 可复用的脚本和配置

## 🗂️ 目录结构

```
vasp-notes/
├── magnetic-systems/          # 磁性体系
│   ├── spin-polarization.md   # 自旋极化计算
│   ├── magmom-setup.md        # 磁矩设置
│   └── examples/              # 示例输入文件
├── interface-systems/         # 界面体系
│   ├── band-alignment.md      # 能带对齐
│   ├── electrostatic.md       # 静电势对齐
│   └── examples/              # 示例计算
├── visualization/             # 可视化
│   ├── vaspvis-guide.md       # vaspvis 使用指南
│   ├── band-structure.md      # 能带结构绘制
│   ├── dos-plotting.md        # 态密度绘图
│   └── scripts/               # 绘图脚本
├── troubleshooting/           # 问题排查
│   ├── convergence.md         # 收敛问题
│   └── common-errors.md       # 常见错误
└── templates/                 # 模板文件
    ├── INCAR_magnetic         # 磁性计算模板
    ├── INCAR_interface        # 界面计算模板
    └── plotting_templates/    # 绘图模板
```

## 🚀 快速开始

### 磁性体系计算示例

```bash
# 铁磁计算基本设置
ISPIN = 2                    # 开启自旋极化
MAGMOM = 4*5.0 4*0.0        # 初始磁矩设置
```

### 能带对齐计算示例

```bash
# 计算平均静电势
# 用于界面能带对齐分析
LVHAR = .TRUE.
LVTOT = .TRUE.
```

### vaspvis 可视化示例

```python
from vaspvis import Band

# 绘制能带结构
band = Band(folder='path/to/calculation')
band.plot_bands()
```

## 🔗 相关资源

- [VASP 官方文档](https://www.vasp.at/wiki/)
- [VASP Wiki](https://www.vasp.at/wiki/index.php/The_VASP_Manual)
- [vaspvis 文档](https://github.com/DerekDardzinski/vaspvis)
- [Materials Project](https://materialsproject.org/)

## 📝 贡献指南

这是个人学习笔记，欢迎提出建议和补充！

## ⚖️ 许可证

MIT License

## 🙏 致谢

感谢 VASP 开发团队和开源社区的贡献。
