# OgreInterface Tutorial

A tutorial for creating and optimizing epitaxial interfaces using OgreInterface.

---

## Installation

### 1. Create Conda Environment

```bash
conda create -n ogre_py311 python=3.11
conda activate ogre_py311
```

### 2. Install Dependencies

```bash
pip install torch
pip install torch-dftd
pip install mace-torch
pip install cuequivariance-torch
pip install chgnet
```

### 3. Install OgreInterface

```bash
pip install "OgreInterface[notebook] @ git+https://github.com/caizefeng/OgreInterface.git"
```

> **Note:** Contact Zefeng Cai for repository access if needed.

### 4. Configure Jupyter Kernel

Create a named kernel for easy identification:

```bash
python -m ipykernel install --user --name ogre_py311 --display-name "Python 3.11 (ogre_py311)"
```

**Why?** Creating a conda environment doesn't automatically create a Jupyter kernel with the same name. Without this step, all environments appear as "Python 3" in Jupyter, making them indistinguishable.

### 5. Verify Installation

```bash
python -c "import OgreInterface; print('✓ OgreInterface installed')"
jupyter --version
```

---

## Launch Jupyter Notebook

```bash
conda activate ogre_py311
cd "/path/to/ogreinterface-tutorial"
jupyter lab
```

Access the notebook at: `http://localhost:8888/lab?token=xxxxxxxxxx`

**Important:** Select **"Python 3.11 (ogre_py311)"** as the kernel in Jupyter.

---

## Workflow Overview

OgreInterface predicts and optimizes epitaxial interface structures between two materials through four main steps:

```
Input CIF files → Lattice Matching → Interface Generation → Surface Matching
```

---

### 1. Lattice Matching

**Purpose:** Identify interface orientations with minimal strain and optimal area.

**Method:** Based on lattice parameters only, without considering atomic interactions.

**Input:**
- CIF files of two materials (A and B)
- Constraints: max Miller index, max strain, max area

**Output:**
A visualization plot where:
- **Circle size** = interface area (smaller is better)
- **Circle color** = strain magnitude (bluer/smaller is better)
- **Shaded region** = actual interface area

**Example:** For a perovskite-semiconductor interface, the smallest and bluest circle represents the optimal match (e.g., (210) perovskite surface + (221) PbS surface).

---

### 2. Interface Generation

**Purpose:** Visualize supercells and select configurations for optimization.

**Key Parameters:**
```python
hkl_A = [1, 0, 0]       # Miller indices for material A
hkl_B = [0, 1, 0]       # Miller indices for material B
strain_fraction = 0.5   # Strain partition coefficient
```

**Strain Partition:**
- `0.0` = all strain in material B (film) - typical for epitaxial growth
- `1.0` = all strain in material A (substrate)
- `0.5` = strain equally shared - typical for nanocrystals

---

### 3. Surface Matching and Ranking

**Purpose:** Optimize atomic configurations and calculate interface energies.

**Method:**
Lattice matching ignores atomic interactions. Surface matching:
1. **Enumerates surface terminations** for the material pair
2. **Combines terminations** to form different interface configurations
3. **Calculates adhesion energy** between the two materials
4. **Computes interface energy** for each configuration
5. **Ranks interfaces** by stability

**Optimization Algorithms:**
- **PSO** (Particle Swarm Optimization) - default
- **BO** (Bayesian Optimization)

**Key Parameters:**
```python
supercell_choice = 0              # 0 = best interface
minimum_slab_thickness = 18       # Minimum slab thickness (Å)
opt_method = "PSO"                # Optimization method
n_particles_PSO = 30              # More particles = more accurate but slower
filter_on_charge = True           # Filter by surface charge
mlip = None                       # Optional: use ML potential
```

---

## Workflow Diagram

```
┌────────────────────────────────────────┐
│  Input: CIF files of materials A & B   │
└────────────────────────────────────────┘
                 ↓
┌────────────────────────────────────────┐
│  Lattice Matching                      │
│  - Based on lattice parameters only    │
│  - Output: possible orientations       │
│  - Select: minimal strain + area       │
└────────────────────────────────────────┘
                 ↓
┌────────────────────────────────────────┐
│  Interface Generation                  │
│  - Visualize supercells                │
│  - Set strain partition                │
│  - Select configuration                │
└────────────────────────────────────────┘
                 ↓
┌────────────────────────────────────────┐
│  Surface Matching                      │
│  - Enumerate terminations              │
│  - Calculate adhesion/interface energy │
│  - Optimize atomic positions (PSO/BO)  │
│  - Output: ranked stable interfaces    │
└────────────────────────────────────────┘
```

---

## Usage

### Directory Structure

```
ogreinterface-tutorial/
├── README.md
├── tutorial.ipynb
└── Reference CIFs/
    ├── CsPbBr3_cubic.cif
    └── Pb4S3Br2.cif
```

### Running the Tutorial

1. **Launch Jupyter**
   ```bash
   conda activate ogre_py311
   cd "/path/to/ogreinterface-tutorial"
   jupyter lab
   ```

2. **Open `tutorial.ipynb`**

3. **Run cells sequentially:**
   - **Initialize:** Load OgreInterface modules
   - **Lattice Matching:** Analyze lattice compatibility
   - **Interface Generation:** Visualize supercells
   - **Surface Matching:** Optimize and rank interfaces

---

## References

Toso et al. (2025) - "Structure Prediction of Ionic Epitaxial Interfaces with Ogre Demonstrated for Colloidal Heterostructures"

---

**Credits:**
Derek Dardzinski (original author) | Zefeng Cai (maintainer)

---

## P.S. Using Conda ogre Environment as NERSC Jupyter Kernel

### Create Jupyter Kernel

```bash
conda activate ogre_py311
python -m ipykernel install --user --name ogre_py311 --display-name "ogre_py311"
```

### Create Kernel Helper Script

```bash
cat > ~/.local/share/jupyter/kernels/ogre_py311/kernel-helper.sh << 'EOF'
#!/bin/bash
# Load Python module (NERSC specific)
module load python

# Activate conda environment
conda activate ogre_py311

# Execute the kernel
exec "$@"
EOF

# Set executable permissions
chmod u+x ~/.local/share/jupyter/kernels/ogre_py311/kernel-helper.sh

# Fix line endings
sed -i 's/\r$//' ~/.local/share/jupyter/kernels/ogre_py311/kernel-helper.sh
```

### Update kernel.json

Modify `~/.local/share/jupyter/kernels/ogre_py311/kernel.json` to use **absolute path**:

```json
{
 "argv": [
  "/global/u1/c/YOUR_USERNAME/.local/share/jupyter/kernels/ogre_py311/kernel-helper.sh",
  "python",
  "-Xfrozen_modules=off",
  "-m",
  "ipykernel_launcher",
  "-f",
  "{connection_file}"
 ],
 "display_name": "ogre_py311",
 "language": "python",
 "metadata": {
  "debugger": true
 }
}
```

**Important**:
- Use **absolute path** (not `{resource_dir}`) to ensure NERSC JupyterHub can find the script
- Ensure Unix line endings (LF) not Windows (CRLF)
- Helper script must have executable permissions

### Kernel Directory Structure

```
~/.local/share/jupyter/kernels/ogre_py311/
├── kernel.json          # Kernel configuration file
├── kernel-helper.sh     # Helper startup script (executable)
├── logo-32x32.png       # Small icon
├── logo-64x64.png       # Medium icon
└── logo-svg.svg         # Vector icon
```

### Troubleshooting

**Kernel fails to start:**
```bash
# Check permissions
ls -l ~/.local/share/jupyter/kernels/ogre_py311/kernel-helper.sh

# Fix line endings
sed -i 's/\r$//' ~/.local/share/jupyter/kernels/ogre_py311/kernel-helper.sh

# Test helper script
~/.local/share/jupyter/kernels/ogre_py311/kernel-helper.sh python --version
```

For detailed instructions, see: https://docs.nersc.gov/services/jupyter/how-to-guides/
