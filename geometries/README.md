# Lattice Geometry Files

This directory contains geometry definitions for various lattice structures used in quantum simulations. Each `.geom` file defines a crystal structure with its symmetries, orbital positions, and hopping parameters.

## Directory Structure

geometries/
├── Basic Geometries
│ ├── 1d.geom # 1D chain
│ ├── 2Dinterface.geom # 2D interface system
│ ├── 2Dinterface_small.geom
│ ├── cubic.geom # 3D cubic lattice
│ ├── honeycomb.geom # 2D honeycomb lattice
│ ├── layers.geom # Layered structure
│ ├── lieb.geom # 2D Lieb lattice
│ ├── square.geom # 2D square lattice
│ ├── strip.geom # Quasi-1D strip geometry
│ └── triangle.geom # 2D triangular lattice
└── templates/ # Parameterized template files
├── chain_template # 1D chain with tunable parameters
├── cubic_template # 3D cubic lattice template
├── honeycomb_template # Honeycomb lattice template
├── input_template # Base template for new geometries
├── kagome_anisotropic_template
├── lieb_template
├── miyahara_flat_band_template
├── one_fifth_depleted_template
├── 9_16_depleted_template
├── square_anisotropic_template
├── square_template
└── triangle_template



## File Format

Each geometry file follows a standardized format with mandatory and optional sections:

### Required Sections

- `#NDIM`: Number of dimensions (1, 2, or 3)
- `#PRIM`: Primary lattice vectors (one per dimension)
- `#SUPER`: Supercell size in each dimension
- `#ORB`: Orbital positions (format: `label x y z [#comment]`)
- `#HAMILT`: Hopping parameters and interactions
  - Format: `orb1 orb2 dx dy dz t↑ t↓ U`
  - Where:
    - `orb1, orb2`: Orbital indices
    - `dx, dy, dz`: Hopping vector
    - `t↑, t↓`: Hopping parameters for up/down spins
    - `U`: On-site interaction
- `#SYMM`: Symmetry operations
- `#END`: End of file marker

### Optional Sections

- `#PHASE`: Phase factors for symmetry operations
- `#BONDS`: Defines bonds between sites for visualization and analysis
- `#PAIR`: Pairing symmetries for superconducting calculations

## Template Files

Template files contain parameterized versions of the geometries with placeholders:
- `{nx}, {ny}, {nz}`: Supercell dimensions
- `{t}`: Primary hopping parameter
- `{tp}`: Secondary hopping parameter
- `{u}`: On-site interaction strength

### Notable Templates

#### Kagome Anisotropic Template
- 2D kagome lattice with 3 sites per unit cell
- Tunable nearest-neighbor (t) and next-nearest-neighbor (tp) hoppings
- Supports anisotropic interactions

#### 9/16 Depleted Template
- 2D depleted square lattice
- 7 sites per unit cell
- Used for studying frustrated magnetism

#### Miyahara Flat Band Template
- Specialized geometry for flat band physics
- Supports the study of strongly correlated states

## Symmetry Operations

The `#SYMM` section defines the symmetry operations:
- `d`: Mirror/reflection
- `c2`: 2-fold rotation
- `c3`: 3-fold rotation
- `c4`: 4-fold rotation
- `s4`, `s6`: Improper rotations
- `i`: Inversion

Format: `operation x0 y0 z0 nx ny nz`
- `x0,y0,z0`: Center of rotation/reflection
- `nx,ny,nz`: Normal vector (for reflections) or rotation axis

## Usage Examples

### Basic Square Lattice

#NDIM
2
#PRIM
1.0 0.0
0.0 1.0
#SUPER
4 4

### Creating New Geometries

1. Start with `input_template`
2. Define dimensionality and lattice vectors
3. Add orbital positions
4. Define hopping parameters
5. Add symmetry operations
6. (Optional) Add bonds and pairing definitions

## Common Applications

- Single-band Hubbard models
- Multi-orbital systems
- Frustrated magnetism
- Topological phases
- Superconductivity studies
- Interface physics
