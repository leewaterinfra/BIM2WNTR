# BIM2WNTR

**Automated Conversion of Building Information Models to Hydraulic Network Models for Premise Plumbing Analysis**

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

BIM2WNTR is an open-source Python tool that automates the conversion of Industry Foundation Classes (IFC) building information models to Water Network Tool for Resilience (WNTR) hydraulic models. This enables hydraulic analysis of premise plumbing systems directly from BIM data, supporting digital twin development for building water systems.

**Supplementary material for:**
> Lee, J. (2025). "Automated Conversion of Building Information Models to Hydraulic Network Models for Premise Plumbing Analysis." *ASCE Journal of Computing in Civil Engineering*. (Technical Note)

## Features

- **IFC Data Extraction**: Parse plumbing elements (`IfcFlowSegment`, `IfcFlowFitting`, `IfcFlowTerminal`)
- **Topology Reconstruction**: Use port-based connectivity (`IfcRelConnectsPorts`) to preserve network structure
- **Real 3D Coordinates**: Extract actual geometry coordinates from IFC representations
- **System Filtering**: Isolate Domestic Cold Water (DCW) systems
- **Fixture-Based Demands**: Assign demands based on IPC/UPC fixture flow rates
- **Diversity Factor**: Apply Hunter's curve methodology for realistic demand scenarios
- **EPANET Export**: Generate simulation-ready INP files with proper coordinates

## Installation

### Requirements

- Python 3.8 or higher
- Dependencies listed in `requirements.txt`

### Install Dependencies

```bash
pip install -r requirements.txt
```

Or install individually:

```bash
pip install ifcopenshell wntr networkx matplotlib seaborn numpy
```

### Note on IfcOpenShell

IfcOpenShell may require additional installation steps depending on your platform:

```bash
# Using conda (recommended)
conda install -c conda-forge ifcopenshell

# Using pip (may need pre-built wheels)
pip install ifcopenshell
```

## Usage

### Command Line

```bash
python bim2wntr.py <path_to_ifc_file> [output_directory]
```

**Example:**

```bash
python bim2wntr.py building.ifc ./output
```

### Python API

```python
from bim2wntr import run_pipeline

# Run conversion pipeline
models = run_pipeline('building.ifc', './output')

# Access WNTR model
wn = models[0]['wn']
print(f"Junctions: {wn.num_junctions}")
print(f"Pipes: {wn.num_pipes}")
```

## Output Files

| File | Description |
|------|-------------|
| `DCW_*.inp` | EPANET INP file with real coordinates |
| `ifc_extracted.json` | Raw extracted IFC data |
| `simulation_results.json` | Hydraulic simulation summary |
| `figure_network_layout.png/pdf` | Network layout visualization |
| `figure_network_analysis.png/pdf` | Network statistics (diameter, length distributions) |
| `figure_hydraulic_results.png/pdf` | Simulation results (pressure, velocity, flow) |

## Methodology

### 1. IFC Data Extraction

The tool extracts plumbing elements from IFC files:

| IFC Entity | WNTR Mapping |
|------------|--------------|
| `IfcFlowSegment` | Pipes |
| `IfcFlowFitting` | Junctions (tees, elbows) |
| `IfcFlowTerminal` | Junctions with demand (fixtures) |

### 2. Topology Reconstruction

Network topology is reconstructed using port-based connectivity:

```
IfcDistributionPort → IfcRelConnectsPorts → Connected Elements
```

This preserves the actual piping connections rather than relying on geometric proximity.

### 3. Coordinate Extraction

Real 3D coordinates are extracted from `IfcExtrudedAreaSolid` geometry representations using `ifcopenshell.geom`. Coordinates are transformed for EPANET visualization:

- X_EPANET = (IFC_X + IFC_Y) × scale (isometric view)
- Y_EPANET = IFC_Z × scale (elevation)

### 4. Demand Allocation

Fixture demands are assigned based on IPC/UPC flow rates:

| Fixture Type | Flow Rate (GPM) | Fixture Units |
|--------------|-----------------|---------------|
| Lavatory | 2.0 | 1.0 |
| Shower | 2.5 | 2.0 |
| Kitchen Sink | 2.5 | 1.5 |
| Bathtub | 4.0 | 2.0 |
| Hose Bib | 2.5 | 2.5 |

### 5. Diversity Factor

A diversity factor based on Hunter's curve is applied to account for the probability that not all fixtures operate simultaneously:

| Number of Fixtures | Diversity Factor |
|--------------------|------------------|
| 1 | 1.00 |
| 2-5 | 0.80 |
| 6-10 | 0.60 |
| 11-20 | 0.40 |
| 21-50 | 0.30 |
| 51-100 | 0.25 |
| >100 | 0.20 |

## IFC Requirements

For successful conversion, the IFC file should include:

- **MEP Connectivity**: Export with "Include MEP connectivity" enabled in Revit
- **Property Sets**: Standard property sets (`Pset_PipeSegmentPipeCommon`)
- **System Assignments**: Elements assigned to plumbing systems (DCW, DHW)

Tested IFC schemas: IFC2X3, IFC4

## Limitations

- Currently processes Domestic Cold Water (DCW) systems only
- Hot water systems require thermal modeling (future work)
- Depends on IFC export quality from BIM authoring tools
- Recirculation loops not automatically detected

## References

1. Lee, J., et al. (2023). *Premise Plumbing Modeling*. IWA Publishing.
2. Hunter, R.B. (1940). "Methods of Estimating Loads in Plumbing Systems." NBS Building Materials and Structures Report BMS65.
3. Klise, K.A., et al. (2017). "Water Network Tool for Resilience (WNTR) User Manual." Sandia National Laboratories.
4. buildingSMART International. *IFC Specifications*. https://www.buildingsmart.org/

## Citation

If you use this tool in your research, please cite:

```bibtex
@article{lee2025bim2wntr,
  title={Automated Conversion of Building Information Models to Hydraulic Network Models for Premise Plumbing Analysis},
  author={Lee, Juneseok},
  journal={ASCE Journal of Computing in Civil Engineering},
  year={2025},
  note={Technical Note}
}
```

## Author

**Dr. Juneseok Lee, PhD, PE, BC.WRE, EWRI Fellow, ASCE Fellow**  
Professor, Civil and Environmental Engineering  
Manhattan University  

R.A.I.S.E. – Resilience, AI, Infrastructure, Systems, and Education

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Water Network Tool for Resilience (WNTR) development team
- IfcOpenShell contributors
- buildingSMART International for IFC standards
