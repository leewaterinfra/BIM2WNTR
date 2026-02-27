#!/usr/bin/env python3
"""
Example Usage of BIM2WNTR
=========================

This script demonstrates how to use the BIM2WNTR module for converting
IFC files to WNTR hydraulic models.

Usage:
    python example_usage.py
"""

from pathlib import Path

# Import BIM2WNTR components
from bim2wntr import (
    run_pipeline,
    IFCExtractor,
    TopologyAnalyzer,
    WNTRBuilder,
    HydraulicSimulator,
    Visualizer,
    Config,
    DiversityFactor,
    PipeDiameterConverter  # New in v1.1.0
)


def example_basic():
    """Basic usage: Run the complete pipeline."""
    print("=" * 60)
    print("EXAMPLE 1: Basic Pipeline Usage")
    print("=" * 60)
    
    # Specify input and output
    ifc_file = "your_building.ifc"  # Replace with your IFC file
    output_dir = "./output"
    
    # Check if file exists
    if not Path(ifc_file).exists():
        print(f"File not found: {ifc_file}")
        print("Please update the ifc_file path to your IFC file.")
        return
    
    # Run pipeline
    models = run_pipeline(ifc_file, output_dir)
    
    # Access results
    if models:
        model = models[0]
        wn = model['wn']
        print(f"\nModel created successfully!")
        print(f"  Junctions: {wn.num_junctions}")
        print(f"  Pipes: {wn.num_pipes}")
        print(f"  INP file: {model['inp_path']}")


def example_step_by_step():
    """Step-by-step usage: Manual control of each step."""
    print("=" * 60)
    print("EXAMPLE 2: Step-by-Step Usage")
    print("=" * 60)
    
    ifc_file = "your_building.ifc"  # Replace with your IFC file
    
    if not Path(ifc_file).exists():
        print(f"File not found: {ifc_file}")
        return
    
    # Step 1: Extract IFC data
    print("\nStep 1: Extracting IFC data...")
    extractor = IFCExtractor(ifc_file)
    data = extractor.extract_all()
    
    print(f"  Pipes: {len(data['pipes'])}")
    print(f"  Fittings: {len(data['fittings'])}")
    print(f"  Terminals: {len(data['terminals'])}")
    
    # Step 2: Analyze topology
    print("\nStep 2: Analyzing topology...")
    analyzer = TopologyAnalyzer(data, extractor)
    components = analyzer.get_components(dcw_only=True)
    
    print(f"  DCW components found: {len(components)}")
    
    # Step 3: Build WNTR model
    print("\nStep 3: Building WNTR model...")
    if components:
        comp = components[0]  # Use first DCW component
        builder = WNTRBuilder(
            analyzer.all_elements, 
            data['connections'],
            extractor.element_coords
        )
        wn, junction_map = builder.build_from_component(comp['elements'])
        
        print(f"  Junctions: {wn.num_junctions}")
        print(f"  Pipes: {wn.num_pipes}")
        
        # Step 4: Run simulation
        print("\nStep 4: Running hydraulic simulation...")
        simulator = HydraulicSimulator(
            wn, analyzer.all_elements, builder.junction_to_element
        )
        
        demand_info = simulator.add_fixture_demands(apply_diversity=True)
        print(f"  Active fixtures: {demand_info['num_fixtures']}")
        print(f"  Diversity factor: {demand_info['diversity_factor']:.2f}")
        
        if simulator.run_steady_state():
            results = simulator.get_results_summary()
            print(f"\n  Pressure range: {results['pressure']['min_psi']:.1f} - "
                  f"{results['pressure']['max_psi']:.1f} psi")
            print(f"  Max velocity: {results['velocity']['max_ms']:.3f} m/s")


def example_custom_config():
    """Example: Customizing configuration parameters."""
    print("=" * 60)
    print("EXAMPLE 3: Custom Configuration")
    print("=" * 60)
    
    # Modify configuration
    Config.DEFAULT_SUPPLY_HEAD = 40  # Higher pressure (40m instead of 30m)
    Config.ROUGHNESS_COPPER = 140    # Different roughness value
    
    # Add custom fixture type
    Config.FIXTURE_DEMANDS['custom_fixture'] = (0.000200, 3.0)  # 3.2 GPM, 3 FU
    
    print("Configuration updated:")
    print(f"  Supply head: {Config.DEFAULT_SUPPLY_HEAD} m")
    print(f"  Roughness: {Config.ROUGHNESS_COPPER}")
    print(f"  Custom fixture added: 'custom_fixture'")


def example_diversity_factor():
    """Example: Using diversity factor calculator (updated in v1.1.0)."""
    print("=" * 60)
    print("EXAMPLE 4: Diversity Factor Calculation (v1.1.0)")
    print("=" * 60)
    
    fixture_counts = [1, 5, 10, 20, 50, 100, 200]
    
    print("\nFixtures | Diversity Factor")
    print("-" * 30)
    for n in fixture_counts:
        df = DiversityFactor.calculate(n)
        print(f"    {n:4d} | {df:.2f}")
    
    print("\nNote: v1.1.0 updated diversity factors based on industry feedback")
    print("      (actual simultaneous use often 15-20% for residential)")
    
    # Hunter's curve estimation
    print("\nFixture Units | Peak Demand (GPM)")
    print("-" * 35)
    for fu in [10, 50, 100, 200, 500]:
        gpm = DiversityFactor.hunter_curve_gpm(fu)
        print(f"        {fu:4d} | {gpm:.1f}")


def example_diameter_conversion():
    """Example: OD to ID conversion using copper Type L specs (new in v1.1.0)."""
    print("=" * 60)
    print("EXAMPLE 5: Pipe Diameter Conversion (v1.1.0)")
    print("=" * 60)
    
    print("\nCopper Type L - Outside Diameter to Inside Diameter Conversion")
    print("-" * 65)
    print(f"{'Nominal':>10} | {'OD (mm)':>10} | {'Wall (mm)':>10} | {'ID (mm)':>10}")
    print("-" * 65)
    
    # Test common pipe sizes
    test_ods_mm = [15.875, 22.225, 28.575, 34.925, 53.975, 79.375]
    nominal_sizes = ['1/2"', '3/4"', '1"', '1-1/4"', '2"', '3"']
    
    for nominal, od_mm in zip(nominal_sizes, test_ods_mm):
        od_m = od_mm / 1000
        id_m = PipeDiameterConverter.od_to_id(od_m)
        wall_m = PipeDiameterConverter.get_wall_thickness(od_m)
        id_mm = id_m * 1000
        wall_mm = wall_m * 1000
        print(f"{nominal:>10} | {od_mm:>10.2f} | {wall_mm:>10.2f} | {id_mm:>10.2f}")
    
    print("\nNote: IFC geometry represents pipe OUTSIDE diameter.")
    print("      WNTR/EPANET hydraulic calculations use INSIDE diameter.")
    print("      v1.1.0 automatically converts ODâ†’ID using Type L copper specs.")


if __name__ == "__main__":
    # Run examples
    # Uncomment the example you want to run
    
    # example_basic()
    # example_step_by_step()
    example_custom_config()
    example_diversity_factor()
    example_diameter_conversion()  # New in v1.1.0
