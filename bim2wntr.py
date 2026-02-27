#!/usr/bin/env python3
"""
BIM2WNTR: Automated Conversion of IFC Building Information Models 
           to WNTR Hydraulic Network Models for Premise Plumbing Analysis
================================================================================

Author: Dr. Juneseok Lee, Manhattan University
        R.A.I.S.E. - Resilience, AI, Infrastructure, Systems, and Education

Purpose: Supplementary code for Technical Note submitted to 
         ASCE Journal of Computing in Civil Engineering

Version: 1.1.0

Description:
    This script automates the conversion of Industry Foundation Classes (IFC) 
    files to Water Network Tool for Resilience (WNTR) hydraulic models. It 
    extracts plumbing system data from BIM models, reconstructs network topology 
    using port-based connectivity, and generates simulation-ready EPANET INP files.

Features:
    - IFC data extraction (pipes, fittings, terminals)
    - Port-based topology reconstruction (IfcRelConnectsPorts)
    - Real 3D coordinate extraction from IFC geometry
    - System filtering (DCW - Domestic Cold Water)
    - Fixture-based demand allocation (IPC/UPC flow rates)
    - Diversity factor application (Hunter's curve)
    - Publication-quality visualization
    - OD to ID conversion using copper Type L wall thickness (v1.1.0)

Dependencies:
    pip install ifcopenshell wntr networkx matplotlib seaborn numpy

Usage:
    python bim2wntr.py <path_to_ifc_file> [output_directory]

Example:
    python bim2wntr.py building.ifc ./output

References:
    - Lee et al. (2023). Premise Plumbing Modeling. IWA Publishing.
    - Hunter (1940). Methods of Estimating Loads in Plumbing Systems. NBS.
    - Klise et al. (2017). Water Network Tool for Resilience (WNTR).

License: MIT License
================================================================================
"""

import ifcopenshell
import ifcopenshell.geom
import wntr
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
import sys
import warnings

warnings.filterwarnings('ignore')

__version__ = "1.1.0"
__author__ = "Dr. Juneseok Lee"


# =============================================================================
# CONFIGURATION
# =============================================================================
class Config:
    """
    Configuration parameters for IFC-to-WNTR conversion.
    
    Attributes:
        ROUGHNESS_COPPER: Hazen-Williams roughness coefficient for copper pipes
        DEFAULT_SUPPLY_HEAD: Default reservoir head in meters
        MIN_PIPE_LENGTH: Minimum pipe length in meters
        DEFAULT_DIAMETER: Default pipe diameter in meters
        TARGET_SYSTEMS: Keywords to identify cold water systems
        FIXTURE_DEMANDS: Dictionary of fixture demands and fixture units
    """
    
    # Hydraulic parameters
    ROUGHNESS_COPPER = 130          # Hazen-Williams C for copper
    DEFAULT_SUPPLY_HEAD = 30        # meters (~44 psi)
    MIN_PIPE_LENGTH = 0.1           # meters
    DEFAULT_DIAMETER = 0.0254       # meters (1 inch)
    
    # System filter - only process cold water systems
    TARGET_SYSTEMS = ['DCW', 'COLD', 'CW', 'DOMESTIC COLD']
    
    # Copper Type L Wall Thickness Lookup Table
    # Format: nominal_size_inches -> (OD_mm, wall_thickness_mm, ID_mm)
    # Reference: ASTM B88 Standard Specification for Seamless Copper Water Tube
    # Type L is most common for residential/commercial plumbing
    COPPER_TYPE_L_SPECS = {
        0.375: (12.7, 0.762, 11.176),      # 3/8"
        0.500: (15.875, 0.889, 14.097),    # 1/2"
        0.625: (19.05, 1.016, 17.018),     # 5/8"
        0.750: (22.225, 1.143, 19.939),    # 3/4"
        1.000: (28.575, 1.270, 26.035),    # 1"
        1.125: (28.575, 1.270, 26.035),    # 1-1/8" (same as 1" tube)
        1.250: (34.925, 1.397, 32.131),    # 1-1/4"
        1.500: (41.275, 1.588, 38.099),    # 1-1/2"
        2.000: (53.975, 1.778, 50.419),    # 2"
        2.500: (66.675, 2.032, 62.611),    # 2-1/2"
        3.000: (79.375, 2.286, 74.803),    # 3"
        3.125: (79.375, 2.286, 74.803),    # 3-1/8" (same as 3" tube)
        3.500: (92.075, 2.540, 86.995),    # 3-1/2"
        4.000: (104.775, 2.794, 99.187),   # 4"
        5.000: (130.175, 3.175, 123.825),  # 5"
        6.000: (155.575, 3.556, 148.463),  # 6"
    }
    
    # Fixture demands based on IPC/UPC maximum flow rates
    # Format: keyword -> (demand_m3s, fixture_units)
    # Reference: Lee et al. (2023), Premise Plumbing Modeling, Chapter 3
    FIXTURE_DEMANDS = {
        'lavatory': (0.000126, 1.0),      # 2.0 GPM
        'lav': (0.000126, 1.0),
        'sink': (0.000158, 1.5),          # 2.5 GPM
        'kitchen': (0.000158, 1.5),
        'shower': (0.000158, 2.0),        # 2.5 GPM
        'bath': (0.000252, 2.0),          # 4.0 GPM
        'tub': (0.000252, 2.0),
        'water closet': (0.0, 2.5),       # Intermittent (flush)
        'wc': (0.0, 2.5),
        'toilet': (0.0, 2.5),
        'urinal': (0.0, 2.0),
        'washer': (0.000252, 2.0),        # 4.0 GPM
        'laundry': (0.000252, 2.0),
        'hose bib': (0.000158, 2.5),      # 2.5 GPM
        'hose': (0.000158, 2.5),
        'bibb': (0.000158, 2.5),
        'drinking': (0.000063, 0.5),      # 1.0 GPM
        'water heater': (0.0, 0.0),
        'heater': (0.0, 0.0),
    }
    
    # Default fixture demand if type not recognized
    DEFAULT_FIXTURE_DEMAND = (0.000126, 1.0)  # 2.0 GPM, 1 FU


# =============================================================================
# PIPE DIAMETER CONVERTER
# =============================================================================
class PipeDiameterConverter:
    """
    Convert pipe outside diameter (OD) to inside diameter (ID).
    
    IFC/BIM models typically store geometry as OD (what you see in 3D).
    Hydraulic calculations require ID for accurate friction loss calculations.
    
    Uses Copper Type L wall thickness lookup table.
    Reference: ASTM B88
    """
    
    @staticmethod
    def od_to_id(od_m):
        """
        Convert outside diameter to inside diameter using Type L copper specs.
        
        Args:
            od_m: Outside diameter in meters
            
        Returns:
            float: Inside diameter in meters
        """
        if od_m is None or od_m <= 0:
            return Config.DEFAULT_DIAMETER
        
        od_mm = od_m * 1000  # Convert to mm for lookup
        
        # Find closest matching OD in lookup table
        best_match = None
        best_diff = float('inf')
        
        for nominal, (od_spec, wall, id_spec) in Config.COPPER_TYPE_L_SPECS.items():
            diff = abs(od_mm - od_spec)
            if diff < best_diff:
                best_diff = diff
                best_match = (od_spec, wall, id_spec)
        
        if best_match and best_diff < 5.0:  # Within 5mm tolerance
            od_spec, wall_thickness, id_mm = best_match
            return id_mm / 1000  # Convert back to meters
        else:
            # If no match found, estimate wall thickness as ~5% of OD
            # (conservative estimate for copper tubing)
            estimated_wall = od_mm * 0.05
            id_mm = od_mm - (2 * estimated_wall)
            return max(id_mm / 1000, 0.010)  # Minimum 10mm ID
    
    @staticmethod
    def get_wall_thickness(od_m):
        """
        Get wall thickness for a given OD.
        
        Args:
            od_m: Outside diameter in meters
            
        Returns:
            float: Wall thickness in meters
        """
        if od_m is None or od_m <= 0:
            return 0.001  # Default 1mm
        
        od_mm = od_m * 1000
        
        # Find closest matching OD
        best_match = None
        best_diff = float('inf')
        
        for nominal, (od_spec, wall, id_spec) in Config.COPPER_TYPE_L_SPECS.items():
            diff = abs(od_mm - od_spec)
            if diff < best_diff:
                best_diff = diff
                best_match = wall
        
        if best_match and best_diff < 5.0:
            return best_match / 1000  # Convert to meters
        else:
            return od_mm * 0.05 / 1000  # Estimate 5% of OD


# =============================================================================
# DIVERSITY FACTOR CALCULATOR
# =============================================================================
class DiversityFactor:
    """
    Calculate diversity factor based on Hunter's Curve concept.
    
    As the number of fixtures increases, the probability that all operate
    simultaneously decreases. This class provides methods to estimate
    realistic simultaneous demand.
    
    Reference: Hunter (1940), NBS Building Materials and Structures Report BMS65
    
    Note: Updated in v1.1.0 to use more realistic (lower) diversity factors
    based on industry feedback indicating actual simultaneous use is often
    15-20% for residential applications.
    """
    
    @staticmethod
    def calculate(num_fixtures, fixture_units_total=None):
        """
        Calculate diversity factor based on number of fixtures.
        
        Args:
            num_fixtures: Number of active fixtures in system
            fixture_units_total: Total fixture units (optional)
        
        Returns:
            float: Diversity factor between 0.0 and 1.0
        
        Note: v1.1.0 - Updated values based on industry feedback (Jeremy Williams)
              that actual diversity can be as low as 15-20% for residential.
        """
        if num_fixtures <= 1:
            return 1.0
        elif num_fixtures <= 5:
            return 0.7
        elif num_fixtures <= 10:
            return 0.5
        elif num_fixtures <= 20:
            return 0.30
        elif num_fixtures <= 50:
            return 0.20   # Updated from 0.30 (v1.0.0)
        elif num_fixtures <= 100:
            return 0.18   # Updated from 0.25 (v1.0.0)
        else:
            return 0.15   # Updated from 0.20 (v1.0.0)
    
    @staticmethod
    def hunter_curve_gpm(fixture_units):
        """
        Estimate peak demand from Hunter's curve (simplified).
        
        Args:
            fixture_units: Total fixture units
        
        Returns:
            float: Estimated peak demand in GPM
        """
        if fixture_units <= 0:
            return 0
        elif fixture_units <= 10:
            return fixture_units * 1.0
        elif fixture_units <= 50:
            return 10 + (fixture_units - 10) * 0.8
        elif fixture_units <= 200:
            return 42 + (fixture_units - 50) * 0.6
        else:
            return 132 + (fixture_units - 200) * 0.4


# =============================================================================
# IFC DATA EXTRACTION
# =============================================================================
class IFCExtractor:
    """
    Extract plumbing system data from IFC files.
    
    This class handles parsing of IFC files, extraction of plumbing elements
    (pipes, fittings, terminals), port-based connectivity mapping, and
    3D coordinate extraction from IFC geometry representations.
    
    Attributes:
        ifc_path: Path to IFC file
        ifc_file: Loaded IFC file object
        port_to_element: Mapping of port GlobalIds to element GlobalIds
        element_to_system: Mapping of element GlobalIds to system names
        element_coords: Dictionary of element coordinates from geometry
    """
    
    def __init__(self, ifc_path):
        """
        Initialize IFC extractor.
        
        Args:
            ifc_path: Path to IFC file
        """
        self.ifc_path = ifc_path
        self.ifc_file = ifcopenshell.open(ifc_path)
        self.port_to_element = {}
        self.element_to_ports = {}
        self.element_to_system = {}
        self.element_coords = {}
        self.all_elements = {}
        
        # Setup geometry settings for coordinate extraction
        self.geom_settings = ifcopenshell.geom.settings()
        self.geom_settings.set(self.geom_settings.USE_WORLD_COORDS, True)
    
    def get_file_info(self):
        """Get basic IFC file information."""
        return {
            'schema': self.ifc_file.schema,
            'total_entities': len(list(self.ifc_file))
        }
    
    def get_systems(self):
        """Extract system definitions from IFC."""
        systems = self.ifc_file.by_type("IfcSystem")
        return [{'name': s.Name, 'description': s.Description} for s in systems]
    
    def build_system_membership(self):
        """
        Build mapping of elements to their assigned MEP systems.
        
        Uses IfcRelAssignsToGroup relationships to determine which
        system (DCW, DHW, etc.) each element belongs to.
        """
        systems = self.ifc_file.by_type("IfcSystem")
        system_info = {}
        
        for sys in systems:
            members = []
            for rel in self.ifc_file.by_type("IfcRelAssignsToGroup"):
                if rel.RelatingGroup == sys:
                    for obj in rel.RelatedObjects:
                        members.append(obj.GlobalId)
                        self.element_to_system[obj.GlobalId] = sys.Name
            system_info[sys.Name] = {'members': members, 'count': len(members)}
        
        return system_info
    
    def is_cold_water_system(self, system_name):
        """Check if system name indicates cold water (DCW)."""
        if system_name is None:
            return False
        system_upper = system_name.upper()
        return any(kw in system_upper for kw in Config.TARGET_SYSTEMS)
    
    def _get_element_center(self, element):
        """
        Extract center coordinates from IFC geometry representation.
        
        Uses ifcopenshell.geom to process the element's geometry and
        calculate the centroid of all vertices.
        
        Args:
            element: IFC element object
            
        Returns:
            numpy.ndarray: [x, y, z] coordinates or None if extraction fails
        """
        try:
            shape = ifcopenshell.geom.create_shape(self.geom_settings, element)
            verts = np.array(shape.geometry.verts).reshape(-1, 3)
            return verts.mean(axis=0)
        except Exception:
            return None
    
    def _get_element_properties(self, element):
        """Extract property set values from IFC element."""
        props = {}
        for rel in element.IsDefinedBy:
            if rel.is_a("IfcRelDefinesByProperties"):
                pset = rel.RelatingPropertyDefinition
                if pset.is_a("IfcPropertySet"):
                    for prop in pset.HasProperties:
                        if prop.is_a("IfcPropertySingleValue") and prop.NominalValue:
                            props[prop.Name] = prop.NominalValue.wrappedValue
        return props
    
    def _get_pipe_geometry(self, element):
        """
        Extract diameter and length from pipe geometry.
        
        Parses IfcExtrudedAreaSolid representation to get:
        - Diameter from IfcCircleProfileDef radius (OUTSIDE diameter)
        - Length from extrusion depth
        
        Note: v1.1.0 - Now converts OD to ID using copper Type L wall thickness
        
        Returns:
            tuple: (inside_diameter_m, length_m, outside_diameter_m) or (None, None, None)
        """
        od_m = None
        length_m = None
        
        if element.Representation:
            for rep in element.Representation.Representations:
                for item in rep.Items:
                    if item.is_a("IfcExtrudedAreaSolid"):
                        # Length from extrusion depth (assumes feet -> meters)
                        length_m = item.Depth * 0.3048
                        profile = item.SweptArea
                        if hasattr(profile, 'Radius'):
                            # This is OUTSIDE radius from IFC geometry
                            od_m = profile.Radius * 2 * 0.3048
        
        # Convert OD to ID using copper Type L wall thickness
        if od_m is not None:
            id_m = PipeDiameterConverter.od_to_id(od_m)
            return id_m, length_m, od_m
        
        return None, length_m, None
    
    def _identify_fixture_type(self, name, obj_type):
        """Identify fixture type and demand from terminal name/type."""
        search_text = f"{name or ''} {obj_type or ''}".lower()
        
        for keyword, (demand, fu) in Config.FIXTURE_DEMANDS.items():
            if keyword in search_text:
                return keyword, demand, fu
        
        return 'unknown', Config.DEFAULT_FIXTURE_DEMAND[0], Config.DEFAULT_FIXTURE_DEMAND[1]
    
    def build_port_mapping(self):
        """
        Build mapping between distribution ports and elements.
        
        Uses IfcRelConnectsPortToElement relationships to create
        bidirectional mapping for topology reconstruction.
        """
        for rel in self.ifc_file.by_type("IfcRelConnectsPortToElement"):
            port = rel.RelatingPort
            element = rel.RelatedElement
            self.port_to_element[port.GlobalId] = element.GlobalId
            
            if element.GlobalId not in self.element_to_ports:
                self.element_to_ports[element.GlobalId] = []
            self.element_to_ports[element.GlobalId].append(port.GlobalId)
        
        return len(self.port_to_element)
    
    def get_connections(self):
        """
        Extract element connectivity from port connections.
        
        Uses IfcRelConnectsPorts to determine which elements are
        physically connected through their distribution ports.
        
        Returns:
            list: List of (element1_id, element2_id) tuples
        """
        connections = []
        for rel in self.ifc_file.by_type("IfcRelConnectsPorts"):
            port1 = rel.RelatingPort
            port2 = rel.RelatedPort
            elem1_id = self.port_to_element.get(port1.GlobalId)
            elem2_id = self.port_to_element.get(port2.GlobalId)
            if elem1_id and elem2_id:
                connections.append((elem1_id, elem2_id))
        return connections
    
    def extract_pipes(self):
        """Extract all pipe segments (IfcFlowSegment)."""
        pipes = []
        for seg in self.ifc_file.by_type("IfcFlowSegment"):
            props = self._get_element_properties(seg)
            id_m, length_m, od_m = self._get_pipe_geometry(seg)
            system_name = self.element_to_system.get(seg.GlobalId, "UNASSIGNED")
            
            # Get real coordinates
            coords = self._get_element_center(seg)
            if coords is not None:
                self.element_coords[seg.GlobalId] = coords
            
            pipe_data = {
                'id': seg.GlobalId,
                'name': seg.Name,
                'diameter_m': id_m,        # Inside diameter (for hydraulics)
                'diameter_od_m': od_m,     # Outside diameter (from IFC)
                'length_m': length_m,
                'elevation_m': props.get('InvertElevation', 0) * 0.3048,
                'system': system_name
            }
            pipes.append(pipe_data)
            self.all_elements[seg.GlobalId] = {'type': 'pipe', 'data': pipe_data}
        
        return pipes
    
    def extract_fittings(self):
        """Extract all fittings (IfcFlowFitting)."""
        fittings = []
        for fit in self.ifc_file.by_type("IfcFlowFitting"):
            system_name = self.element_to_system.get(fit.GlobalId, "UNASSIGNED")
            
            coords = self._get_element_center(fit)
            if coords is not None:
                self.element_coords[fit.GlobalId] = coords
            
            fitting_data = {
                'id': fit.GlobalId,
                'name': fit.Name,
                'type': fit.ObjectType,
                'system': system_name
            }
            fittings.append(fitting_data)
            self.all_elements[fit.GlobalId] = {'type': 'fitting', 'data': fitting_data}
        
        return fittings
    
    def extract_terminals(self):
        """Extract all terminals/fixtures (IfcFlowTerminal)."""
        terminals = []
        for term in self.ifc_file.by_type("IfcFlowTerminal"):
            system_name = self.element_to_system.get(term.GlobalId, "UNASSIGNED")
            fixture_type, demand, fu = self._identify_fixture_type(term.Name, term.ObjectType)
            
            coords = self._get_element_center(term)
            if coords is not None:
                self.element_coords[term.GlobalId] = coords
            
            terminal_data = {
                'id': term.GlobalId,
                'name': term.Name,
                'type': term.ObjectType,
                'system': system_name,
                'fixture_type': fixture_type,
                'base_demand_m3s': demand,
                'fixture_units': fu
            }
            terminals.append(terminal_data)
            self.all_elements[term.GlobalId] = {'type': 'terminal', 'data': terminal_data}
        
        return terminals
    
    def extract_all(self):
        """
        Run complete IFC extraction pipeline.
        
        Returns:
            dict: Extracted data including pipes, fittings, terminals, connections
        """
        print("=" * 60)
        print("IFC DATA EXTRACTION")
        print("=" * 60)
        
        file_info = self.get_file_info()
        print(f"Schema: {file_info['schema']}")
        print(f"Entities: {file_info['total_entities']:,}")
        
        # Build system membership first
        system_info = self.build_system_membership()
        print(f"\nSystems found: {len(system_info)}")
        for name, info in system_info.items():
            dcw = "[DCW]" if self.is_cold_water_system(name) else ""
            print(f"  {name}: {info['count']} elements {dcw}")
        
        # Build port mapping
        num_ports = self.build_port_mapping()
        print(f"\nPort mapping: {num_ports} ports")
        
        # Extract elements
        pipes = self.extract_pipes()
        fittings = self.extract_fittings()
        terminals = self.extract_terminals()
        connections = self.get_connections()
        
        print(f"\nElements extracted:")
        print(f"  Pipes: {len(pipes)}")
        print(f"  Fittings: {len(fittings)}")
        print(f"  Terminals: {len(terminals)}")
        print(f"  Connections: {len(connections)}")
        
        # Count elements with coordinates
        coords_count = len(self.element_coords)
        total_elems = len(pipes) + len(fittings) + len(terminals)
        print(f"  With coordinates: {coords_count}/{total_elems} ({100*coords_count/max(1,total_elems):.0f}%)")
        
        # Report diameter statistics (now showing ID)
        diameters = [p['diameter_m'] for p in pipes if p['diameter_m']]
        if diameters:
            print(f"\nPipe inside diameters (ID):")
            print(f"  Range: {min(diameters)*1000:.1f} - {max(diameters)*1000:.1f} mm")
            print(f"  Note: Converted from OD using Type L copper wall thickness")
        
        return {
            'file_info': file_info,
            'systems': system_info,
            'pipes': pipes,
            'fittings': fittings,
            'terminals': terminals,
            'connections': connections
        }


# =============================================================================
# TOPOLOGY ANALYZER
# =============================================================================
class TopologyAnalyzer:
    """
    Analyze network topology from extracted IFC data.
    
    Uses NetworkX to identify connected components and filter
    for DCW (Domestic Cold Water) systems.
    """
    
    def __init__(self, extracted_data, extractor):
        self.extracted_data = extracted_data
        self.extractor = extractor
        self.all_elements = extractor.all_elements
        self.graph = self._build_graph()
    
    def _build_graph(self):
        """Build NetworkX graph from connections."""
        G = nx.Graph()
        
        for elem_id, elem in self.all_elements.items():
            G.add_node(elem_id, **elem)
        
        for e1, e2 in self.extracted_data['connections']:
            if e1 in G and e2 in G:
                G.add_edge(e1, e2)
        
        return G
    
    def get_components(self, dcw_only=True, dcw_threshold=90.0):
        """
        Get connected components, optionally filtered for DCW systems.
        
        Args:
            dcw_only: If True, only return components that are ≥dcw_threshold% DCW
            dcw_threshold: Minimum percentage of DCW elements to qualify
        
        Returns:
            list: Component info dictionaries
        """
        print("\n" + "=" * 60)
        print("TOPOLOGY ANALYSIS")
        print("=" * 60)
        
        components = list(nx.connected_components(self.graph))
        print(f"Connected components: {len(components)}")
        
        component_info = []
        dcw_components = []
        
        for i, comp in enumerate(sorted(components, key=len, reverse=True)):
            # Count element types and systems
            pipe_count = 0
            fitting_count = 0
            terminal_count = 0
            dcw_count = 0
            total_fu = 0.0
            system_counts = defaultdict(int)
            
            for elem_id in comp:
                elem = self.all_elements.get(elem_id)
                if not elem:
                    continue
                    
                elem_type = elem['type']
                elem_data = elem['data']
                
                if elem_type == 'pipe':
                    pipe_count += 1
                elif elem_type == 'fitting':
                    fitting_count += 1
                elif elem_type == 'terminal':
                    terminal_count += 1
                    total_fu += elem_data.get('fixture_units', 0)
                
                system = elem_data.get('system', 'UNKNOWN')
                system_counts[system] += 1
                if self.extractor.is_cold_water_system(system):
                    dcw_count += 1
            
            dcw_pct = (dcw_count / len(comp) * 100) if comp else 0
            
            info = {
                'index': i,
                'elements': comp,
                'pipes': pipe_count,
                'fittings': fitting_count,
                'terminals': terminal_count,
                'dcw_percentage': dcw_pct,
                'system_breakdown': dict(system_counts),
                'fixture_units': total_fu
            }
            component_info.append(info)
            
            dcw_marker = "[DCW]" if dcw_pct >= dcw_threshold else "[Mixed]"
            print(f"\n  Component {i+1}: {len(comp)} elements {dcw_marker}")
            print(f"    Pipes: {pipe_count}, Fittings: {fitting_count}, Terminals: {terminal_count}")
            print(f"    Fixture units: {total_fu:.1f}, DCW: {dcw_pct:.0f}%")
            
            if dcw_pct >= dcw_threshold:
                dcw_components.append(info)
        
        if dcw_only:
            print(f"\n  DCW components selected: {len(dcw_components)}")
            return dcw_components
        
        return component_info


# =============================================================================
# WNTR MODEL BUILDER
# =============================================================================
class WNTRBuilder:
    """
    Build WNTR hydraulic models from IFC extracted data.
    
    Creates WNTR network with:
    - Junctions at each element location
    - Pipes based on port connectivity
    - Real 3D coordinates from IFC geometry
    - Reservoir at system inlet
    """
    
    def __init__(self, all_elements, connections, element_coords):
        self.all_elements = all_elements
        self.connections = connections
        self.element_coords = element_coords
        self.junction_to_element = {}
    
    def build_from_component(self, component):
        """
        Build WNTR model from a connected component.
        
        Args:
            component: Set of element GlobalIds in the component
            
        Returns:
            tuple: (WaterNetworkModel, junction_map)
        """
        wn = wntr.network.WaterNetworkModel()
        
        comp_set = set(component)
        comp_connections = [(e1, e2) for e1, e2 in self.connections 
                           if e1 in comp_set and e2 in comp_set]
        
        junction_map = {}
        self.junction_to_element = {}
        node_counter = 0
        
        # Create junctions with real coordinates
        for elem_id in component:
            elem = self.all_elements[elem_id]
            jname = f"J{node_counter}"
            junction_map[elem_id] = jname
            self.junction_to_element[jname] = elem_id
            
            # Get elevation from IFC coordinates or property
            elevation = 0
            if elem_id in self.element_coords:
                elevation = self.element_coords[elem_id][2]  # Z coordinate
            elif elem['type'] == 'pipe':
                elevation = elem['data'].get('elevation_m', 0)
            
            wn.add_junction(jname, base_demand=0, elevation=elevation)
            
            # Assign real coordinates for EPANET visualization
            if elem_id in self.element_coords:
                coords = self.element_coords[elem_id]
                # Use actual IFC coordinates (meters)
                # X: horizontal position (IFC X + Y for plan projection)
                # Y: elevation (IFC Z)
                x_epanet = coords[0] + coords[1]  # Combined X+Y in meters
                y_epanet = coords[2]               # Z elevation in meters
                wn.get_node(jname).coordinates = (x_epanet, y_epanet)
            
            node_counter += 1
        
        # Create pipes
        pipe_counter = 0
        added_pipes = set()
        
        for e1, e2 in comp_connections:
            j1, j2 = junction_map.get(e1), junction_map.get(e2)
            if not j1 or not j2 or j1 == j2:
                continue
            
            pipe_key = tuple(sorted([j1, j2]))
            if pipe_key in added_pipes:
                continue
            
            # Get pipe properties
            length_m = Config.MIN_PIPE_LENGTH
            diameter_m = Config.DEFAULT_DIAMETER
            
            # Calculate length from coordinates if available
            if e1 in self.element_coords and e2 in self.element_coords:
                c1, c2 = self.element_coords[e1], self.element_coords[e2]
                length_m = max(np.linalg.norm(c1 - c2), Config.MIN_PIPE_LENGTH)
            
            # Get diameter from pipe data (now using ID from v1.1.0)
            for elem_id in [e1, e2]:
                if self.all_elements[elem_id]['type'] == 'pipe':
                    pipe_data = self.all_elements[elem_id]['data']
                    if pipe_data.get('diameter_m'):
                        diameter_m = pipe_data['diameter_m']  # This is now ID
                    if pipe_data.get('length_m'):
                        length_m = max(pipe_data['length_m'], Config.MIN_PIPE_LENGTH)
                    break
            
            pname = f"P{pipe_counter}"
            wn.add_pipe(pname, j1, j2, 
                       length=length_m, 
                       diameter=diameter_m,
                       roughness=Config.ROUGHNESS_COPPER)
            added_pipes.add(pipe_key)
            pipe_counter += 1
        
        # Add reservoir at source (lowest elevation leaf node)
        degree_1_nodes = [n for n in wn.junction_name_list 
                         if len(wn.get_links_for_node(n)) == 1]
        
        if degree_1_nodes:
            source_junction = min(degree_1_nodes,
                                 key=lambda n: wn.get_node(n).elevation)
            wn.add_reservoir('R1', base_head=Config.DEFAULT_SUPPLY_HEAD)
            
            # Position reservoir below the network (in meters)
            source_coords = wn.get_node(source_junction).coordinates
            wn.get_node('R1').coordinates = (source_coords[0], source_coords[1] - 5)
            
            wn.add_pipe('P_source', 'R1', source_junction,
                       length=1, diameter=0.05, roughness=Config.ROUGHNESS_COPPER)
        
        return wn, junction_map


# =============================================================================
# HYDRAULIC SIMULATOR
# =============================================================================
class HydraulicSimulator:
    """
    Run hydraulic simulations with fixture-based demands.
    
    Applies fixture demands based on terminal types and diversity
    factors based on Hunter's curve methodology.
    """
    
    def __init__(self, wn, all_elements, junction_to_element):
        self.wn = wn
        self.all_elements = all_elements
        self.junction_to_element = junction_to_element
        self.results = None
        self.demand_summary = {}
    
    def add_fixture_demands(self, apply_diversity=True):
        """
        Add demands to fixture nodes based on fixture type.
        
        Args:
            apply_diversity: Whether to apply diversity factor
            
        Returns:
            dict: Summary of demand allocation
        """
        terminal_junctions = []
        total_base_demand = 0
        total_fixture_units = 0
        fixture_counts = defaultdict(int)
        
        for jname in self.wn.junction_name_list:
            elem_id = self.junction_to_element.get(jname)
            if not elem_id:
                continue
                
            elem = self.all_elements.get(elem_id)
            if elem and elem['type'] == 'terminal':
                terminal_data = elem['data']
                base_demand = terminal_data.get('base_demand_m3s', 0)
                fixture_units = terminal_data.get('fixture_units', 0)
                fixture_type = terminal_data.get('fixture_type', 'unknown')
                
                if base_demand > 0:
                    terminal_junctions.append({
                        'junction': jname,
                        'base_demand': base_demand,
                        'fixture_units': fixture_units,
                        'fixture_type': fixture_type
                    })
                    total_base_demand += base_demand
                    total_fixture_units += fixture_units
                    fixture_counts[fixture_type] += 1
        
        num_active = len(terminal_junctions)
        
        # Calculate diversity factor (updated in v1.1.0)
        if apply_diversity and num_active > 1:
            diversity_factor = DiversityFactor.calculate(num_active, total_fixture_units)
        else:
            diversity_factor = 1.0
        
        # Apply demands
        total_applied = 0
        self.demand_junctions = {}  # Store for INP file post-processing
        for term in terminal_junctions:
            junction = self.wn.get_node(term['junction'])
            applied_demand = term['base_demand'] * diversity_factor
            junction.demand_timeseries_list[0].base_value = applied_demand
            self.demand_junctions[term['junction']] = applied_demand
            total_applied += applied_demand
        
        self.demand_summary = {
            'num_fixtures': num_active,
            'fixture_units': total_fixture_units,
            'diversity_factor': diversity_factor,
            'total_base_demand_Ls': total_base_demand * 1000,
            'total_applied_demand_Ls': total_applied * 1000,
            'fixture_breakdown': dict(fixture_counts)
        }
        
        return self.demand_summary
    
    def run_steady_state(self):
        """Run steady-state hydraulic simulation."""
        self.wn.options.time.duration = 0
        self.wn.options.hydraulic.demand_model = 'DD'
        
        try:
            sim = wntr.sim.EpanetSimulator(self.wn)
            self.results = sim.run_sim()
            return True
        except Exception as e:
            print(f"EPANET simulation failed: {e}")
            try:
                sim = wntr.sim.WNTRSimulator(self.wn)
                self.results = sim.run_sim()
                return True
            except Exception as e2:
                print(f"WNTR simulation failed: {e2}")
                return False
    
    def get_results_summary(self):
        """Get summary of simulation results."""
        if self.results is None:
            return None
        
        pressure = self.results.node['pressure'].loc[0]
        flow = self.results.link['flowrate'].loc[0]
        velocity = self.results.link['velocity'].loc[0]
        
        junc_pressure = pressure[self.wn.junction_name_list]
        pipe_flow = flow[self.wn.pipe_name_list]
        pipe_velocity = velocity[self.wn.pipe_name_list]
        
        return {
            'pressure': {
                'min_m': float(junc_pressure.min()),
                'max_m': float(junc_pressure.max()),
                'mean_m': float(junc_pressure.mean()),
                'min_psi': float(junc_pressure.min() * 1.42),
                'max_psi': float(junc_pressure.max() * 1.42),
                'mean_psi': float(junc_pressure.mean() * 1.42),
            },
            'flow': {
                'max_Ls': float(pipe_flow.abs().max() * 1000),
                'mean_Ls': float(pipe_flow.abs().mean() * 1000),
            },
            'velocity': {
                'max_ms': float(pipe_velocity.abs().max()),
                'mean_ms': float(pipe_velocity.abs().mean()),
            },
            'demand': self.demand_summary
        }
    
    def update_inp_file_demands(self, inp_path):
        """
        Update INP file with demands and fix units for SI compatibility.
        
        WNTR writes UNITS=GPM with lengths in feet and diameters in inches.
        This method:
        1. Changes UNITS to LPS (SI-compatible)
        2. Converts lengths from feet to meters
        3. Converts diameters from inches to millimeters
        4. Writes demands in L/s
        """
        if not hasattr(self, 'demand_junctions'):
            return
        
        with open(inp_path, 'r') as f:
            content = f.read()
        
        # Update units from GPM to LPS
        content = content.replace('Units              	GPM', 'Units              	LPS')
        
        # Update JUNCTIONS section with demands
        lines = content.split('\n')
        new_lines = []
        in_junctions = False
        in_pipes = False
        
        for line in lines:
            if line.strip().startswith('[JUNCTIONS]'):
                in_junctions = True
                new_lines.append(line)
                continue
            elif line.strip().startswith('[') and in_junctions:
                in_junctions = False
            
            if line.strip().startswith('[PIPES]'):
                in_pipes = True
                new_lines.append(line)
                continue
            elif line.strip().startswith('[') and in_pipes:
                in_pipes = False
            
            if in_junctions and not line.strip().startswith(';') and line.strip():
                parts = line.split()
                if len(parts) >= 2 and parts[0] in self.demand_junctions:
                    # Update demand in L/s
                    demand_Ls = self.demand_junctions[parts[0]] * 1000
                    new_line = f"{parts[0]:<16}{parts[1]:<12}{demand_Ls:.6f}"
                    new_lines.append(new_line)
                    continue
            
            if in_pipes and not line.strip().startswith(';') and line.strip():
                parts = line.split()
                if len(parts) >= 6:
                    # Convert length from feet to meters, diameter from inches to mm
                    try:
                        length_ft = float(parts[3])
                        diameter_in = float(parts[4])
                        length_m = length_ft * 0.3048
                        diameter_mm = diameter_in * 25.4
                        new_line = f"{parts[0]:<16}{parts[1]:<16}{parts[2]:<16}{length_m:<12.3f}{diameter_mm:<12.1f}{parts[5]}"
                        if len(parts) > 6:
                            new_line += "  " + "  ".join(parts[6:])
                        new_lines.append(new_line)
                        continue
                    except (ValueError, IndexError):
                        pass
            
            new_lines.append(line)
        
        with open(inp_path, 'w') as f:
            f.write('\n'.join(new_lines))


# =============================================================================
# VISUALIZER
# =============================================================================
class Visualizer:
    """Publication-quality visualization for BIM2WNTR results."""
    
    def __init__(self, output_dir):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set publication style
        plt.rcParams.update({
            'font.size': 10,
            'font.family': 'sans-serif',
            'axes.labelsize': 10,
            'axes.titlesize': 11,
            'xtick.labelsize': 9,
            'ytick.labelsize': 9,
            'legend.fontsize': 9,
            'figure.dpi': 150,
            'savefig.dpi': 300,
            'savefig.bbox': 'tight',
        })
        
        # Color palette - using seaborn's muted palette
        self.colors = sns.color_palette("muted", 8)
        self.cmap = 'viridis'
    
    def plot_network_layout(self, wn, system_name='DCW', filename='figure_network_layout.png'):
        """Plot network layout colored by elevation."""
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Get node data
        junctions = wn.junction_name_list
        elevations = [wn.get_node(j).elevation for j in junctions]
        coords = [wn.get_node(j).coordinates for j in junctions]
        x = [c[0] for c in coords]
        y = [c[1] for c in coords]
        
        # Plot pipes
        for pname in wn.pipe_name_list:
            pipe = wn.get_link(pname)
            n1, n2 = pipe.start_node_name, pipe.end_node_name
            c1 = wn.get_node(n1).coordinates
            c2 = wn.get_node(n2).coordinates
            ax.plot([c1[0], c2[0]], [c1[1], c2[1]], 'gray', linewidth=0.8, alpha=0.6)
        
        # Plot junctions colored by elevation
        scatter = ax.scatter(x, y, c=elevations, cmap=self.cmap, s=30, zorder=5)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.15)
        cbar = plt.colorbar(scatter, cax=cax, label='Elevation (m)')
        
        # Plot reservoir
        if 'R1' in wn.node_name_list:
            r_coords = wn.get_node('R1').coordinates
            ax.scatter([r_coords[0]], [r_coords[1]], c='blue', s=100, 
                      marker='s', zorder=10, label='Reservoir')
        
        ax.set_xlabel('Horizontal Position (m)')
        ax.set_ylabel('Elevation (m)')
        ax.set_title(f'Network Layout - {system_name} Distribution System')
        ax.set_aspect('equal')
        ax.legend()
        
        plt.tight_layout()
        filepath = self.output_dir / filename
        plt.savefig(filepath, dpi=300, bbox_inches='tight', facecolor='white')
        plt.savefig(filepath.with_suffix('.pdf'), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {filepath}")
        return filepath
    
    def plot_network_analysis(self, wn, system_name='DCW', filename='figure_network_analysis.png'):
        """Plot network characteristics (4-panel figure)."""
        fig = plt.figure(figsize=(12, 10))
        
        # (a) Side elevation view
        ax1 = fig.add_subplot(2, 2, 1)
        junctions = wn.junction_name_list
        coords = [wn.get_node(j).coordinates for j in junctions]
        x = [c[0] for c in coords]
        y = [c[1] for c in coords]  # This is elevation
        
        ax1.scatter(x, y, c=self.colors[0], s=15, alpha=0.6)
        ax1.set_xlabel('Horizontal Position (m)')
        ax1.set_ylabel('Elevation (m)')
        ax1.set_title('(a) Side Elevation View')
        
        # (b) Pipe diameter distribution (now showing ID)
        ax2 = fig.add_subplot(2, 2, 2)
        diameters_mm = [wn.get_link(p).diameter * 1000 for p in wn.pipe_name_list]
        ax2.hist(diameters_mm, bins=15, color=self.colors[1], edgecolor='black', alpha=0.8)
        ax2.set_xlabel('Inside Diameter (mm)')
        ax2.set_ylabel('Number of Pipes')
        ax2.set_title('(b) Pipe Inside Diameter Distribution')
        
        # (c) Pipe length distribution
        ax3 = fig.add_subplot(2, 2, 3)
        lengths_m = [wn.get_link(p).length for p in wn.pipe_name_list]
        ax3.hist(lengths_m, bins=25, color=self.colors[2], edgecolor='black', alpha=0.8)
        ax3.set_xlabel('Length (m)')
        ax3.set_ylabel('Number of Pipes')
        mean_length = np.mean(lengths_m)
        ax3.axvline(mean_length, color='red', linestyle='--', label=f'Mean = {mean_length:.2f} m')
        ax3.legend()
        ax3.set_title('(c) Pipe Length Distribution')
        
        # (d) Summary statistics
        ax4 = fig.add_subplot(2, 2, 4)
        ax4.axis('off')
        
        total_length = sum(lengths_m)
        stats_data = [
            ['Parameter', 'Value'],
            ['Junctions', f"{len(junctions)}"],
            ['Pipes', f"{len(wn.pipe_name_list)}"],
            ['Total Length', f"{total_length:.1f} m"],
            ['ID Range', f"{min(diameters_mm):.1f} - {max(diameters_mm):.1f} mm"],
            ['Mean Pipe Length', f"{mean_length:.2f} m"],
        ]
        
        table = ax4.table(cellText=stats_data[1:], colLabels=stats_data[0],
                         loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.2, 1.5)
        ax4.set_title('(d) Network Statistics', pad=20)
        
        plt.suptitle(f'Network Characteristics - {system_name}', fontsize=12, fontweight='bold')
        plt.tight_layout()
        filepath = self.output_dir / filename
        plt.savefig(filepath, dpi=300, bbox_inches='tight', facecolor='white')
        plt.savefig(filepath.with_suffix('.pdf'), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {filepath}")
        return filepath
    
    def plot_hydraulic_results(self, wn, results, demand_summary, filename='figure_hydraulic_results.png'):
        """Plot hydraulic simulation results (4-panel figure)."""
        fig = plt.figure(figsize=(12, 10))
        
        # Get simulation data
        pressure = results.node['pressure'].loc[0]
        flow = results.link['flowrate'].loc[0]
        velocity = results.link['velocity'].loc[0]
        
        junc_pressure = pressure[wn.junction_name_list]
        pipe_flow = flow[wn.pipe_name_list]
        pipe_velocity = velocity[wn.pipe_name_list]
        
        # Convert units
        pipe_flow_Ls = pipe_flow.abs() * 1000  # m³/s to L/s
        
        # (a) Pressure distribution
        ax1 = fig.add_subplot(2, 2, 1)
        ax1.hist(junc_pressure.values, bins=25, color=self.colors[0], 
                edgecolor='black', alpha=0.8)
        ax1.axvline(14, color='red', linestyle='--', linewidth=2, 
                   label='Min recommended (14 m ≈ 20 psi)')
        ax1.set_xlabel('Pressure (m)')
        ax1.set_ylabel('Number of Junctions')
        ax1.set_title('(a) Junction Pressure Distribution')
        ax1.legend()
        
        # (b) Flow distribution
        ax2 = fig.add_subplot(2, 2, 2)
        ax2.hist(pipe_flow_Ls, bins=25, color=self.colors[1], 
                edgecolor='black', alpha=0.8)
        ax2.set_xlabel('Flow Rate (L/s)')
        ax2.set_ylabel('Number of Pipes')
        ax2.set_title('(b) Pipe Flow Rate Distribution')
        
        # (c) Velocity distribution
        ax3 = fig.add_subplot(2, 2, 3)
        ax3.hist(pipe_velocity, bins=25, color=self.colors[2], 
                edgecolor='black', alpha=0.8)
        ax3.axvline(2.5, color='red', linestyle='--', linewidth=2, 
                   label='Max recommended (2.5 m/s)')
        ax3.set_xlabel('Velocity (m/s)')
        ax3.set_ylabel('Number of Pipes')
        ax3.set_title('(c) Pipe Velocity Distribution')
        ax3.legend()
        
        # (d) Demand summary
        ax4 = fig.add_subplot(2, 2, 4)
        ax4.axis('off')
        
        demand_data = [
            ['Parameter', 'Value'],
            ['Active Fixtures', f"{demand_summary['num_fixtures']}"],
            ['Fixture Units', f"{demand_summary['fixture_units']:.1f}"],
            ['Diversity Factor', f"{demand_summary['diversity_factor']:.2f}"],
            ['Base Demand', f"{demand_summary['total_base_demand_Ls']:.3f} L/s"],
            ['Applied Demand', f"{demand_summary['total_applied_demand_Ls']:.3f} L/s"],
        ]
        
        table = ax4.table(cellText=demand_data[1:], colLabels=demand_data[0],
                         loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.2, 1.5)
        ax4.set_title('(d) Demand Summary', pad=20)
        
        plt.tight_layout()
        filepath = self.output_dir / filename
        plt.savefig(filepath, dpi=300, bbox_inches='tight', facecolor='white')
        plt.savefig(filepath.with_suffix('.pdf'), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {filepath}")
        return filepath


# =============================================================================
# MAIN PIPELINE
# =============================================================================
def run_pipeline(ifc_path, output_dir='./output'):
    """
    Run complete IFC to WNTR conversion pipeline.
    
    Args:
        ifc_path: Path to input IFC file
        output_dir: Output directory for results
        
    Returns:
        list: Generated WNTR models with metadata
    """
    print("\n" + "=" * 60)
    print(f"BIM2WNTR v{__version__}: IFC to Hydraulic Model Conversion")
    print("DCW Systems with Fixture-Based Demands")
    print("=" * 60)
    print(f"Input: {ifc_path}")
    print(f"Output: {output_dir}")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Extract IFC data
    extractor = IFCExtractor(ifc_path)
    extracted_data = extractor.extract_all()
    
    # Save extracted data
    json_data = {
        'pipes': extracted_data['pipes'],
        'fittings': extracted_data['fittings'],
        'terminals': extracted_data['terminals'],
        'connections': extracted_data['connections'],
    }
    with open(output_path / 'ifc_extracted.json', 'w') as f:
        json.dump(json_data, f, indent=2)
    print(f"\nSaved: {output_path / 'ifc_extracted.json'}")
    
    # Step 2: Analyze topology (DCW only)
    analyzer = TopologyAnalyzer(extracted_data, extractor)
    components = analyzer.get_components(dcw_only=True, dcw_threshold=90.0)
    
    if not components:
        print("\nWARNING: No DCW-only components found!")
        components = analyzer.get_components(dcw_only=False)
    
    # Step 3: Build WNTR models
    builder = WNTRBuilder(analyzer.all_elements, extracted_data['connections'], 
                         extractor.element_coords)
    
    models = []
    for comp_info in components:
        if comp_info['pipes'] < 5:
            continue
        
        wn, junction_map = builder.build_from_component(comp_info['elements'])
        
        primary_system = max(comp_info['system_breakdown'].items(), 
                            key=lambda x: x[1])[0]
        
        # Save INP file
        inp_filename = f"DCW_{comp_info['index']+1}_{primary_system.replace(' ', '_')}.inp"
        inp_path = output_path / inp_filename
        wntr.network.io.write_inpfile(wn, str(inp_path))
        print(f"Saved: {inp_path}")
        
        models.append({
            'wn': wn,
            'system': primary_system,
            'info': comp_info,
            'inp_path': inp_path,
            'junction_to_element': builder.junction_to_element.copy()
        })
    
    # Step 4: Run simulation on all DCW systems
    if models:
        print("\n" + "=" * 60)
        print("HYDRAULIC SIMULATION")
        print("=" * 60)
        
        for model in models:
            print(f"\nSystem: {model['system']}")
            print(f"  Pipes: {model['info']['pipes']}")
            print(f"  Fixture units: {model['info']['fixture_units']:.1f}")
            
            simulator = HydraulicSimulator(
                model['wn'], analyzer.all_elements, model['junction_to_element']
            )
            
            demand_summary = simulator.add_fixture_demands(apply_diversity=True)
            print(f"  Fixtures with demand: {demand_summary['num_fixtures']}")
            print(f"  Diversity factor: {demand_summary['diversity_factor']:.2f}")
            print(f"  Applied demand: {demand_summary['total_applied_demand_Ls']:.3f} L/s")
            
            # Re-save INP file WITH demands
            wntr.network.io.write_inpfile(model['wn'], str(model['inp_path']))
            simulator.update_inp_file_demands(str(model['inp_path']))
            print(f"  Updated: {model['inp_path']}")
            
            if simulator.run_steady_state():
                results_summary = simulator.get_results_summary()
                print(f"  Pressure: {results_summary['pressure']['min_m']:.1f} - "
                      f"{results_summary['pressure']['max_m']:.1f} m")
                print(f"  Max velocity: {results_summary['velocity']['max_ms']:.3f} m/s")
                
                model['results'] = results_summary
                model['simulator'] = simulator
        
        # Save results from largest model
        main_model = max(models, key=lambda m: m['info']['pipes'])
        if 'results' in main_model:
            with open(output_path / 'simulation_results.json', 'w') as f:
                json.dump(main_model['results'], f, indent=2)
            
            # Generate figures
            print("\n" + "=" * 60)
            print("GENERATING FIGURES")
            print("=" * 60)
            
            viz = Visualizer(output_path)
            viz.plot_network_layout(main_model['wn'], system_name=main_model['system'])
            viz.plot_network_analysis(main_model['wn'], system_name=main_model['system'])
            if 'simulator' in main_model and main_model['simulator'].results is not None:
                viz.plot_hydraulic_results(
                    main_model['wn'], 
                    main_model['simulator'].results, 
                    main_model['results']['demand']
                )
    
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    print(f"Output: {output_path}")
    
    return models


# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"""
BIM2WNTR v{__version__} - IFC to WNTR Hydraulic Model Converter
{__author__}

Usage:
    python bim2wntr.py <ifc_file> [output_dir]

Example:
    python bim2wntr.py building.ifc ./output

Features:
    - Domestic Cold Water (DCW) system extraction
    - Real 3D coordinates from IFC geometry
    - Fixture-based demands (IPC/UPC)
    - Diversity factor (Hunter's curve, updated v1.1.0)
    - OD to ID conversion using Type L copper specs (v1.1.0)
    - EPANET INP file generation
    - Publication-quality figures

Changes in v1.1.0:
    - Added copper Type L wall thickness lookup for OD to ID conversion
    - Updated diversity factor values (0.20 for 21-50 fixtures)
    - Pipe diameters now reported as inside diameter (ID)
        """)
        sys.exit(1)
    
    ifc_path = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else './output'
    
    if not Path(ifc_path).exists():
        print(f"Error: File not found: {ifc_path}")
        sys.exit(1)
    
    run_pipeline(ifc_path, output_dir)
