from dataclasses import dataclass
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.Basics.FoamFileGenerator import FoamFileGenerator
import os
import copy
import re
import math
import gmsh
import numpy

create_case = True
create_mesh = False

case_folder = "/mnt/d/Temp/pyfoam_test"
mesh_file = "/mnt/d/Temp/hydro.stp"
#mesh_file = "/mnt/d/Computer Aided Design/ahmed_body_something/ahmed_body.stp"

#-------------------------------------------------------------------------------
# Utilities
#-------------------------------------------------------------------------------

@dataclass
class FoamFile:
    header : dict
    data : dict
    footer: str = None

def make_header(version="2.0", format="ascii", file_class="dictionary", location=None, object=None):
    header = {
        "version": version,
        "format": format,
        "class": file_class,
    }
    if location:
        header["location"] = location
    if object:
        header["object"] = object
    return header

def foam_file_string(foam_file : FoamFile):
    header = copy.deepcopy(foam_file.header)
    if header and header["location"]:
        header["location"] = "\"" + str(header["location"]) + "\""
    file = FoamFileGenerator(foam_file.data, header=header)
    s = file.makeString()
    s = re.sub("_IncludeEtc_[0-9]+([^\\n]+);", "#includeEtc\\1", s)
    return s + foam_file.footer if foam_file.footer else s

def foam_uniform(value):
    if isinstance(value, list):
        return "uniform (" + " ".join([str(v) for v in value]) + ")"
    return "uniform " + str(value)

def foam_dimension(value):
    return "[" + " ".join([str(v) for v in value]) + "]"

def make_boundaries(field_classes, field_dimensions, field_internals, field_boundaries):
    files = []
    for field_name, field_class in field_classes.items():
        dimension = field_dimensions[field_name]
        internal = field_internals[field_name]
        boundary = field_boundaries[field_name]
        files.append(FoamFile(
            make_header(file_class=field_class, location="0", object=field_name),
            {
                "dimensions": dimension,
                "internalField": internal,
                "boundaryField": boundary,
            },
            "boundaryField {\n  #includeEtc \"caseDicts/setConstraintTypes\"\n}"
        ))
    return files

class OpenFoamCase:
    files = {}

    def add_file(self, foam_file : FoamFile, path=None):
        if not path:
            if foam_file.header and foam_file.header["location"] and foam_file.header["object"]:
                path = foam_file.header["location"] + "/" + foam_file.header["object"]
            else:
                raise ValueError()
        self.files[path] = foam_file

    def write(self, folder):
        if not os.path.exists(folder):
            raise ValueError("folder does not exist")
        for path, foam_file in self.files.items():
            full_path = os.path.join(folder, path)
            full_dir = os.path.dirname(full_path)
            if not os.path.exists(full_dir):
                os.makedirs(full_dir)
            with open(full_path, "w") as file:
                file.write(foam_file_string(foam_file))
        
    def print(self):
        for path, foam_file in self.files.items():
            print("//---------------------------------------")
            print("// " + path)
            print("//---------------------------------------\n")
            print(foam_file_string(foam_file))
            print("\n")

#-------------------------------------------------------------------------------
# System dicts
#-------------------------------------------------------------------------------

controlDict = FoamFile(
    make_header(location="system", object="controlDict"),
    {
        "application": "simpleFoam",
        "startTime": 0,
        "endTime": 200,
        "deltaT": 1,
        "startFrom": "startTime",
        "stopAt": "endTime",
        "writeControl": "timeStep",
        "writeInterval": 100,
        "purgeWrite": 0,
        "writeFormat": "ascii",
        "writePrecision": 6,
        "writeCompression": "off",
        "timeFormat": "general",
        "timePrecision": 6,
        "graphFormat": "raw",
        "runtTimeModifiable": "true",
    }
)

fvSolution = FoamFile(
    make_header(location="system", object="fvSolution"),
    {
        "solvers": {
            "\"(p|Phi)\"": {
                "preconditioner": "DIC",
                "solver": "PCG",
                "smoother": "GaussSeidel",
                "tolerance": 1e-06,
                "relTol": 0.01,
                "nCellsInCoursestLevel": 20,
                "minIter": 3
            },
            "\"(U|k|omega)\"": {
                "preconditioner": "DILU",
                "solver": "PBiCGStab",
                "smoother": "GaussSeidel",
                "tolerance": 1e-08,
                "relTol": 0.01,
                "nSweeps": 1,
                "minIter": 3
            }
        },
        "SIMPLE": {
            "nNonOrthogonalCorrectors": 0,
            "consistent": "yes",
        },
        "potentialFlow": {
            "nNonOrthogonalCorrectors": 10,
        },
        "relaxationFactors": {
            "equations": {
                "U": 0.4,
                "k": 0.35,
                "omega": 0.35,
                "p": 0.8,
            }
        }
    }
)

gradSchemeLimited = "cellLimited Gauss linear 1"

fvSchemes = FoamFile(
    make_header(location="system", object="fvSchemes"),
    {
        "ddtSchemes": {
            "default": "steadyState"
        },
        "gradSchemes": {
            "default": gradSchemeLimited,
            "grad(U)": gradSchemeLimited,
            "grad(k)": gradSchemeLimited,
            "grad(omega)": gradSchemeLimited,
        },
        "divSchemes": {
            "default": "none",
            "div(phi,U)": "Gauss linearUpwindV grad(U)",
            "div(phi,k)": "Gauss linearUpwind default",
            "div(phi,omega)": "Gauss linearUpwind default",
            "div((nuEff*dev2(T(grad(U)))))": "Gauss linear",
        },
        "laplacianSchemes": {
            "default": "Gauss linear limited 1",
        },
        "interpolationSchemes": {
            "default": "linear",
        },
        "snGradSchemes": {
            "default": "limited 1",
        },
        "wallDist": {
            "method": "meshWave",
        }
    }
)

fvConstraints = FoamFile(
    make_header(location="system", object="fvConstraints"),
    {
        "limitp": {
            "type": "limitPressure",
            "minFactor": 0.1,
            "maxFactor": 2.0,
        }
    }
)


#-------------------------------------------------------------------------------
# Constant dicts
#-------------------------------------------------------------------------------

momentumTransport = FoamFile(
    make_header(location="constant", object="momentumTransport"),
    {
        "simulationType": "RAS",
        "RAS": {
            "model": "kOmegaSST",
            "turbulence": "on",
            "printCoeffs": "on",
        }
    }
)

transportProperties = FoamFile(
    make_header(location="constant", object="transportProperties"),
    {
        "transportModel": "Newtonian",
        "nu": "[0 2 -1 0 0 0 0] 1.5e-05",
    }
)

#-------------------------------------------------------------------------------
# Initial and boundary conditions
#-------------------------------------------------------------------------------

pressure = foam_uniform(101325)
velocity = foam_uniform([0, -30, 0])
nut = foam_uniform(0)
turbulent_k = foam_uniform(0.24)
turbulent_omega = foam_uniform(1.78)

field_classes = {}
field_dimensions = {}
field_internals = {}
field_boundaries = {}

field_classes["U"] = "volVectorField"
field_classes["p"] = "volScalarField"
field_classes["nut"] = "volScalarField"
field_classes["k"] = "volScalarField"
field_classes["omega"] = "volScalarField"

field_dimensions["U"] = foam_dimension([0, 1, -1, 0, 0, 0, 0])
field_dimensions["p"] = foam_dimension([0, 2, -2, 0, 0, 0, 0])
field_dimensions["nut"] = foam_dimension([0, 2, -1, 0, 0, 0, 0])
field_dimensions["k"] = foam_dimension([0, 2, -2, 0, 0, 0, 0])
field_dimensions["omega"] = foam_dimension([0, 0, -1, 0, 0, 0, 0])

field_internals["U"] = velocity
field_internals["p"] = pressure
field_internals["nut"] = nut
field_internals["k"] = turbulent_k
field_internals["omega"] = turbulent_omega

field_boundaries["U"] = {}
field_boundaries["p"] = {}
field_boundaries["nut"] = {}
field_boundaries["k"] = {}
field_boundaries["omega"] = {}

field_boundaries["U"]["walls"] = { "type": "slip" }
field_boundaries["U"]["inlet"] = { "type": "fixedValue", "value": velocity }
field_boundaries["U"]["outlet"] = { "type": "inletOutlet", "value": velocity, "inletValue": foam_uniform([0,0,0]) }
#field_boundaries["U"]["ground"] = { "type": "fixedValue", "value": velocity }
field_boundaries["U"]["ground"] = { "type": "slip" }
field_boundaries["U"]["model"] = { "type": "noSlip" }

field_boundaries["p"]["walls"] = { "type": "slip" }
field_boundaries["p"]["inlet"] = { "type": "zeroGradient" }
field_boundaries["p"]["outlet"] = { "type": "fixedValue", "value": pressure }
#field_boundaries["p"]["ground"] = { "type": "zeroGradient" }
field_boundaries["p"]["ground"] = { "type": "slip" }
field_boundaries["p"]["model"] = { "type": "zeroGradient" }

field_boundaries["nut"]["walls"] = { "type": "calculated", "value": nut }
field_boundaries["nut"]["inlet"] = { "type": "calculated", "value": nut }
field_boundaries["nut"]["outlet"] = { "type": "calculated", "value": nut }
#field_boundaries["nut"]["ground"] = { "type": "nutkWallFunction", "value": nut }
field_boundaries["nut"]["ground"] = { "type": "calculated", "value": nut }
field_boundaries["nut"]["model"] = { "type": "nutkWallFunction", "value": nut }

field_boundaries["k"]["walls"] = { "type": "slip" }
field_boundaries["k"]["inlet"] = { "type": "fixedValue", "value": turbulent_k }
field_boundaries["k"]["outlet"] = { "type": "inletOutlet", "value": turbulent_k, "inletValue": turbulent_k }
#field_boundaries["k"]["ground"] = { "type": "kqRWallFunction", "value": turbulent_k }
field_boundaries["k"]["ground"] = { "type": "slip" }
field_boundaries["k"]["model"] = { "type": "kqRWallFunction", "value": turbulent_k }

field_boundaries["omega"]["walls"] = { "type": "slip" }
field_boundaries["omega"]["inlet"] = { "type": "fixedValue", "value": turbulent_omega }
field_boundaries["omega"]["outlet"] = { "type": "inletOutlet", "value": turbulent_omega, "inletValue": turbulent_omega }
#field_boundaries["omega"]["ground"] = { "type": "omegaWallFunction", "value": turbulent_omega }
field_boundaries["omega"]["ground"] = { "type": "slip", "value": turbulent_omega }
field_boundaries["omega"]["model"] = { "type": "omegaWallFunction", "value": turbulent_omega }

boundary_files = make_boundaries(field_classes, field_dimensions, field_internals, field_boundaries)

#-------------------------------------------------------------------------------
# Create case
#-------------------------------------------------------------------------------

# Create case
if create_case:
    case = OpenFoamCase()

    case.add_file(controlDict)
    case.add_file(fvSolution)
    case.add_file(fvSchemes)
    #case.add_file(fvConstraints)

    case.add_file(momentumTransport)
    case.add_file(transportProperties)

    for boundary_file in boundary_files:
        case.add_file(boundary_file)

    case.write(case_folder)

#-------------------------------------------------------------------------------
# Create mesh
#-------------------------------------------------------------------------------
if create_mesh:
    gmsh.initialize()
    gmsh.option.setString("Geometry.OCCTargetUnit", "M")
    gmsh.option.setNumber("Mesh.Algorithm", 1) # meshadapt
    gmsh.option.setNumber("Mesh.Algorithm3D", 1) # delaunay
    gmsh.option.setNumber("General.NumThreads", 8)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.option.setNumber("Mesh.AllowSwapAngle", 5)
    gmsh.model.add("Wind Tunnel")

    model_step = gmsh.model.occ.importShapes(mesh_file, True)
    model = model_step[0][1]

    model_aabb = gmsh.model.occ.getBoundingBox(3, model)
    model_x = (model_aabb[0], model_aabb[3])
    model_y = (model_aabb[1], model_aabb[4])
    model_z = (model_aabb[2], model_aabb[5])
    model_extents = [
        model_x[1] - model_x[0],
        model_y[1] - model_y[0],
        model_z[1] - model_z[0],
    ]

    model_min_frontal_extent = min(model_extents[0], model_extents[2])
    model_frontal_area = model_extents[0] * model_extents[2]
    model_max_extent = numpy.max(model_extents)

    enclosure_frontal_extent = math.sqrt(model_frontal_area / 2) * 16
    enclosure_x = (-enclosure_frontal_extent/2, enclosure_frontal_extent/2)
    enclosure_y = (model_z[0] - 1.5*enclosure_frontal_extent, model_z[1] + 0.75*enclosure_frontal_extent)
    enclosure_z = (0.0025, enclosure_frontal_extent/2)

    enclosure_extents = [
        enclosure_x[1] - enclosure_x[0],
        enclosure_y[1] - enclosure_y[0],
        enclosure_z[1] - enclosure_z[0],
    ]

    # Mesh sizes
    resolution = 0.05
    boundary_mesh_size = resolution * model_min_frontal_extent
    inner_boundary_size = 0.30 * boundary_mesh_size
    enable_inner_boundary = True
    near_mesh_size = 1.8 * boundary_mesh_size
    wake_mesh_size = 3 * boundary_mesh_size
    far_mesh_size = 6 * boundary_mesh_size
    base_mesh_size = far_mesh_size * 3

    near_region_offset = 0.7 * model_min_frontal_extent
    near_region_transition = 1.3 * near_region_offset
    far_region_offset = 2.1 * model_min_frontal_extent
    far_region_transition = 1.3 * far_region_offset

    wake_region_length = max(1.5*model_extents[1], 8.0*model_min_frontal_extent)
    far_region_length = 1.8*wake_region_length

    # subtract model from enclosure
    enclosure = gmsh.model.occ.addBox(enclosure_x[0], enclosure_y[0], enclosure_z[0], enclosure_extents[0], enclosure_extents[1], enclosure_extents[2])
    gmsh.model.occ.cut([(3, enclosure)], [(3, model)])
    gmsh.model.occ.synchronize()

    # set name for fluid volume
    volumes = gmsh.model.getEntities(dim=3)
    tag_fluid = gmsh.model.addPhysicalGroup(3, [volumes[0][1]])
    gmsh.model.setPhysicalName(3, tag_fluid, "fluid")

    # set name for boundaries
    surfaces = gmsh.model.getEntities(dim=2)
    tag_walls = []
    tag_models = []
    for surface in surfaces:
        tag_surface = surface[1]
        (mx, my, mz) = gmsh.model.occ.getCenterOfMass(2, tag_surface)
        if abs(my - enclosure_y[1]) < 1e-6:
            tag_inlet = tag_surface
        elif abs(my - enclosure_y[0]) < 1e-6:
            tag_outlet = tag_surface
        elif abs(mz - enclosure_z[0]) < 1e-6:
            tag_ground = tag_surface
        elif abs(mx - enclosure_x[0]) < 1e-6 or abs(mx - enclosure_x[1]) < 1e-6 or abs(mz - enclosure_z[1]) < 1e-6:
            tag_walls.append(tag_surface)
        else:
            tag_models.append(tag_surface)

    tag_inlet_group = gmsh.model.addPhysicalGroup(2, [tag_inlet])
    gmsh.model.setPhysicalName(2, tag_inlet_group, "inlet")

    tag_outlet_group = gmsh.model.addPhysicalGroup(2, [tag_outlet])
    gmsh.model.setPhysicalName(2, tag_outlet_group, "outlet")

    tag_ground_group = gmsh.model.addPhysicalGroup(2, [tag_ground])
    gmsh.model.setPhysicalName(2, tag_ground_group, "ground")

    tag_wall_group = gmsh.model.addPhysicalGroup(2, tag_walls)
    gmsh.model.setPhysicalName(2, tag_wall_group, "walls")

    tag_model_group = gmsh.model.addPhysicalGroup(2, tag_models)
    gmsh.model.setPhysicalName(2, tag_model_group, "model")

    # set mesh sizes
    distance_from_model = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(distance_from_model, "SurfacesList", tag_models)

    boundary_field = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(boundary_field, "IField", distance_from_model)
    gmsh.model.mesh.field.setNumber(boundary_field, "SizeMin", boundary_mesh_size)
    gmsh.model.mesh.field.setNumber(boundary_field, "SizeMax", near_mesh_size)
    gmsh.model.mesh.field.setNumber(boundary_field, "DistMin", 5*boundary_mesh_size)
    gmsh.model.mesh.field.setNumber(boundary_field, "DistMax", 5*near_mesh_size)
    gmsh.model.mesh.field.setNumber(boundary_field, "StopAtDistMax", 1)

    inner_boundary_field = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(inner_boundary_field, "IField", distance_from_model)
    gmsh.model.mesh.field.setNumber(inner_boundary_field, "SizeMin", inner_boundary_size)
    gmsh.model.mesh.field.setNumber(inner_boundary_field, "SizeMax", boundary_mesh_size)
    gmsh.model.mesh.field.setNumber(inner_boundary_field, "DistMin", 5*inner_boundary_size)
    gmsh.model.mesh.field.setNumber(inner_boundary_field, "DistMax", 5*boundary_mesh_size)
    gmsh.model.mesh.field.setNumber(inner_boundary_field, "StopAtDistMax", 1)

    near_region = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(near_region, "Thickness", far_region_transition)
    gmsh.model.mesh.field.setNumber(near_region, "VIn", near_mesh_size)
    gmsh.model.mesh.field.setNumber(near_region, "VOut", base_mesh_size)
    gmsh.model.mesh.field.setNumber(near_region, "XMin", model_x[0] - near_region_offset)
    gmsh.model.mesh.field.setNumber(near_region, "YMin", model_y[0] - near_region_offset)
    gmsh.model.mesh.field.setNumber(near_region, "ZMin", model_z[0] - near_region_offset)
    gmsh.model.mesh.field.setNumber(near_region, "XMax", model_x[1] + near_region_offset)
    gmsh.model.mesh.field.setNumber(near_region, "YMax", model_y[1] + near_region_offset)
    gmsh.model.mesh.field.setNumber(near_region, "ZMax", model_z[1] + near_region_offset)

    wake_region = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(wake_region, "Thickness", far_region_transition)
    gmsh.model.mesh.field.setNumber(wake_region, "VIn", wake_mesh_size)
    gmsh.model.mesh.field.setNumber(wake_region, "VOut", base_mesh_size)
    gmsh.model.mesh.field.setNumber(wake_region, "XMin", model_x[0] - near_region_offset)
    gmsh.model.mesh.field.setNumber(wake_region, "YMin", model_y[0] - wake_region_length)
    gmsh.model.mesh.field.setNumber(wake_region, "ZMin", model_z[0] - near_region_offset)
    gmsh.model.mesh.field.setNumber(wake_region, "XMax", model_x[1] + near_region_offset)
    gmsh.model.mesh.field.setNumber(wake_region, "YMax", model_y[0])
    gmsh.model.mesh.field.setNumber(wake_region, "ZMax", model_z[1] + near_region_offset)

    far_region = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(far_region, "Thickness", far_region_transition)
    gmsh.model.mesh.field.setNumber(far_region, "VIn", far_mesh_size)
    gmsh.model.mesh.field.setNumber(far_region, "VOut", base_mesh_size)
    gmsh.model.mesh.field.setNumber(far_region, "XMin", model_x[0] - far_region_offset)
    gmsh.model.mesh.field.setNumber(far_region, "YMin", model_y[0] - far_region_length)
    gmsh.model.mesh.field.setNumber(far_region, "ZMin", model_z[0] - far_region_offset)
    gmsh.model.mesh.field.setNumber(far_region, "XMax", model_x[1] + far_region_offset)
    gmsh.model.mesh.field.setNumber(far_region, "YMax", model_y[1] + 0.4*far_region_length)
    gmsh.model.mesh.field.setNumber(far_region, "ZMax", model_z[1] + far_region_offset)

    min_field = gmsh.model.mesh.field.add("Min")
    combined_fields = [boundary_field, near_region, wake_region, far_region]
    if enable_inner_boundary:
        combined_fields.append(inner_boundary_field)
    gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", combined_fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(min_field)

    for surface in gmsh.model.getEntities(2):
        gmsh.model.mesh.setSmoothing(2, surface[1], 100)

    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)

    gmsh.write(os.path.join(case_folder, "mesh.neu"))

    gmsh.finalize()
