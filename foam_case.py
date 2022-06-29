from dataclasses import dataclass
import pathlib
from typing import Sequence
from PyFoam.Basics.FoamFileGenerator import FoamFileGenerator
import os
import copy
import re

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

#-------------------------------------------------------------------------------
# System dicts
#-------------------------------------------------------------------------------

def create_control_dict(iterations, pressure, cofr, rho):
    return FoamFile(
        make_header(location="system", object="controlDict"),
        {
            "application": "simpleFoam",
            "startTime": 0,
            "endTime": iterations,
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
            "functions": {
                "modelForces": {
                    "type": "forces",
                    "libs": ['"libforces.so"'],
                    "patches": ["model"],
                    "pRef": pressure,
                    "CofR": cofr,
                    "rho": "rhoInf",
                    "rhoInf": rho,
                }
            }
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
                "p": 0.7,
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

def create_decompose_dict(parallel: int):
    return FoamFile(
        make_header(location="system", object="decomposeParDict"),
        {
            "numberOfSubdomains": parallel,
            "method": "scotch",
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

def create_boundary_files(pressure: float, velocity: Sequence[float], nut: float=0.0, turbulent_k: float=0.24, turbulent_omega: float=1.78):
    pressure_str = foam_uniform(pressure)
    velocity_str = foam_uniform(velocity)
    nut_str = foam_uniform(nut)
    turbulent_k_str = foam_uniform(turbulent_k)
    turbulent_omega_str = foam_uniform(turbulent_omega)

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

    field_internals["U"] = velocity_str
    field_internals["p"] = pressure_str
    field_internals["nut"] = nut_str
    field_internals["k"] = turbulent_k_str
    field_internals["omega"] = turbulent_omega_str

    field_boundaries["U"] = {}
    field_boundaries["p"] = {}
    field_boundaries["nut"] = {}
    field_boundaries["k"] = {}
    field_boundaries["omega"] = {}

    field_boundaries["U"]["walls"] = { "type": "slip" }
    field_boundaries["U"]["inlet"] = { "type": "fixedValue", "value": velocity_str }
    field_boundaries["U"]["outlet"] = { "type": "inletOutlet", "value": velocity_str, "inletValue": foam_uniform([0,0,0]) }
    field_boundaries["U"]["ground"] = { "type": "slip" }
    field_boundaries["U"]["model"] = { "type": "noSlip" }

    field_boundaries["p"]["walls"] = { "type": "slip" }
    field_boundaries["p"]["inlet"] = { "type": "zeroGradient" }
    field_boundaries["p"]["outlet"] = { "type": "fixedValue", "value": pressure_str }
    field_boundaries["p"]["ground"] = { "type": "slip" }
    field_boundaries["p"]["model"] = { "type": "zeroGradient" }

    field_boundaries["nut"]["walls"] = { "type": "calculated", "value": nut_str }
    field_boundaries["nut"]["inlet"] = { "type": "calculated", "value": nut_str }
    field_boundaries["nut"]["outlet"] = { "type": "calculated", "value": nut_str }
    field_boundaries["nut"]["ground"] = { "type": "calculated", "value": nut_str }
    field_boundaries["nut"]["model"] = { "type": "nutkWallFunction", "value": nut_str }

    field_boundaries["k"]["walls"] = { "type": "slip" }
    field_boundaries["k"]["inlet"] = { "type": "fixedValue", "value": turbulent_k_str }
    field_boundaries["k"]["outlet"] = { "type": "inletOutlet", "value": turbulent_k_str, "inletValue": turbulent_k_str }
    field_boundaries["k"]["ground"] = { "type": "slip" }
    field_boundaries["k"]["model"] = { "type": "kqRWallFunction", "value": turbulent_k_str }

    field_boundaries["omega"]["walls"] = { "type": "slip" }
    field_boundaries["omega"]["inlet"] = { "type": "fixedValue", "value": turbulent_omega_str }
    field_boundaries["omega"]["outlet"] = { "type": "inletOutlet", "value": turbulent_omega_str, "inletValue": turbulent_omega_str }
    field_boundaries["omega"]["ground"] = { "type": "slip", "value": turbulent_omega_str }
    field_boundaries["omega"]["model"] = { "type": "omegaWallFunction", "value": turbulent_omega_str }

    boundary_files = make_boundaries(field_classes, field_dimensions, field_internals, field_boundaries)
    return boundary_files


#-------------------------------------------------------------------------------
# Create case
#-------------------------------------------------------------------------------


class OpenFoamCase:
    files = {}

    def add_file(self, foam_file : FoamFile, path=None):
        if not path:
            if foam_file.header and foam_file.header["location"] and foam_file.header["object"]:
                path = foam_file.header["location"] + "/" + foam_file.header["object"]
            else:
                raise ValueError()
        self.files[path] = foam_file

    def write(self, folder: pathlib.Path):
        if not folder.exists():
            raise ValueError("Root folder must exist!")
        for path, foam_file in self.files.items():
            full_path = folder / path
            full_dir = full_path.parent
            full_dir.mkdir(parents=True, exist_ok=True)
            if not os.path.exists(full_dir):
                os.makedirs(full_dir)
            full_path.write_text(foam_file_string(foam_file))
        
    def print(self):
        for path, foam_file in self.files.items():
            print("//---------------------------------------")
            print("// " + path)
            print("//---------------------------------------\n")
            print(foam_file_string(foam_file))
            print("\n")


def create_openfoam_case(iterations: int, pressure: float, velocity: Sequence[float], parallel: int = 1, cofr=[0, 0, 0], rho=1.2041):
    case = OpenFoamCase()

    case.add_file(create_control_dict(iterations, pressure=pressure, cofr=cofr, rho=rho))
    case.add_file(fvSolution)
    case.add_file(fvSchemes)
    case.add_file(create_decompose_dict(parallel))

    case.add_file(momentumTransport)
    case.add_file(transportProperties)

    for boundary_file in create_boundary_files(pressure, velocity):
        case.add_file(boundary_file)

    return case
