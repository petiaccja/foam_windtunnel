import math
from sys import argv
from foam_case import create_openfoam_case
from mesher import create_mesh, Transform
import pathlib
from PyFoam.Applications.PlotRunner import PlotRunner
from PyFoam.Applications.Execute import Execute
import re
import shutil
import matplotlib.pyplot as plt


case_dir = pathlib.Path("/mnt/d/Computer_Simulations/hydro_aero/case")
run_dir = pathlib.Path("/mnt/d/Computer_Simulations/hydro_aero/run")
mesh_file = "/mnt/d/Computer_Simulations/hydro_aero/model.stp"

iterations = 600
parallel = 24
pressure = 0
velocity = [0, -30, 0]
cofr = [0.000, -0.100, 0.025]
mesh_resolution = 0.04
transform = Transform((0, 0, -0.001), ((1, 0, 0), math.radians(0)))


def write_case():
    openfoam_case = create_openfoam_case(iterations=iterations,
                                         pressure=pressure,
                                         velocity=velocity,
                                         parallel=parallel,
                                         cofr = cofr)
    openfoam_case.write(case_dir)


def write_mesh():
    create_mesh(output_file=(case_dir / "mesh.neu"), 
                model=mesh_file,
                transform=transform,
                resolution=mesh_resolution)


def update_boundary_patches():
    file = case_dir / "constant" / "polyMesh" / "boundary"
    text = file.read_text()
    new = re.sub(r"type(\s+)patch;", r"type\1wall;", text)
    file.write_text(new)


def clean():
    patterns = [
        "0",
        "constant",
        "system",
        "mesh.*",
        "postProcessing",
        "processor*",
        "PyFoam*",
        "Gnuplot*"
    ]
    for pattern in patterns:
        for subdir in case_dir.glob(pattern):
            print(f"Removing: {subdir}")
            if subdir.is_dir():
                shutil.rmtree(subdir)
            else:
                subdir.unlink()


def prepare(create_case: bool, create_mesh: bool):
    if create_case:
        write_case()
    if create_mesh:
        write_mesh()
        Execute(args=["gambitToFoam", "-case", str(case_dir), str(case_dir / "mesh.neu")])
        update_boundary_patches()
        Execute(args=["polyDualMesh", "-case", str(case_dir), "-overwrite", "70"])
        Execute(args=["decomposePar", "-case", str(case_dir)])


def run():
    PlotRunner(args=["--procnr={}".format(parallel), "simpleFoam", "-case", str(case_dir)])


def post_process():
    forces_file = case_dir / "postProcessing" / "modelForces" / "0" / "forces.dat"
    forces_text = forces_file.read_text()
    drags = []
    lifts = []
    torques = []
    for line in forces_text.splitlines():
        if line[0] == "#":
            continue
        items = list(filter(lambda s : len(s) > 0, line.replace("(", " ").replace(")", " ").replace("\t", " ").split(" ")))
        iter = int(items[0])
        pressure = [float(items[1]), float(items[2]), float(items[3])]
        viscous = [float(items[4]), float(items[5]), float(items[6])]
        pressure_moment = [float(items[7]), float(items[8]), float(items[9])]
        viscous_moment = [float(items[10]), float(items[11]), float(items[12])]
        drag = -(pressure[1] + viscous[1])
        lift = (pressure[2] + viscous[2])
        torque = pressure_moment[0] + viscous_moment[0]
        drags.append(drag)
        lifts.append(lift)
        torques.append(torque)
    return drags, lifts, torques


if __name__ == "__main__":
    clean()
    prepare(create_case=True, create_mesh=True)
    run()
    drags, lifts, torques = post_process()
    print("Drag: {} N".format(drags[-1]))
    print("Lift: {} N".format(lifts[-1]))
    print("Torque X: {} Nm".format(torques[-1]))