from sys import argv
from foam_case import create_openfoam_case
from mesher import create_mesh, Transform
import pathlib
from PyFoam.Applications.PlotRunner import PlotRunner
from PyFoam.Applications.Execute import Execute
import re
import shutil


case_dir = pathlib.Path("/mnt/d/Computer_Simulations/hydro_aero/case")
run_dir = pathlib.Path("/mnt/d/Computer_Simulations/hydro_aero/run")
mesh_file = "/mnt/d/Computer_Simulations/hydro_aero/model.stp"

iterations = 500
parallel = 24
pressure = 101325
velocity = [0, -30, 0]
mesh_resolution = 0.1
transform = Transform((0, 0, -0.002), (0, 0, 0))


def write_case():
    openfoam_case = create_openfoam_case(iterations=iterations,
                                         pressure=pressure,
                                         velocity=velocity,
                                         parallel=parallel)
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


def prepare(create_case: bool, create_mesh: bool):
    if create_case:
        write_case()
    if create_mesh:
        write_mesh()
        Execute(args=["gambitToFoam", "-case", str(case_dir), str(case_dir / "mesh.neu")])
        update_boundary_patches()
        Execute(args=["polyDualMesh", "-case", str(case_dir), "-overwrite", "70"])


def run():
    Execute(args=["decomposePar", "-case", str(case_dir)])
    PlotRunner(args=["--procnr={}".format(parallel), "simpleFoam", "-case", str(case_dir)])


if __name__ == "__main__":
    #shutil.rmtree(case_dir / "constant")
    #shutil.rmtree(case_dir / "0")
    prepare(create_case=True, create_mesh=False)
    run()