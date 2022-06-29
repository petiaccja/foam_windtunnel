from dataclasses import dataclass
import pathlib
import math
import gmsh
from typing import Tuple


@dataclass
class Transform:
    pos: Tuple[float, float, float]
    rot: Tuple[Tuple[float, float, float], float]


def create_mesh(output_file: pathlib.Path, model: pathlib.Path, transform: Transform, resolution: float):
    gmsh.initialize()
    gmsh.option.setString("Geometry.OCCTargetUnit", "M")
    gmsh.option.setNumber("Mesh.Algorithm", 1) # meshadapt
    gmsh.option.setNumber("Mesh.Algorithm3D", 1) # delaunay
    gmsh.option.setNumber("General.NumThreads", 8)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.option.setNumber("Mesh.AllowSwapAngle", 5)
    gmsh.model.add("Wind Tunnel")

    model_step = gmsh.model.occ.importShapes(str(model), True)
    gmsh.model.occ.rotate(model_step, 0, 0, 0, transform.rot[0][0], transform.rot[0][1], transform.rot[0][2], transform.rot[1])
    gmsh.model.occ.translate(model_step, transform.pos[0], transform.pos[1], transform.pos[2])
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

    enclosure_frontal_extent = max([2*model_extents[0], 2*model_extents[2], math.sqrt(model_frontal_area / 2) * 36])
    enclosure_x = (-enclosure_frontal_extent/2, enclosure_frontal_extent/2)
    enclosure_y = (model_z[0] - 1.5*enclosure_frontal_extent, model_z[1] + 0.75*enclosure_frontal_extent)
    enclosure_z = (0, enclosure_frontal_extent/2)

    enclosure_extents = [
        enclosure_x[1] - enclosure_x[0],
        enclosure_y[1] - enclosure_y[0],
        enclosure_z[1] - enclosure_z[0],
    ]

    # Mesh sizes
    boundary_mesh_size = resolution * model_min_frontal_extent
    inner_boundary_size = 0.30 * boundary_mesh_size
    enable_inner_boundary = True
    near_mesh_size = 1.8 * boundary_mesh_size
    wake_mesh_size = 3 * boundary_mesh_size
    far_mesh_size = 6 * boundary_mesh_size
    base_mesh_size = far_mesh_size * 3

    near_region_offset = 0.7 * model_min_frontal_extent
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

    edges = gmsh.model.getEntities(dim=1)
    tag_model_edges = []
    for edge in edges:
        tag_edge = edge[1]
        (mx, my, mz) = gmsh.model.occ.getCenterOfMass(1, tag_edge)
        if enclosure_x[0] < mx and mx < enclosure_x[1]\
            and enclosure_y[0] < my and my < enclosure_y[1]\
            and enclosure_z[0] < mz and mz < enclosure_z[1]:
            tag_model_edges.append(tag_edge)

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

    gmsh.write(str(output_file))

    gmsh.finalize()