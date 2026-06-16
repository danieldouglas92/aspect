
from mpi4py import MPI
import importlib.util
import os
import sys
import numpy as np
import landlab
from landlab.components import LinearDiffuser, AdvectionSolverTVD

# ---------------------------------------------------------------------------
# Import LandLabTemplate base class, which is located in contrib/python/scripts/.
# ---------------------------------------------------------------------------
from landlab_template import LandLabTemplate


# ---------------------------------------------------------------------------
# Custom model — override only the methods you want to change.
# ---------------------------------------------------------------------------
class MyAspectLandlabModel(LandLabTemplate):

    def update_until(self, end_time, ASPECT_dim, ASPECT_fields_at_Landlab_nodes_dict):
        dt = end_time - self.current_time
        self.timestep += 1

        deposition_erosion = np.zeros(self.model_grid.number_of_nodes)

        vertical_velocity                = self.determine_uplift_velocity(ASPECT_dim, ASPECT_fields_at_Landlab_nodes_dict)
        self.horizontal_velocity         = self.determine_horizontal_velocity(ASPECT_dim, ASPECT_fields_at_Landlab_nodes_dict)
        self.horizontal_surface_advector = AdvectionSolverTVD(self.model_grid, fields_to_advect=self.elevation)

        # Extract compositions along the ASPECT surface.
        slice_weak_composition   = ASPECT_fields_at_Landlab_nodes_dict["weak"]
        slice_strong_composition = ASPECT_fields_at_Landlab_nodes_dict["strong"]

        # Project composition values from y=0 to all Landlab nodes.
        strong_composition = np.zeros(self.model_grid.number_of_nodes)
        weak_composition   = np.zeros(self.model_grid.number_of_nodes)

        unique_x_values = np.unique(self.model_grid.x_of_node)
        for x in unique_x_values:
            strong_composition[self.model_grid.x_of_node == x] = slice_strong_composition[unique_x_values == x]
            weak_composition[self.model_grid.x_of_node == x]   = slice_weak_composition[unique_x_values == x]

        # Modify diffusivity based on composition.
        self.Diffusivity[strong_composition >= 0.5] = 1e-10 / self.s2yr
        self.Diffusivity[weak_composition >= 0.5]   = 10 / self.s2yr

        if dt > 0:
            n_substeps = 10
            sub_dt = dt / n_substeps
            for _ in range(n_substeps):
                elevation_before = self.elevation.copy()
                self.linear_diffuser.run_one_step(sub_dt)
                self.elevation += vertical_velocity * sub_dt
                deposition_erosion += self.elevation - elevation_before

        self.current_time = end_time
        print("Max elevation:", np.max(self.elevation), "Min elevation:", np.min(self.elevation))

        # Return the change in topography along y=0, where the ASPECT mesh is located.
        dimensional_deposition_erosion = self.dimensional_deposition_erosion(ASPECT_dim, deposition_erosion)
        return dimensional_deposition_erosion

    def set_mesh_information(self, grid_dictionary):
        if self.model_grid is not None:
            return

        print("* Creating RasterModelGrid ...")
        x_extent = 100e3
        y_extent = 100e3
        spacing  = 500.0

        nrows = int(y_extent / spacing) + 1
        ncols = int(x_extent / spacing) + 1

        self.model_grid = landlab.RasterModelGrid(
            (nrows, ncols),
            xy_spacing=(spacing, spacing),
            xy_of_lower_left=(0, -y_extent / 2),
        )

        print("* Creating topographic elevation ...")
        self.elevation = self.model_grid.add_zeros("topographic__elevation", at="node")
        self.horizontal_velocity = self.model_grid.add_zeros("advection__velocity", at="link")

        # Triangular mountain in the centre of the domain.
        topo_height  = 20e3
        left_x_arr   = np.array([25e3, 50e3])
        left_y_arr   = np.array([0.0, topo_height])
        right_x_arr  = np.array([50e3, 75e3])
        right_y_arr  = np.array([topo_height, 0.0])

        left_m,  left_b  = np.polyfit(left_x_arr,  left_y_arr,  deg=1)
        right_m, right_b = np.polyfit(right_x_arr, right_y_arr, deg=1)

        self.elevation[self.model_grid.x_of_node <= 50e3] = (
            left_m * self.model_grid.x_of_node[self.model_grid.x_of_node <= 50e3] + left_b
        )
        self.elevation[self.model_grid.x_of_node > 50e3] = (
            right_m * self.model_grid.x_of_node[self.model_grid.x_of_node > 50e3] + right_b
        )
        self.elevation[self.model_grid.x_of_node < np.min(left_x_arr)]  = 0.0
        self.elevation[self.model_grid.x_of_node > np.max(right_x_arr)] = 0.0

        self.model_grid.set_closed_boundaries_at_grid_edges(
            right_is_closed=True,
            left_is_closed=True,
            top_is_closed=True,
            bottom_is_closed=True,
        )
        print("\tnumber of nodes:", self.model_grid.number_of_nodes)

        self.initialize_landlab_components()
        print("* Done")
        return self.model_grid

    def initialize_landlab_components(self):
        D = 1e-10 / self.s2yr
        self.Diffusivity = self.model_grid.add_zeros("linear_diffusivity", at="node")
        self.Diffusivity += D
        self.linear_diffuser = LinearDiffuser(self.model_grid, linear_diffusivity=self.Diffusivity)


# ---------------------------------------------------------------------------
# Module-level instance — ASPECT calls these functions by name.
# ---------------------------------------------------------------------------
model = MyAspectLandlabModel()
model.export_aspect_callbacks(model, globals())

