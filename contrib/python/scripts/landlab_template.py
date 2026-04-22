"""
Template class used for creating a LandLab script that is then used in the ASPECT-Landlab 
coupling. This can be imported and the functions can be modified to tailor the needs of the 
user for productions models or simple tests. 
"""

from mpi4py import MPI
import inspect
import json
import os
import numpy as np
import landlab
from landlab.io.native_landlab import save_grid, load_grid
from landlab.io.legacy_vtk import write_legacy_vtk
from landlab.components import AdvectionSolverTVD


class LandLabTemplate:
    """Template class for LandLab scripts used in ASPECT-Landlab coupling."""

    # Functions that actually get called within ASPECT. These functions are checked to ensure
    # that the signatures of these function are not modified when instantiating this class.
    _signature_checked_methods = (
        "initialize",
        "finalize",
        "set_mesh_information",
        "initialize_landlab_components",
        "update_until",
        "checkpoint_model_grid",
        "load_model_grid",
        "get_initial_topography",
        "write_output",
        "get_grid_x",
        "get_grid_y",
        "get_grid_z",
    )

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        for method_name in LandLabTemplate._signature_checked_methods:
            # Only validate methods that the subclass overrides directly.
            if method_name not in cls.__dict__:
                continue

            base_method = getattr(LandLabTemplate, method_name)
            sub_method = cls.__dict__[method_name]

            base_signature = inspect.signature(base_method)
            sub_signature = inspect.signature(sub_method)

            if sub_signature != base_signature:
                raise TypeError(
                    f"Invalid override for '{method_name}' in class "
                    f"'{cls.__name__}': expected signature {base_signature}, "
                    f"got {sub_signature}."
                )

    def __init__(self):
        self.current_time = 0.0
        self.comm = None

        self.model_grid = None
        self.elevation  = None

        self.horizontal_velocity         = None
        self.horizontal_surface_advector = None


        self.s2yr = 60 * 60 * 24 * 365.25
        self.timestep = 0
        self.vtks = []

    def initialize(self, comm_handle):
        """Called once by ASPECT at startup."""
        if comm_handle is not None:
            self.comm = MPI.Comm.f2py(comm_handle)
            rank = self.comm.Get_rank()
            size = self.comm.Get_size()
            print(f"Python: Hello from Rank {rank} of {size}")

            globalsum = self.comm.allreduce(1, op=MPI.SUM)
            if self.comm.rank == 0:
                print(f"\tPython: testing communication; sum {globalsum}")
        else:
            print("Python: running sequentially!")

    def finalize(self):
        """Currently unused."""
        pass

    # ---------------------------------------------------------------------------
    # The following functions will need to be modified based on the needs of the user:
    # - set_mesh_information: create the model grid and initialize the elevation field.
    #   Here the user can define the geometry and resolution of the landlab mesh.
    #
    # - initialize_landlab_components: initialize any LandLab components. Here the user
    #   chooses which Landlab components will influence the evolution of the topograhy, and
    #   initializes them with the necessary parameters.
    #
    # - update_until: advance the LandLab model through time based on the ASPECT timestep.
    #   Here the user needs to define how the Landlab components influence the evolution of
    #   the topography, as well as how fields from ASPECT (composition, temperature etc.,)
    #   influence the Landlab model.
    # ---------------------------------------------------------------------------
    def set_mesh_information(self, grid_dictionary):
        """
        Create the Landlab model grid, initialize the elevation field, and any other grid fields needed for
        the Landlab model. Call the function initialize_landlab_components() at the end of this function to 
        initialize the Landlab components after the grid and fields have been created.
        ** Derived classes must override this function in custom scripts. **
        """
        raise NotImplementedError("LandLabTemplate.set_mesh_information() must be overridden in your custom script." \
                                  "See the doc string in the template file for information on how to write this function.")



    def initialize_landlab_components(self):
        """
        Initialize the Landlab components that will be used to evolve the topography over time.
        ** Derived classes must override this function in custom scripts. **
        """
        raise NotImplementedError("LandLabTemplate.initialize_landlab_components() must be overridden in your custom script." \
                                  "See the doc string in the template file for information on how to write this function.")



    def update_until(self, end_time, ASPECT_dim, ASPECT_fields_at_Landlab_nodes_dict):
        """
        Run the Landlab model for the duration of the ASPECT timestep and return the change in topography at each node.
        This function is where the LandLab components will be called. The change in topography is used to determine the 
        surface velocity of the ASPECT mesh.
        ** Derived classes must override this function in custom scripts. **

        Parameters:
        - end_time: the time to advance the Landlab model to.

        - ASPECT_dim: the dimensionality of the ASPECT model (2 or 3).

        - ASPECT_fields_at_Landlab_nodes_dict: a dictionary mapping ASPECT variable names to their values at each node on
                                               the LandLab mesh. This dictionary will include entries for "x velocity" and
                                               "y velocity" (and "z velocity" for 3D models), as well as the compositional
                                                fields, and the pressure and temperature.
        """

        raise NotImplementedError("LandLabTemplate.update_until() must be overridden in your custom script." \
                                  "See the doc string in the template file for information on how to write this function.")



    # ---------------------------------------------------------------------------
    # The rest of these functions likely will not require any modification.
    # ---------------------------------------------------------------------------
    def determine_uplift_velocity(self, ASPECT_dim, dict_variable_name_to_value_in_nodes):
        """Determine uplift velocity from ASPECT variables."""
        if ASPECT_dim == 2:
            slice_y_velocity = dict_variable_name_to_value_in_nodes["y velocity"]

            vertical_velocity = np.zeros(self.model_grid.number_of_nodes)
            unique_x_values = np.unique(self.model_grid.x_of_node)
            for x in unique_x_values:
                vertical_velocity[self.model_grid.x_of_node == x] = slice_y_velocity[unique_x_values == x]
        elif ASPECT_dim == 3:
            vertical_velocity = dict_variable_name_to_value_in_nodes["z velocity"]

        return vertical_velocity
    
    def determine_horizontal_velocity(self, ASPECT_dim, dict_variable_name_to_value_in_nodes):
        """Determine horizontal velocity from ASPECT variables."""
        x_velocity = dict_variable_name_to_value_in_nodes["x velocity"]

        if ASPECT_dim == 2:
            projected_x_velocity = np.zeros(self.model_grid.number_of_nodes)
            unique_x_values   = np.unique(self.model_grid.x_of_node)
            for x in unique_x_values:
                projected_x_velocity[self.model_grid.x_of_node == x] = x_velocity[unique_x_values == x]
            
            x_vel_at_links = self.model_grid.map_mean_of_link_nodes_to_link(projected_x_velocity)
            y_vel_at_links = self.model_grid.map_mean_of_link_nodes_to_link(np.zeros(self.model_grid.number_of_nodes)) # y velocity is zero since the ASPECT model is 2D.

        elif ASPECT_dim == 3:
            y_velocity = dict_variable_name_to_value_in_nodes["y velocity"]

            x_vel_at_links = self.model_grid.map_mean_of_link_nodes_to_link(x_velocity)
            y_vel_at_links = self.model_grid.map_mean_of_link_nodes_to_link(y_velocity) # y velocity is zero since the ASPECT model is 2D.

        self.horizontal_velocity[self.model_grid.horizontal_links] = x_vel_at_links[self.model_grid.horizontal_links]
        self.horizontal_velocity[self.model_grid.vertical_links]   = y_vel_at_links[self.model_grid.vertical_links]

        return self.horizontal_velocity

    def dimensional_deposition_erosion(self, ASPECT_dim, deposition_erosion):
        if ASPECT_dim == 2:
            deposition_erosion_2d = np.zeros(len(np.unique(self.model_grid.x_of_node)))
            unique_x_values = np.unique(self.model_grid.x_of_node)

            for x in unique_x_values:
                deposition_erosion_2d[unique_x_values == x] = np.average(deposition_erosion[self.model_grid.x_of_node == x])
            return deposition_erosion_2d
        
        elif ASPECT_dim == 3:
            return deposition_erosion

    def checkpoint_model_grid(self):
        filename = "./landlab_model_grid_checkpoint.grid"
        save_grid(self.model_grid, filename, clobber=True)
        pass

    def load_model_grid(self):
        filename = "./landlab_model_grid_checkpoint.grid"
        self.model_grid = load_grid(filename)
        self.elevation = self.model_grid.at_node["topographic__elevation"]
        self.initialize_landlab_components(None)
        pass

    def get_initial_topography(self, ASPECT_dim):
        if ASPECT_dim == 2:
            return self.elevation[self.model_grid.y_of_node == 0]
        elif ASPECT_dim == 3:
            return self.elevation
        
    def write_output(self, postprocess_dictionary):
        step = postprocess_dictionary["ASPECT timestep"]
        time = postprocess_dictionary["ASPECT time"]
        output_directory = postprocess_dictionary["ASPECT output directory"]

        output_frequency = 10
        if step % output_frequency != 0:
            return

        filename = f"{output_directory}/landlab_{str(step).zfill(3)}.vtk"
        write_legacy_vtk(path=filename, grid=self.model_grid, clobber=True)
        self.vtks.append((time, filename))

        with open(f"{output_directory}/landlab.vtk.series", "w") as f:
            series = {
                "file-series-version": "1.0",
                "files": [
                    {"name": os.path.basename(vtk_name), "time": vtk_time}
                    for vtk_time, vtk_name in self.vtks
                ],
            }
            json.dump(series, f, indent=2)

        pass

    def get_grid_x(self, ASPECT_dim):
        if ASPECT_dim == 2:
            return np.unique(self.model_grid.x_of_node)
        elif ASPECT_dim == 3:
            return self.model_grid.x_of_node
    
    def get_grid_y(self, ASPECT_dim):
        if ASPECT_dim == 2:
            return np.unique(self.model_grid.y_of_node)
        elif ASPECT_dim == 3:
            return self.model_grid.y_of_node
        
    def get_grid_z(self, ASPECT_dim):
        if ASPECT_dim == 3:
            return self.model_grid.z_of_node
        else:
            raise ValueError("get_grid_z is only applicable for 3D ASPECT models.")

    @staticmethod
    def export_aspect_callbacks(model, namespace):
        """
        Export the functions of the derived class as module-level functions in the provided namespace. 
        By using this function, the user can simply define a class that inherits from LandLabTemplate, 
        override any desired functions, and then call this export function to make the functions available 
        to ASPECT without needing to always write boilerplate functions.

        Parameters:
        - model: an instance of a class that inherits from LandLabTemplate.
        - namespace: the namespace to which the functions should be added. This is typically done by
                     passing in globals() from the script that defines the derived class.

        Example usage:
        model = MyAspectLandlabModel()
        model.export_aspect_callbacks(model, globals())
        """

        # Name of all the functions that ASPECT requires.
        callback_names = (
            "initialize",
            "finalize",
            "update_until",
            "set_mesh_information",
            "get_grid_x",
            "get_grid_y",
            "get_initial_topography",
            "checkpoint_model_grid",
            "load_model_grid",
            "write_output",
        )

        namespace["model"] = model
        for name in callback_names:
            namespace[name] = getattr(model, name)
