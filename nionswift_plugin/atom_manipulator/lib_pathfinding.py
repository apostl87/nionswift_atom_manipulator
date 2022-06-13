"""
Pathfinding algorithm and library.
- Determines the shortest total path from all sources to all targets.
- Based on pathfinding algorithms published in
  [1] Hart, P., Nilsson, N. & Raphael, B. IEEE Trans. Syst. Sci. Cybern. 4, 100–107 (1968).
  [2] Kalff, F. E., et al. Nature nanotechnology 11.11, 926-929 (2016).
"""

import gettext
import threading
import numpy as np

import copy
import time
import logging

# Nion libraries
from nion.utils import Geometry

# Custom libraries
from .classes import atoms_and_bonds as aab, paths
from . import lib_utils

_ = gettext.gettext


# Main pathfinding function.
def find_paths(manipulator, auto_manipulate=False):

    if (manipulator.sites == []) and not auto_manipulate:
            logging.info(lib_utils.log_message("No sites found. Pathfinder aborted."))
            return
    if manipulator.t5 is not None and manipulator.t5.is_alive():
            logging.info(lib_utils.log_message("Pathfinder still working. Wait until finished."))
            return 
    
    def do_this():
        stop = False
        while not stop and not (auto_manipulate and manipulator.manipulation_module.stop_auto_manipulate_event.is_set()):
            if auto_manipulate:
                if not manipulator.structure_recognition_module.rdy.wait(0.1):
                    continue
                else:
                    manipulator.structure_recognition_module.rdy.clear()
            else:
                stop = True

            # Aliases.
            pdi = manipulator.processed_data_item
            
            # Set bonds.
            t = time.time()

            logging.info(lib_utils.log_message("Setting bonds..."))
            max_bond_length_px = manipulator.source_xdata.dimensional_calibrations[0].convert_from_calibrated_size(
                    manipulator.pathfinding_module.max_bond_length/10)
            if manipulator.simulation_mode: # Fix for wrong conversion in nionswift-usim fork
                max_bond_length_px *= 10 # Wrong conversion im usim fork
            manipulator.bonds = aab.Bonds(manipulator.sites, max_bond_length_px)
            
            t = time.time()-t
            logging.info(lib_utils.log_message(f"Setting bonds finished after {t:.5f} seconds"))
            
            # Call pathfinder.
            t = time.time()

            logging.info(lib_utils.log_message("Pathfinder called."))
            try:
                manipulator.paths = paths.Paths(manipulator.sources, manipulator.targets)
            except ValueError as e:
                print(e)
                return
            else:
                manipulator.paths.determine_paths_nooverlap()
            
            # Display paths.
            while not manipulator.rdy_init_pdi.wait(1) or not manipulator.rdy_update_pdi.wait(1):
                pass
            tmp_image = copy.copy(pdi.data)
            lib_utils.plot_paths(tmp_image, manipulator.paths)
            manipulator.rdy_update_pdi.clear()
            lib_utils.update_pdi(manipulator, tmp_image)
    
            t = time.time()-t
            logging.info(lib_utils.log_message(f"Pathfinder finished after {t:.5f} seconds"))
            
            # Move probe if in "Auto Manipulation" operation mode.
            if auto_manipulate:
                move_probe(manipulator)
            
            # Trigger ready-event.
            manipulator.pathfinding_module.rdy.set()

        if auto_manipulate:
            logging.info(lib_utils.log_message("Pathfinder stopped."))
        
    manipulator.t5 = threading.Thread(target = do_this)
    manipulator.t5.start()


# Wrapper for adding/removing foreign atoms / target sites
def add_or_remove_foreign_atoms_or_target_sites(manipulator, mode=None, startstop=False):
    # Mode description:
    # 0 ... Add foreign atoms
    # 1 ... Remove foreign atoms
    # 2 ... Add target sites
    # 3 ... Remove target sites

    # Aliases.
    dc = manipulator.document_controller

    if not hasattr(manipulator, 'original_mouse_clicked'):
        manipulator.original_mouse_clicked = dc._document_controller.selected_display_panel.canvas_widget.on_mouse_clicked

    if startstop: # Check boolean (start/stop)
        def on_mouse_clicked_custom(x, y, modifiers):
            manipulator.original_mouse_clicked(x, y, modifiers)
            if dc._document_controller.selected_display_panel.data_item == manipulator.processed_data_item._data_item:
                canvas_item = dc._document_controller.selected_display_panel.\
                    root_container._RootCanvasItem__mouse_tracking_canvas_item
                if canvas_item:
                    canvas_item_point = dc._document_controller.selected_display_panel.\
                        root_container.map_to_canvas_item(Geometry.IntPoint(y=y, x=x), canvas_item)
                    if hasattr(canvas_item, 'map_widget_to_image'):
                        image_point = canvas_item.map_widget_to_image(canvas_item_point)
                        # image_point contains the mouse click position on the image in px relative to top-left corner.
                        if mode == 0 or mode == 2:
                            add_atom_or_site_near_image_point(manipulator, image_point, mode)
                        elif mode == 1 or mode == 3:
                            remove_atom_or_site_near_image_point(manipulator, image_point, mode)
                        lib_utils.refresh_GUI(manipulator, ['foreigns', 'targets'])
        
        dc._document_controller.selected_display_panel.canvas_widget.on_mouse_clicked = on_mouse_clicked_custom
    else:
        dc._document_controller.selected_display_panel.canvas_widget.on_mouse_clicked = manipulator.original_mouse_clicked


# Add foreign atom or target site by grabbing the position of the mouse cursor.
def add_atom_or_site_near_image_point(manipulator, image_point, mode):
    
    if manipulator.sites:
        # Find nearest atom site.
        all_coords = map(lambda x: x.coords, manipulator.sites)
        nearest_site = manipulator.sites[ np.linalg.norm(
            np.array(list(all_coords))-image_point, axis=1).argmin() ]
        
        # Style.
        relative_size = 0.05
        shape = manipulator.source_xdata.data_shape

        if mode == 0: # Add foreign atom.
            manipulator.sources.append( aab.Atom(nearest_site, 'pseudo-element', defined_by_user=True))
            # Insert region.
            manipulator.rectangle_regions.append( manipulator.processed_data_item.add_rectangle_region(
                nearest_site.coords[0]/shape[0], nearest_site.coords[1]/shape[1], relative_size, relative_size))
            # Mutual variable assignment.
            manipulator.rectangle_regions[-1].atom = manipulator.sources[-1]
            manipulator.sources[-1].graphic = manipulator.rectangle_regions[-1]
            lib_utils.add_listener_graphic_changed(manipulator, manipulator.sources[-1].graphic)
        
        elif mode == 2: # Add target site.
            manipulator.targets.append( nearest_site )
            # Insert region.
            manipulator.ellipse_regions.append( manipulator.processed_data_item.add_ellipse_region(
                nearest_site.coords[0]/shape[0], nearest_site.coords[1]/shape[1], relative_size, relative_size))
            # Mutual variable assignment.
            manipulator.ellipse_regions[-1].site = nearest_site
            manipulator.targets[-1].graphic = manipulator.ellipse_regions[-1]
            lib_utils.add_listener_graphic_changed(manipulator, manipulator.targets[-1].graphic)
    else:
        logging.info(lib_utils.log_message("No atom sites found."))


# Remove foreign atom or target site by grabbing the position of the mouse cursor.
def remove_atom_or_site_near_image_point(manipulator, image_point, mode):
    
    if mode == 1: # Remove foreign atom.
        if manipulator.sources:
            # Find nearest foreign atom.
            all_coords = map(lambda x: x.site.coords, manipulator.sources)
            nearest_source = manipulator.sources[ np.linalg.norm(
                np.array(list(all_coords))-image_point, axis=1).argmin() ]
            # Remove object and region.
            manipulator.sources.remove(nearest_source)
            try:
                manipulator.rectangle_regions.remove(nearest_source.graphic)
            except:
                manipulator.rectangle_regions_auto.remove(nearest_source.graphic)    
            manipulator.processed_data_item.remove_region(nearest_source.graphic)
        else:
            logging.info(lib_utils.log_message("No foreign atoms found."))
    
    if mode == 3: # Remove target site
        if manipulator.targets:
            # Find nearest target site.
            all_coords = map(lambda x: x.coords, manipulator.targets)
            nearest_target = manipulator.targets[ np.linalg.norm(
                np.array(list(all_coords))-image_point, axis=1).argmin() ]
            # Remove object and region.
            manipulator.targets.remove(nearest_target)
            manipulator.ellipse_regions.remove(nearest_target.graphic)
            manipulator.processed_data_item.remove_region(nearest_target.graphic)
        else:
            logging.info(lib_utils.log_message("No target sites found."))

  
# Set probe position.
def move_probe(manipulator):

    if not hasattr(manipulator.paths, 'members'):
        logging.info(lib_utils.log_message("No paths found. Probe not repositioned."))
    
    # Choose next path that is longer than 1 site
    # (1 site >=> atom is already at target site).
    yx = None
    for path_number, path in enumerate(manipulator.paths.members):
        if len(path.sitelist) >= 2:
            yx = path.sitelist[1].coords
            break
        else:
            continue
    if manipulator.superscan._hardware_source.probe_position is not None:
        if yx is not None:
            yx_frame = manipulator.superscan.get_frame_parameters()["size"]
            yx_frac = yx / yx_frame
            manipulator.superscan._hardware_source.probe_position = list(yx_frac)
            logging.info(lib_utils.log_message("Probe repositioned."))
        else:
            logging.info(lib_utils.log_message("No paths found. Probe not repositioned."))
    else:
        pass