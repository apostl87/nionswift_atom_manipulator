"""
Path finding
- Determines shortest total path from sources to targets 
"""

import copy
import time
import threading
import logging
import numpy as np
from .classes import atoms_and_bonds as aab, paths
from nion.utils import Geometry
from . import lib_utils

# Path finding
def find_paths(manip_obj, auto_manipulate=False):
    if (manip_obj.sites == []) and not auto_manipulate: # or (manip_obj.targets is None):
            logging.info(lib_utils.log_message("No sites found. Path finder aborted."))
            return
    if manip_obj.t5 is not None and manip_obj.t5.is_alive():
            logging.info(lib_utils.log_message("Path finder still working. Wait until finished."))
            return 
    
    def thread_this():
        stop = False
        while not stop and not (auto_manipulate and manip_obj.manipulation_module.stop_auto_manipulate_event.is_set()):
            if auto_manipulate:
                if not manip_obj.structure_recognition_module.rdy.wait(0.1):
                    continue
                else:
                    manip_obj.structure_recognition_module.rdy.clear()
            else:
                stop = True

            # Aliases
            pdi = manip_obj.processed_data_item
            
            # Set bonds
            t = time.time()
            logging.info(lib_utils.log_message("Setting bonds ..."))
            max_bond_length_px = manip_obj.source_xdata.dimensional_calibrations[0].convert_from_calibrated_size(
                    manip_obj.path_finding_module.max_bond_length/10)
            if manip_obj.simulation_mode:
                max_bond_length_px *= 10 # Wrong conversion im usim fork
            manip_obj.bonds = aab.Bonds(manip_obj.sites, max_bond_length_px)
            
            t = time.time()-t
            logging.info(lib_utils.log_message(f"Setting bonds finished after {t:.5f} seconds"))
            
            # Call path finder
            t = time.time()
            logging.info(lib_utils.log_message("Path finder called."))
            
            try:
                manip_obj.paths = paths.Paths(manip_obj.sources, manip_obj.targets)
            except ValueError as e:
                print(e)
                return
            else:
                manip_obj.paths.determine_paths_nooverlap()
            
            # Display paths
            while not manip_obj.rdy_init_pdi.wait(1) or not manip_obj.rdy_update_pdi.wait(1):
                pass
            tmp_image = copy.copy(pdi.data)
            lib_utils.plot_paths(tmp_image, manip_obj.paths)
            manip_obj.rdy_update_pdi.clear()
            lib_utils.update_pdi(manip_obj, tmp_image)
    
            t = time.time()-t
            logging.info(lib_utils.log_message(f"Path finder finished after {t:.5f} seconds"))
            
            if auto_manipulate:
                move_probe(manip_obj)
            manip_obj.path_finding_module.rdy.set()

        if auto_manipulate:
            logging.info(lib_utils.log_message("Path finder stopped."))
        
    manip_obj.t5 = threading.Thread(target = thread_this)
    manip_obj.t5.start()

# Add foreign atoms or targets with grabbing the mouse cursor
def add_or_remove_foreign_atoms_or_target_sites(manip_obj, mode=None, startstop=False):
    # Modes:
    # 0 ... Add foreign atoms
    # 1 ... Remove foreign atoms
    # 2 ... Add target sites
    # 3 ... Remove target sites
    dc = manip_obj.document_controller
    if not hasattr(manip_obj, 'original_mouse_clicked'):
        manip_obj.original_mouse_clicked = dc._document_controller.selected_display_panel.canvas_widget.on_mouse_clicked

    if startstop: # Check boolean (start/stop)
        def on_mouse_clicked_custom(x, y, modifiers):
            manip_obj.original_mouse_clicked(x, y, modifiers)
            if dc._document_controller.selected_display_panel.data_item == manip_obj.processed_data_item._data_item:
                canvas_item = dc._document_controller.selected_display_panel.root_container._RootCanvasItem__mouse_tracking_canvas_item
                if canvas_item:
                    canvas_item_point = dc._document_controller.selected_display_panel.root_container.map_to_canvas_item(Geometry.IntPoint(y=y, x=x), canvas_item)
                    if hasattr(canvas_item, 'map_widget_to_image'):
                        image_point = canvas_item.map_widget_to_image(canvas_item_point)
                        # image_point contains the mouse click position in the image in pixels relative to top-left corner
                        if mode == 0 or mode == 2:
                            add_atom_or_site_near_image_point(manip_obj, image_point, mode)
                        elif mode == 1 or mode == 3:
                            remove_atom_or_site_near_image_point(manip_obj, image_point, mode)
                        lib_utils.refresh_GUI(manip_obj, ['foreigns', 'targets'])
                        #print(image_point)
        
        dc._document_controller.selected_display_panel.canvas_widget.on_mouse_clicked = on_mouse_clicked_custom
    else:
        dc._document_controller.selected_display_panel.canvas_widget.on_mouse_clicked = manip_obj.original_mouse_clicked

def add_atom_or_site_near_image_point(manip_obj, image_point, mode):
    # Find nearest atom site
    if manip_obj.sites:
        all_coords = map(lambda x: x.coords, manip_obj.sites)
        nearest_site = manip_obj.sites[ np.linalg.norm(
            np.array(list(all_coords))-image_point, axis=1).argmin() ]
        relative_size = 0.05
        shape = manip_obj.source_xdata.data_shape
        if mode == 0:
            # Add foreign atom
            manip_obj.sources.append( aab.Atom(nearest_site, 'pseudo-element', defined_by_user=True))
            # Insert region
            manip_obj.rectangle_regions.append( manip_obj.processed_data_item.add_rectangle_region(
                nearest_site.coords[0]/shape[0], nearest_site.coords[1]/shape[1], relative_size, relative_size))
            # Mutual assignment
            manip_obj.rectangle_regions[-1].atom = manip_obj.sources[-1]
            manip_obj.sources[-1].graphic = manip_obj.rectangle_regions[-1]
            lib_utils.add_listener_graphic_changed(manip_obj, manip_obj.sources[-1].graphic)
        elif mode == 2:
            # Add target site
            manip_obj.targets.append( nearest_site )
            # Insert region
            manip_obj.ellipse_regions.append( manip_obj.processed_data_item.add_ellipse_region(
                nearest_site.coords[0]/shape[0], nearest_site.coords[1]/shape[1], relative_size, relative_size))
            # Mutual assignment
            manip_obj.ellipse_regions[-1].site = nearest_site
            manip_obj.targets[-1].graphic = manip_obj.ellipse_regions[-1]
            lib_utils.add_listener_graphic_changed(manip_obj, manip_obj.targets[-1].graphic)
    else:
        logging.info(lib_utils.log_message("No atom sites found."))
        
def remove_atom_or_site_near_image_point(manip_obj, image_point, mode):
    if mode == 1:
        if manip_obj.sources:
            # Find nearest foreign atom
            all_coords = map(lambda x: x.site.coords, manip_obj.sources)
            nearest_source = manip_obj.sources[ np.linalg.norm(
                np.array(list(all_coords))-image_point, axis=1).argmin() ]
            # Remove foreign atom
            manip_obj.sources.remove(nearest_source)
            # Remove region
            try:
                manip_obj.rectangle_regions.remove(nearest_source.graphic)
            except:
                manip_obj.rectangle_regions_auto.remove(nearest_source.graphic)    
            manip_obj.processed_data_item.remove_region(nearest_source.graphic)
        else:
            logging.info(lib_utils.log_message("No foreign atoms found."))
    if mode == 3:
        if manip_obj.targets:
            # Find nearest target site
            all_coords = map(lambda x: x.coords, manip_obj.targets)
            nearest_target = manip_obj.targets[ np.linalg.norm(
                np.array(list(all_coords))-image_point, axis=1).argmin() ]
            # Remove foreign atom
            manip_obj.targets.remove(nearest_target)
            # Remove region
            manip_obj.ellipse_regions.remove(nearest_target.graphic)
            manip_obj.processed_data_item.remove_region(nearest_target.graphic)
        else:
            logging.info(lib_utils.log_message("No target sites found."))
    
# Set probe position
def move_probe(manip_obj):
    # Choose next path that is longer than 1 site (1 site >=> atom is already at target site)
    if not hasattr(manip_obj.paths, 'members'):
        logging.info(lib_utils.log_message("No paths found. Probe not repositioned."))
    yx = None
    for path_number, path in enumerate(manip_obj.paths.members):
        if len(path.sitelist) >= 2:
            yx = path.sitelist[1].coords
            break
        else:
            continue
    if manip_obj.superscan._hardware_source.probe_position is not None:
        if yx is not None:
            yx_frame = manip_obj.superscan.get_frame_parameters()["size"]
            yx_frac = yx / yx_frame
            manip_obj.superscan._hardware_source.probe_position = list(yx_frac)
            logging.info(lib_utils.log_message("Probe repositioned."))
        else:
            logging.info(lib_utils.log_message("No paths found. Probe not repositioned."))
    else:
        pass