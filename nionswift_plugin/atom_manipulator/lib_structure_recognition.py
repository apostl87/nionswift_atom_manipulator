"""
Structure recognition algorithm and library.
- Performed by a deep convolutional neural network (NN).
    1) Calibration of the image scale with FourierScaleCalibrator or RealSpaceCalibrator.
    2) Identification and separation of pristine lattice and contaminated areas.
    3) Recognition of atom positions in a 2D material.
    4) Detection of substitutional heteroatoms/dopants.
- Elemental identification of dopants: integration over specified radius and evaluation with a specified Z-radius.
"""

import gettext
import threading
import numpy as np

import copy
import time
import logging

import math

# Non-standard packages
from fourier_scale_calibration.fourier_scale_calibration import FourierSpaceCalibrator

# Custom libraries
from .classes import atoms_and_bonds as aab
from . import lib_utils
from . import lib_pathfinding 

_ = gettext.gettext
   

# Main structure recognition function.
def analyze_and_show(structure_recognition_module, auto_manipulate=False, live_analysis=False):
    if structure_recognition_module.manipulator.t1 is not None and structure_recognition_module.manipulator.t1.is_alive():
            logging.info(lib_utils.log_message("Structure recognition still working. Wait until finished"))
            return
        
    # Aliases.
    manipulator = structure_recognition_module.manipulator

    # Assign parameters.
    imgsrc = structure_recognition_module.image_source_combo_box.current_item[0:5].upper()
    if auto_manipulate:
        structure_recognition_module.stop_live_analysis_event.clear()
        live_analysis = False
        if imgsrc == "SELEC": # Fix variable imgsrc if user has not changed it to MAADF or HAADF.
            imgsrc = "FIRST"
        
    # Create processed data item if none exists.
    if manipulator.processed_data_item is None or manipulator.processed_data_item not in manipulator.api.library.data_items:
        manipulator.rdy_create_pdi.clear()
        lib_utils.create_pdi(manipulator)

    def do_this():
        while not (not auto_manipulate and structure_recognition_module.stop_live_analysis_event.is_set()) \
              and \
              not (auto_manipulate and manipulator.manipulation_module.stop_auto_manipulate_event.is_set()):
                  
                if not live_analysis:
                    structure_recognition_module.stop_live_analysis_event.set()

                wait_time = 2
                if manipulator.rdy_create_pdi.wait(wait_time): # wait for a maximum of {wait_time} seconds.
                    pass
                else:
                    logging.info(lib_utils.log_message("Waiting for UI thread for more than {wait_time} seconds."
                                                       "Possible code performance issue or a crashed thread."))
             
                if auto_manipulate:
                    if not manipulator.tractor_beam_module.rdy.wait(0.1):
                        #print("waiting for tb module")
                        continue
                manipulator.tractor_beam_module.rdy.clear()
                
                if "SELEC" in imgsrc:
                    structure_recognition_module.stop_live_analysis_event.set()
                    tdi = manipulator.document_controller.target_data_item
                    manipulator.source_xdata = tdi.xdata
                    manipulator.source_title = tdi.title
                    manipulator.scan_parameters_changed = True
                else:
                    logging.info(lib_utils.log_message("Grabbing next STEM image ..."))
                    if (not manipulator.superscan.is_playing and not live_analysis) or auto_manipulate:
                        # Start and stop scanning when these conditions are met.
                        manipulator.superscan.start_playing()
                        time.sleep(0.05) # Small delay is necessary due to a delayed response of the Nion Swift code.
                        manipulator.superscan.stop_playing()
                        
                    last_record = manipulator.superscan.grab_next_to_finish()
                    for item in last_record:
                        if imgsrc == "FIRST" or item.metadata['hardware_source']['channel_name'] == imgsrc:
                            scan_parameters = manipulator.superscan.get_frame_parameters()
                            scan_parameters = np.array((scan_parameters['fov_nm'], 'placeholder'))
                            if manipulator.scan_parameters is None:
                                manipulator.scan_parameters = scan_parameters
                            if all(scan_parameters==manipulator.scan_parameters):
                                manipulator.scan_parameters_changed = False
                            else:
                                manipulator.scan_parameters_changed = True
                                manipulator.scan_parameters = scan_parameters
                            manipulator.source_xdata = item
                            manipulator.source_title = \
                                item.metadata['hardware_source']['hardware_source_name'] + \
                                ' (' + imgsrc + ")"
                            break
                        
                if auto_manipulate:
                    structure_recognition_module.new_image.set() # For TractorBeam module.
                    
                # Aliases.
                pdi = manipulator.processed_data_item
                shape = np.array(manipulator.source_xdata.data_shape)

                if manipulator.source_xdata is None:
                    logging.info(lib_utils.log_message("No source data."))
                    structure_recognition_module.stop_live_analysis_event.set()
                    break
                
                # Calibrates the image scale based on a Fourier transform of the lattice.
                if structure_recognition_module.scale_calibration_mode == 1:
                    t = time.time()
                    logging.info(lib_utils.log_message("FourierSpaceCalibrator called."))
                
                    calibrator = FourierSpaceCalibrator('hexagonal', 2.46)
                    structure_recognition_module.sampling = calibrator(manipulator.source_xdata.data)
                    structure_recognition_module.fov = [structure_recognition_module.sampling*s for s in manipulator.source_xdata.data.shape]
                
                    t = time.time()-t
                    logging.info(lib_utils.log_message(f"FourierSpaceCalibrator finished after {t:.5f} seconds."))
                else:
                    # structure_recognition_module.sampling must have been written before.
                    if structure_recognition_module.sampling is None:
                        logging.info(lib_utils.log_message("Saved value for 'sampling' is None. Stopping..."))
                        return None # Stop manipulator.

                # Call deep convolutional neural network.
                t = time.time()
                logging.info(lib_utils.log_message("Structure recognition called."))
                
                structure_recognition_module.nn_output = structure_recognition_module.model(manipulator.source_xdata.data, structure_recognition_module.sampling)
                
                t = time.time()-t
                logging.info(lib_utils.log_message(f"Structure recognition finished after {t:.5f} seconds."))

                manipulator.rdy_init_pdi.clear()
                lib_utils.init_pdi(manipulator) # This is done here to give the user a possibility to look at the paths.

                # Conditioning NN output
                if structure_recognition_module.nn_output is not None:
                    manipulator.maxima_locations = np.fliplr(structure_recognition_module.nn_output['points'])
                    number_maxima = len(manipulator.maxima_locations)
                else:
                    manipulator.maxima_locations = None
                    number_maxima = 0
                    
                logging.info(lib_utils.log_message(f"{number_maxima:d} atoms were found."))
                
                # Call object-oriented backend to draw atom positions and bonds.
                t = time.time()
                manipulator.sites = []
                manipulator.paths = []
                for i in range(number_maxima):
                    loc = manipulator.maxima_locations[i]
                    manipulator.sites.append( aab.Site(
                        loc[0], loc[1], site_id=i) )
                        
                # Refresh GUI.
                lib_utils.refresh_GUI(manipulator, ['atoms', 'sampling'])
                t = time.time()-t
                logging.info(lib_utils.log_message(f"Setting sites (back end) finished after {t:.5f} seconds."))
        
                if structure_recognition_module.visualize_atoms:
                    while not manipulator.rdy_init_pdi.wait(1):
                        pass
                    tmp_image = copy.copy(pdi.data)
                    lib_utils.plot_points(tmp_image, manipulator.maxima_locations)
                    manipulator.rdy_update_pdi.clear()
                    lib_utils.update_pdi(manipulator, tmp_image)
       
                # Try to keep target sites and foreign atoms till the next frame.
                if manipulator.scan_parameters_changed: # re-init
                    clear_user_defined_atoms_and_targets(manipulator)
                
                else: # Reposition foreign atoms, target sites and the corresponding graphics.
                    
                    all_coords = np.full((len(manipulator.sites), 2), np.nan)
                    if len(manipulator.sources) > 0 or len(manipulator.targets) > 0:
                        for i, site in enumerate(manipulator.sites):
                            all_coords[i, :] = site.coords
                    else:
                        pass
                       
                    # Target sites.
                    graphics_pos = np.full((len(manipulator.targets), 2), np.nan)
                    for i, target in enumerate(manipulator.targets):
                        graphics_pos[i, :] = target.graphic.center
                    graphics_pos *= shape
                        
                    distances = np.linalg.norm(graphics_pos - all_coords[:, np.newaxis, :], axis=2)
                    site_indices = np.argmin(distances, axis=0)
                    
                    for i, site_index in enumerate(site_indices):
                        graphic = manipulator.targets[i].graphic
                        manipulator.targets[i] = manipulator.sites[site_index]
                        manipulator.targets[i].graphic = graphic
                    def reposition_target_site_graphics():
                        with manipulator.api.library.data_ref_for_data_item(manipulator.processed_data_item):
                            for target in manipulator.targets:
                                target.graphic.center = target.coords / shape
                    manipulator.api.queue_task(reposition_target_site_graphics)
                     
                    # Foreign atoms (only user-defined).
                    atoms_user_def = []
                    for atom in manipulator.sources:
                        if atom.defined_by_user: atoms_user_def.append(atom)
                        
                    graphics_pos = np.full((len(atoms_user_def), 2), np.nan)
                    for i, atom in enumerate(atoms_user_def):
                        graphics_pos[i, :] = atom.graphic.center
                    graphics_pos *= shape
                        
                    distances = np.linalg.norm(graphics_pos - all_coords[:, np.newaxis, :], axis=2)
                    site_indices = np.argmin(distances, axis=0)
                    
                    for i, site_index in enumerate(site_indices):
                        atoms_user_def[i].site = manipulator.sites[site_index]
                        atoms_user_def[i].origin = manipulator.sources[i].site

                    def reposition_foreign_atom_graphics():
                        with manipulator.api.library.data_ref_for_data_item(manipulator.processed_data_item):
                            for atom in atoms_user_def:
                                atom.graphic.center = atom.origin.coords / shape
                    manipulator.api.queue_task(reposition_foreign_atom_graphics)
                    
                    t = time.time()-t
                    #devel anchor
                    logging.info(lib_utils.log_message(f"Repositioning of foreign atoms, target sites, "
                                                       f"and graphics finished after {t:.5f} seconds."))
                
                # Auto-detect of sources.
                func_auto_detect_foreign_atoms(structure_recognition_module)
                structure_recognition_module.rdy.set() # Signal that this is ready.
                #print("why error when deleting a data item")
                
                # Integrate pathfinding in live analysis.
                if live_analysis and (len(manipulator.sources)>=0 and len(manipulator.targets)>=0):
                    lib_pathfinding.find_paths(manipulator)
                
        if auto_manipulate:
            logging.info(lib_utils.log_message("Structure recognition stopped."))
            
    manipulator.t1 = threading.Thread(target = do_this, name = 'do_this')
    manipulator.t1.start()
    
def create_site_and_connect_to_point_region(manipulator, loc, site_id, point_region):
    manipulator.sites.append( aab.Site(loc[0], loc[1], site_id=site_id) )
    # Mutual assignment
    manipulator.sites[-1].graphic = point_region
    point_region.site = manipulator.sites[-1]
            
# Auto-detect and display foreign atoms.
def func_auto_detect_foreign_atoms(structure_recognition_module):
    # aliases
    manipulator = structure_recognition_module.manipulator
    rra = manipulator.rectangle_regions_auto # "rra" means rectangle regions auto-detected
    pdi = manipulator.processed_data_item
    shape = manipulator.source_xdata.data_shape
    
    # Delete old auto-detected sources.
    atoms = copy.copy(manipulator.sources)
    for atom in atoms:
        if not atom.defined_by_user:
            manipulator.sources.remove(atom)
    
    foreigns_site_id = []
    if structure_recognition_module.auto_detect_foreign_atoms:
        foreigns_site_id = np.nonzero(structure_recognition_module.nn_output['labels'] == 1)[0]
        logging.info(lib_utils.log_message(f"Detected {len(foreigns_site_id):d} foreign atoms."))
      
    relative_size = 0.05
    number_foreigns = len(foreigns_site_id)
    
    while not manipulator.rdy_update_pdi.wait(1):
        logging.info(lib_utils.log_message("Waiting for an update of the processed data item..."))
        pass
    
    new_centers = []    
    for i in range(number_foreigns):
        site_id = foreigns_site_id[i]
        loc = manipulator.maxima_locations[site_id]
        added = False
        
        #debug anchor
        #print(loc[0]/shape[0])
        #print(loc[1]/shape[1])
        
        if i < len(rra):
            # Reposition if a graphic exists
            new_centers.append( (max(0, loc[0]/shape[0]), min(1, loc[1]/shape[1])) )
        else:
            # Insert region
            rra.append(pdi.add_rectangle_region( max(0, loc[0]/shape[0]), min(1, loc[1]/shape[1]), relative_size, relative_size ))
            added = True
                
        manipulator.sources.append( aab.Atom(manipulator.sites[site_id], 'pseudo-element') )
        # Mutual assignment
        rra[i].atom = manipulator.sources[-1]
        manipulator.sources[-1].graphic = rra[i]
        if added:
            lib_utils.add_listener_graphic_changed(manipulator, rra[i])
            
    if len(new_centers) > 0:
        for i, center in enumerate(new_centers):
            rra[i].center = center
        
    if len(rra) > number_foreigns: # Remove graphics
        #print(" ! Removing excess RECTANGLE graphics ")
        manipulator.rectangle_regions_auto = rra[0:number_foreigns]
        def func():
            for k in range(number_foreigns, len(rra)):
                pdi.remove_region(rra[k])
        manipulator.api.queue_task(func)
    
    def do_this(): # The following is threaded, because it is not needed for later processes.
        lib_utils.refresh_GUI(manipulator, ['foreigns'])
        lib_utils.elemental_identification(manipulator)
    threading.Thread(target=do_this).start()

def clear_user_defined_atoms_and_targets(manipulator):
    pdi = manipulator.processed_data_item
    def func():
        with manipulator.api.library.data_ref_for_data_item(pdi):
            for region in manipulator.rectangle_regions:
                pdi.remove_region(region)
            for region in manipulator.ellipse_regions:
                pdi.remove_region(region)
        manipulator.rectangle_regions = []
        manipulator.ellipse_regions = []
    manipulator.api.queue_task(func)
    manipulator.sources = []
    manipulator.targets = []

    logging.info(lib_utils.log_message("Cleared all user-defined foreign atoms and target sites."))
    lib_utils.refresh_GUI(manipulator, ['foreigns', 'targets'])