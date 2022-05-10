"""
Structure recognition is performed by a deep convolutional neural network (NN).
- Identification and separation of pristine lattice and contaminated areas (NN).
- Recognition of atom positions in a 2D material (NN).
- Detection of substitutional heteroatoms (ie. dopants) (NN).
- Elemental identification: integration over specific radius (operational, but almost never works correctly #TODO).
"""

import time
import gettext
import copy
import threading
import logging
import math
import numpy as np
from .classes import atoms_and_bonds as aab
from fourier_scale_calibration.fourier_scale_calibration import FourierSpaceCalibrator

from . import lib_utils
from . import lib_path_finding 

_ = gettext.gettext
   
# Analysis main function.
def analyze_and_show(sr_obj, auto_manipulate=False, live_analysis=False):
    if sr_obj.manip_obj.t1 is not None and sr_obj.manip_obj.t1.is_alive():
            logging.info(lib_utils("Structure recognition still working. Wait until finished"))
            return
        
    # Aliases.
    manip_obj = sr_obj.manip_obj

    # Assign parameters.
    imgsrc = sr_obj.image_source_combo_box.current_item[0:5].upper()
    if auto_manipulate:
        sr_obj.stop_live_analysis_event.clear()
        live_analysis = False
        if imgsrc == "SELEC":
            imgsrc = "FIRST" #Correct that if user has not changed "Image source"
        
    # Create processed data item if none exists.
    if manip_obj.processed_data_item is None or manip_obj.processed_data_item not in manip_obj.api.library.data_items:
        manip_obj.rdy_create_pdi.clear()
        lib_utils.create_pdi(manip_obj)

    def do_this():
        while not (not auto_manipulate and sr_obj.stop_live_analysis_event.is_set()) \
              and \
              not (auto_manipulate and manip_obj.manipulation_module.stop_auto_manipulate_event.is_set()):
                  
                if not live_analysis:
                    sr_obj.stop_live_analysis_event.set()

                wait_time = 2
                if manip_obj.rdy_create_pdi.wait(wait_time): # wait for a maximum of {wait_time} seconds.
                    pass
                else:
                    logging.info(lib_utils.log_message("Waiting for UI thread for more than {wait_time} seconds."
                                                       "Possible code performance issue or a crashed thread."))
             
                if auto_manipulate:
                    if not manip_obj.tractor_beam_module.rdy.wait(0.1):
                        #print("waiting for tb module")
                        continue
                manip_obj.tractor_beam_module.rdy.clear()
                
                if "SELEC" in imgsrc:
                    sr_obj.stop_live_analysis_event.set()
                    tdi = manip_obj.document_controller.target_data_item
                    manip_obj.source_xdata = tdi.xdata
                    manip_obj.source_title = tdi.title
                    manip_obj.scan_parameters_changed = True
                else:
                    logging.info(lib_utils.log_message("Grabbing next STEM image ..."))
                    if (not manip_obj.superscan.is_playing and not live_analysis) or auto_manipulate:
                        # Start and stop scanning when these conditions are met.
                        manip_obj.superscan.start_playing()
                        time.sleep(0.05) # Small delay is necessary due to a delayed response of the Nion Swift code.
                        manip_obj.superscan.stop_playing()
                        
                    last_record = manip_obj.superscan.grab_next_to_finish()
                    for item in last_record:
                        if imgsrc == "FIRST" or item.metadata['hardware_source']['channel_name'] == imgsrc:
                            scan_parameters = manip_obj.superscan.get_frame_parameters()
                            scan_parameters = np.array((scan_parameters['fov_nm'], 'placeholder'))
                            if manip_obj.scan_parameters is None:
                                manip_obj.scan_parameters = scan_parameters
                            if all(scan_parameters==manip_obj.scan_parameters):
                                manip_obj.scan_parameters_changed = False
                            else:
                                manip_obj.scan_parameters_changed = True
                                manip_obj.scan_parameters = scan_parameters
                            manip_obj.source_xdata = item
                            manip_obj.source_title = \
                                item.metadata['hardware_source']['hardware_source_name'] + \
                                ' (' + imgsrc + ")"
                            break
                        
                if auto_manipulate:
                    sr_obj.new_image.set() # For TractorBeam module.
                    
                # Aliases.
                pdi = manip_obj.processed_data_item
                shape = np.array(manip_obj.source_xdata.data_shape)

                if manip_obj.source_xdata is None:
                    logging.info(lib_utils.log_message("No source data."))
                    sr_obj.stop_live_analysis_event.set()
                    break
                
                # Calibrates the image scale based on a Fourier transform of the lattice.
                if sr_obj.scale_calibration_mode == 1:
                    t = time.time()
                    logging.info(lib_utils.log_message("FourierSpaceCalibrator called."))
                
                    calibrator = FourierSpaceCalibrator('hexagonal', 2.46)
                    sr_obj.sampling = calibrator(manip_obj.source_xdata.data)
                    sr_obj.fov = [sr_obj.sampling*s for s in manip_obj.source_xdata.data.shape]
                
                    t = time.time()-t
                    logging.info(lib_utils.log_message(f"FourierSpaceCalibrator finished after {t:.5f} seconds."))
                else:
                    # sr_obj.sampling must have been written before.
                    if sr_obj.sampling is None:
                        logging.info(lib_utils.log_message("Saved value for 'sampling' is None. Stopping..."))
                        return None # Stop manipulator.

                # Call deep convolutional neural network.
                t = time.time()
                logging.info(lib_utils.log_message("Structure recognition called."))
                
                sr_obj.nn_output = sr_obj.model(manip_obj.source_xdata.data, sr_obj.sampling)
                
                t = time.time()-t
                logging.info(lib_utils.log_message(f"Structure recognition finished after {t:.5f} seconds."))

                manip_obj.rdy_init_pdi.clear()
                lib_utils.init_pdi(manip_obj) # This is done here to give the user a possibility to look at the paths.

                # Conditioning NN output
                if sr_obj.nn_output is not None:
                    manip_obj.maxima_locations = np.fliplr(sr_obj.nn_output['points'])
                    number_maxima = len(manip_obj.maxima_locations)
                else:
                    manip_obj.maxima_locations = None
                    number_maxima = 0
                    
                logging.info(lib_utils.log_message(f"{number_maxima:d} atoms were found."))
                
                # Call object-oriented backend to draw atom positions and bonds.
                t = time.time()
                manip_obj.sites = []
                manip_obj.paths = []
                for i in range(number_maxima):
                    loc = manip_obj.maxima_locations[i]
                    manip_obj.sites.append( aab.Site(
                        loc[0], loc[1], site_id=i) )
                        
                # Refresh GUI.
                lib_utils.refresh_GUI(manip_obj, ['atoms', 'sampling'])
                t = time.time()-t
                logging.info(lib_utils.log_message(f"Setting sites (back end) finished after {t:.5f} seconds."))
        
                if sr_obj.visualize_atoms:
                    while not manip_obj.rdy_init_pdi.wait(1):
                        pass
                    tmp_image = copy.copy(pdi.data)
                    lib_utils.plot_points(tmp_image, manip_obj.maxima_locations)
                    manip_obj.rdy_update_pdi.clear()
                    lib_utils.update_pdi(manip_obj, tmp_image)
       
                # Try to keep target sites and foreign atoms till the next frame.
                if manip_obj.scan_parameters_changed: # re-init
                    clear_user_defined_atoms_and_targets(manip_obj)
                
                else: # Reposition foreign atoms, target sites and the corresponding graphics.
                    
                    all_coords = np.full((len(manip_obj.sites), 2), np.nan)
                    if len(manip_obj.sources) > 0 or len(manip_obj.targets) > 0:
                        for i, site in enumerate(manip_obj.sites):
                            all_coords[i, :] = site.coords
                    else:
                        pass
                       
                    # Target sites.
                    graphics_pos = np.full((len(manip_obj.targets), 2), np.nan)
                    for i, target in enumerate(manip_obj.targets):
                        graphics_pos[i, :] = target.graphic.center
                    graphics_pos *= shape
                        
                    distances = np.linalg.norm(graphics_pos - all_coords[:, np.newaxis, :], axis=2)
                    site_indices = np.argmin(distances, axis=0)
                    
                    for i, site_index in enumerate(site_indices):
                        graphic = manip_obj.targets[i].graphic
                        manip_obj.targets[i] = manip_obj.sites[site_index]
                        manip_obj.targets[i].graphic = graphic
                    def reposition_target_site_graphics():
                        with manip_obj.api.library.data_ref_for_data_item(manip_obj.processed_data_item):
                            for target in manip_obj.targets:
                                target.graphic.center = target.coords / shape
                    manip_obj.api.queue_task(reposition_target_site_graphics)
                     
                    # Foreign atoms (only user-defined).
                    atoms_user_def = []
                    for atom in manip_obj.sources:
                        if atom.defined_by_user: atoms_user_def.append(atom)
                        
                    graphics_pos = np.full((len(atoms_user_def), 2), np.nan)
                    for i, atom in enumerate(atoms_user_def):
                        graphics_pos[i, :] = atom.graphic.center
                    graphics_pos *= shape
                        
                    distances = np.linalg.norm(graphics_pos - all_coords[:, np.newaxis, :], axis=2)
                    site_indices = np.argmin(distances, axis=0)
                    
                    for i, site_index in enumerate(site_indices):
                        atoms_user_def[i].site = manip_obj.sites[site_index]
                        atoms_user_def[i].origin = manip_obj.sources[i].site

                    def reposition_foreign_atom_graphics():
                        with manip_obj.api.library.data_ref_for_data_item(manip_obj.processed_data_item):
                            for atom in atoms_user_def:
                                atom.graphic.center = atom.origin.coords / shape
                    manip_obj.api.queue_task(reposition_foreign_atom_graphics)
                    
                    t = time.time()-t
                    #devel anchor
                    logging.info(lib_utils.log_message(f"Repositioning of foreign atoms, target sites, "
                                                       f"and graphics finished after {t:.5f} seconds."))
                
                # Auto-detect of sources.
                func_auto_detect_foreign_atoms(sr_obj)
                sr_obj.rdy.set() # Signal that this is ready.
                #print("why error when deleting a data item")
                
                # Integrate path finder in live analysis.
                if live_analysis and (len(manip_obj.sources)>=0 and len(manip_obj.targets)>=0):
                    lib_path_finding.find_paths(manip_obj)
                
        if auto_manipulate:
            logging.info(lib_utils.log_message("Structure recognition stopped."))
            
    manip_obj.t1 = threading.Thread(target = do_this, name = 'do_this')
    manip_obj.t1.start()
    
def create_site_and_connect_to_point_region(manip_obj, loc, site_id, point_region):
    manip_obj.sites.append( aab.Site(loc[0], loc[1], site_id=site_id) )
    # Mutual assignment
    manip_obj.sites[-1].graphic = point_region
    point_region.site = manip_obj.sites[-1]
            
# Auto-detect and display foreign atoms.
def func_auto_detect_foreign_atoms(sr_obj):
    # aliases
    manip_obj = sr_obj.manip_obj
    rra = manip_obj.rectangle_regions_auto # "rra" means rectangle regions auto-detected
    pdi = manip_obj.processed_data_item
    shape = manip_obj.source_xdata.data_shape
    
    # Delete old auto-detected sources.
    atoms = copy.copy(manip_obj.sources)
    for atom in atoms:
        if not atom.defined_by_user:
            manip_obj.sources.remove(atom)
    
    foreigns_site_id = []
    if sr_obj.auto_detect_foreign_atoms:
        foreigns_site_id = np.nonzero(sr_obj.nn_output['labels'] == 1)[0]
        logging.info(lib_utils.log_message(f"Detected {len(foreigns_site_id):d} foreign atoms."))
      
    relative_size = 0.05
    number_foreigns = len(foreigns_site_id)
    
    while not manip_obj.rdy_update_pdi.wait(1):
        logging.info(lib_utils.log_message("Waiting for an update of the processed data item..."))
        pass
    
    new_centers = []    
    for i in range(number_foreigns):
        site_id = foreigns_site_id[i]
        loc = manip_obj.maxima_locations[site_id]
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
                
        manip_obj.sources.append( aab.Atom(manip_obj.sites[site_id], 'pseudo-element') )
        # Mutual assignment
        rra[i].atom = manip_obj.sources[-1]
        manip_obj.sources[-1].graphic = rra[i]
        if added:
            lib_utils.add_listener_graphic_changed(manip_obj, rra[i])
            
    if len(new_centers) > 0:
        for i, center in enumerate(new_centers):
            rra[i].center = center
        
    if len(rra) > number_foreigns: # Remove graphics
        #print(" ! Removing excess RECTANGLE graphics ")
        manip_obj.rectangle_regions_auto = rra[0:number_foreigns]
        def func():
            for k in range(number_foreigns, len(rra)):
                pdi.remove_region(rra[k])
        manip_obj.api.queue_task(func)
    
    def do_this(): # The following is threaded, because it is not needed for later processes.
        lib_utils.refresh_GUI(manip_obj, ['foreigns'])
        lib_utils.elemental_identification(manip_obj)
    threading.Thread(target=do_this).start()

def clear_user_defined_atoms_and_targets(manip_obj):
    pdi = manip_obj.processed_data_item
    def func():
        with manip_obj.api.library.data_ref_for_data_item(pdi):
            for region in manip_obj.rectangle_regions:
                pdi.remove_region(region)
            for region in manip_obj.ellipse_regions:
                pdi.remove_region(region)
        manip_obj.rectangle_regions = []
        manip_obj.ellipse_regions = []
    manip_obj.api.queue_task(func)
    manip_obj.sources = []
    manip_obj.targets = []

    logging.info(lib_utils.log_message("Cleared all user-defined foreign atoms and target sites."))
    lib_utils.refresh_GUI(manip_obj, ['foreigns', 'targets'])