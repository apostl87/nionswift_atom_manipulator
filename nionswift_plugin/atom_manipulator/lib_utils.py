"""
Utilities library.
"""

import gettext
import numpy as np

import copy
import logging

import math
import string

from matplotlib import colors as mcolors
from skimage import draw

# Non-standard packages
try:
    from periodictable import elements as pt_elements # Optional
except:
    pass
from double_gaussian_blur import dgb

_ = gettext.gettext


# Convert two-dimensional indices to one-dimensional index.
def sub2ind(rows, cols, array_shape):
    return rows * array_shape[1] + cols


# Convert one-dimensional index to two-dimenionsal indices.
def ind2sub(array_shape, ind):
    rows = (np.int32(ind) // array_shape[1])
    cols = (np.int32(ind) % array_shape[1])
    return (rows, cols)


# Template class for modules.
class AtomManipulatorModule(object):

    def __init__(self, ui, api, document_controller):
        self.ui = ui
        self.api = api
        self.document_controller = document_controller

    def create_widgets(self, column):
        raise NotImplementedError()

    def set_preset(self, name):
        raise NotImplementedError()

    def fetch_parameters(self):
        raise NotImplementedError()


# Insert points into image.
def plot_points(image, points, size=3, color="blue"):
    if points is None:
        return image
    color = (255, 165, 0)
    points = np.round(points).astype(int)
    for point in points:
        try:
            image[draw.disk((point[0], point[1]), size, shape=image.shape)] = color
        except IndexError:
            pass
    return image     


# Insert paths into image.
def plot_paths(image, paths):
    color = (0, 165, 255)
    for path in paths.members:
        for i in range(len(path.sitelist)-1):
            y1, x1 = path.sitelist[i].coords
            y2, x2 = path.sitelist[i+1].coords
            # Shorten the display of the bond.
            w = 1/5 # Weight of position A.
            y1 = (y1*(1-w) + y2*w).round().astype(int)
            x1 = (x1*(1-w) + x2*w).round().astype(int)
            y2 = (y1*w + y2*(1-w)).round().astype(int)
            x2 = (x1*w + x2*(1-w)).round().astype(int)
            try:
                image[draw.line(y1, x1, y2, x2)] = color
            except:
                pass
    return image


def refresh_GUI(manipulator, var_strings):
    def func():
        if 'atoms' in var_strings:
            manipulator.structure_recognition_module.N_atoms_label.text = str(len(manipulator.sites))
        if 'foreigns' in var_strings:
            manipulator.pathfinding_module.N_foreign_atoms_label.text = str(len(manipulator.sources))
        if 'targets' in var_strings:
            manipulator.pathfinding_module.N_target_sites_label.text = str(len(manipulator.targets))
        if 'sampling' in var_strings:
            if manipulator.structure_recognition_module.sampling is None:
                manipulator.structure_recognition_module.sampling_label.text = \
                    f"N/A"
            else:
                manipulator.structure_recognition_module.sampling_label.text = \
                    f"{manipulator.structure_recognition_module.sampling:4f} Å/px"
            if manipulator.structure_recognition_module.fov is None:
                manipulator.structure_recognition_module.fov_label.text = \
                    f"N/A"
            else:
                manipulator.structure_recognition_module.fov_label.text = \
                    f"{manipulator.structure_recognition_module.fov[0]:.2f} x " \
                    f"{manipulator.structure_recognition_module.fov[1]:.2f} Å^2"
                    
    manipulator.api.queue_task(func)


# GUI task function for creating a new processed data item without data.
def create_pdi(manipulator):
    if manipulator.t1 and manipulator.t1.is_alive():
        logging.info("Cannot create new data item while AtomManipulator is running")
        return None
    dummy_data = np.zeros((1,1,3), dtype=np.uint8)
    xdata = manipulator.api.create_data_and_metadata(dummy_data)
    def func():
        try:
            manipulator.processed_data_item.title = manipulator.processed_data_item.title[7:]
        except:
            pass                                                                       
        manipulator.processed_data_item = manipulator.document_controller.create_data_item_from_data_and_metadata(
                           xdata, title=_('[LIVE] ') + ('AtomManipulator_') + _('dummy'))
        manipulator.rdy_create_pdi.set()
        manipulator.clear_manipulator_objects()
        refresh_GUI(manipulator, ['atoms', 'foreigns', 'targets'])
    manipulator.api.queue_task(func)


# GUI task function to be called after new image has been read.
def init_pdi(manipulator):
    def func():
        if manipulator.processed_data_item not in manipulator.api.library.data_items:
            manipulator.rdy_create_pdi.clear()
            create_pdi(manipulator)
        while not manipulator.rdy_create_pdi.wait(1):
            pass # Waiting for creation of processed_data_item.
        manipulator.processed_data_item.title = _('[LIVE] ') + 'AtomManipulator_' + manipulator.source_title
        xdata = copy.deepcopy(manipulator.source_xdata)
        xdata.metadata[manipulator.metadata_root_key] = manipulator.metadata_to_append
        
        # Snapshot RAW data if checkbox is checked
        if manipulator.snapshot_counter is not None:
            manipulator.processed_data_item.set_data_and_metadata(xdata)
            with manipulator.api.library.data_ref_for_data_item(manipulator.processed_data_item):
                sdi = manipulator.api.library.snapshot_data_item(manipulator.processed_data_item)
            sdi.title = _('AtomManipulator frame ' + str(manipulator.snapshot_counter) +
                          ' RAW_' + manipulator.source_title)
            manipulator.snapshot_counter += 1

        # Convert data to RGB values, save original data as well as rgb data in data item
        data = np.array(xdata.data)
        manipulator.processed_data_item.original_data = data
        manipulator.processed_data_item.rgb_data = np.tile(
            ((data - data.min()) / data.ptp() * 255).astype(np.uint8)[..., None], (1, 1, 3))
        
        # Set data and metadata of data item
        manipulator.processed_data_item.set_data(manipulator.processed_data_item.rgb_data)
        manipulator.processed_data_item.set_metadata(xdata.metadata)

        manipulator.rdy_init_pdi.set()
    manipulator.api.queue_task(func)


# GUI task function for updating the data in the processed_data_item.
def update_pdi(manipulator, new_data):
    def func():
        if manipulator.processed_data_item not in manipulator.api.library.data_items:
            manipulator.rdy_init_pdi.clear()
            init_pdi()
        while not manipulator.rdy_init_pdi.wait(1):
            pass

        metadata = copy.deepcopy(manipulator.processed_data_item.metadata)
        metadata[manipulator.metadata_root_key] = manipulator.metadata_to_append
        
        manipulator.processed_data_item.set_data(new_data)
        manipulator.processed_data_item.set_metadata(metadata)
        
        manipulator.rdy_update_pdi.set()
    manipulator.api.queue_task(func)


# Support function: add a listener to "graphic changed" events.
def add_listener_graphic_changed(manipulator, graphic):
    def check_site():
        shape = manipulator.source_xdata.data_shape
        all_coords = map(lambda x: x.coords, manipulator.sites)
        nearest_site = manipulator.sites[ np.linalg.norm(
                np.array(list(all_coords))-np.array(graphic.center)*shape, axis=1).argmin() ]
        #print(nearest_site)
        if hasattr(graphic, 'atom') and (nearest_site is not graphic.atom.origin):
            graphic.atom.site = nearest_site
            graphic.atom.origin = nearest_site
            #logging.info(" changed site of foreign atom ")
        elif hasattr(graphic, 'site') and (nearest_site is not graphic.site):
            graphic.site = nearest_site
            #logging.info(" changed target site ")
        else:
            pass
    #graphic = manipulator.api.library.get_data_item_by_uuid(manipulator.processed_data_item.uuid).graphics[-1]
    manipulator.listeners.append( graphic._graphic.property_changed_event.listen(check_site) )


# Element identification.
def element_identification(manipulator):
    # Calculate intensity values.
    while not manipulator.rdy_create_pdi.wait(1) or not manipulator.rdy_init_pdi.wait(1) or not manipulator.rdy_update_pdi.wait(1):
        pass # waiting for tasks on processed_data_item to be completed
    
    # Aliases.
    sampling = manipulator.structure_recognition_module.sampling
    labels = manipulator.structure_recognition_module.nn_output['labels']
    int_radius_A = manipulator.structure_recognition_module.element_id_int_radius
    Z_exponent = manipulator.structure_recognition_module.element_id_exponent
    
    # Get double Gaussian blur.
    if np.isnan(sampling):
        print("Element identfication cannot be perfomed, because there is no sampling value [Angstroem/px] available.")
        return

    sigma1 = 0.25 # in Angstroem
    sigma1 /= sampling
    data = dgb(manipulator.processed_data_item.original_data, sigma1=sigma1, sigma2=3*sigma1, weight2=0.4)
               
    intensities = integrate_intensities(data, manipulator.maxima_locations,
                integration_radius = int_radius_A/sampling
                )
    mean_intensity_carbon = np.nanmean( intensities[labels == 0] )
    Z_carbon = 6
        
    labels = []
    graphics = []
    for atom in manipulator.sources:
        if atom.defined_by_user:
            #continue
            pass
        intensity = intensities[atom.site.id]
        if math.isnan(intensity):
            labels.append( "n.a." )
        else:
            Z_estimator = (intensity / mean_intensity_carbon) ** (1/Z_exponent) * Z_carbon
            try:
                labels.append( pt_elements[round(Z_estimator)] )
            except:
                labels.append( str(round(Z_estimator,1)) )
        graphics.append( atom.graphic )

    def func():
        for label, graphic in zip(labels, graphics):
            #graphic.label = label ## display of Z disabled
            pass
    manipulator.api.queue_task(func)

# Helper function for element identification
def integrate_intensities(data, maxima_locations, integration_radius=1):
    # data ... image data (numpy.ndarray)
    # integration_radius ... (scalar)
    # maxima_locations ... (N x 2 numpy.ndarray)

    # Conditioning of inputs.
    if type(data) is not np.ndarray:
        data = np.array(data)
    if type(maxima_locations) is not np.ndarray:
        maxima_locations = np.array(maxima_locations)

    # Aliases.
    shape = data.shape
    N = maxima_locations.shape[0]
    
    # Init array for intensity values.
    values = np.full(N, np.nan)
    
    # Create mask.
    integration_radius_floor = math.floor(integration_radius) # integer value; NOT making area of integration smaller.
    M = 2*integration_radius_floor + 1
    mask = np.ones((M, M))
    for i in range(-integration_radius_floor, integration_radius_floor+1):
        for j in range(-integration_radius_floor, integration_radius_floor+1):
            if math.sqrt(i**2 + j**2) > integration_radius:
                mask[i+integration_radius_floor, j+integration_radius_floor] = np.nan
    
    # Background subtraction.
    intensity_min = np.min(data)
    data_shifted = data - intensity_min
    
    for i in range(N):
        # Location of a maximum.
        loc = np.array(maxima_locations[i,:]+0.5, dtype=int) # Integer values; pixel position in the image.
        
        # Check if integration range is non-negative and smaller than shape, no try-block.
        lower = loc-integration_radius_floor
        upper = loc+integration_radius_floor+1
        if any(lower < 0) or any(upper > np.array(shape)):
            #TODO For now, just skip that
            continue
        else: # Do the integration.
            tmp = data_shifted[lower[0]:upper[0], lower[1]:upper[1]] # Shift minimum to 0.
            values[i] = np.nanmean( tmp*mask )
            
    return values


# Uniform log messages.
def log_message(input: string):
    return f"Atom Manipulator: {input} "#.center(80, '-')
