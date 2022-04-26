import string
import numpy as np
import copy
from matplotlib import colors as mcolors
from skimage import draw
import logging
import gettext
try: from periodictable import elements as pt_elements
except: pass
import math
from nionswift_plugin.double_gaussian_filter.computation import double_gaussian_blur as dgb

_ = gettext.gettext

def sub2ind(rows, cols, array_shape):
    return rows * array_shape[1] + cols

def ind2sub(array_shape, ind):
    rows = (np.int32(ind) // array_shape[1])
    cols = (np.int32(ind) % array_shape[1])
    return (rows, cols)

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
        
def plot_points(image, points, size=3, color="blue"):
    if points is None:
        return image
    color = (255, 165, 0)
    points = np.round(points).astype(np.int)
    for point in points:
        try:
            image[draw.disk((point[0], point[1]), size, shape=image.shape)] = color
        except IndexError:
            pass
    return image     

def plot_paths(image, paths):
    color = (0, 165, 255)
    for path in paths.members:
        for i in range(len(path.sitelist)-1):
            y1, x1 = path.sitelist[i].coords
            y2, x2 = path.sitelist[i+1].coords
            # shorten the display of the bond
            w = 1/5 # weight of position A
            y1 = (y1*(1-w) + y2*w).round().astype(int)
            x1 = (x1*(1-w) + x2*w).round().astype(int)
            y2 = (y1*w + y2*(1-w)).round().astype(int)
            x2 = (x1*w + x2*(1-w)).round().astype(int)
            try:
                image[draw.line(y1, x1, y2, x2)] = color
            except:
                pass
    return image

def refresh_GUI(manip_obj, var_strings):
    def func():
        if 'atoms' in var_strings:
            manip_obj.structure_recognition_module.N_atoms_label.text = str(len(manip_obj.sites))
        if 'foreigns' in var_strings:
            manip_obj.path_finding_module.N_foreign_atoms_label.text = str(len(manip_obj.sources))
        if 'targets' in var_strings:
            manip_obj.path_finding_module.N_target_sites_label.text = str(len(manip_obj.targets))
        if 'sampling' in var_strings:
            if manip_obj.structure_recognition_module.sampling is None:
                manip_obj.structure_recognition_module.sampling_label.text = \
                    f"N/A"
            else:
                manip_obj.structure_recognition_module.sampling_label.text = \
                    f"{manip_obj.structure_recognition_module.sampling:4f} Å/px"
            if manip_obj.structure_recognition_module.fov is None:
                manip_obj.structure_recognition_module.fov_label.text = \
                    f"N/A"
            else:
                manip_obj.structure_recognition_module.fov_label.text = \
                    f"{manip_obj.structure_recognition_module.fov[0]:.2f} x {manip_obj.structure_recognition_module.fov[1]:.2f} Å^2"
                    
    manip_obj.api.queue_task(func)
    
# GUI task function for creating a new processed data item without data
def create_pdi(manip_obj):
    if manip_obj.t1 and manip_obj.t1.is_alive():
        logging.info("Cannot create new data item while AtomManipulator is running")
        return None
    dummy_data = np.zeros((1,1,3), dtype=np.uint8)
    xdata = manip_obj.api.create_data_and_metadata(dummy_data)
    def func():
        try:
            manip_obj.processed_data_item.title = manip_obj.processed_data_item.title[7:]
        except:
            pass                                                                       
        manip_obj.processed_data_item = manip_obj.document_controller.create_data_item_from_data_and_metadata(
                           xdata, title=_('[LIVE] ') + ('AtomManipulator_') + _('dummy'))
        manip_obj.rdy_create_pdi.set()
        manip_obj.clear_manipulator_objects()
        refresh_GUI(manip_obj, ['atoms', 'foreigns', 'targets'])
    manip_obj.api.queue_task(func)

# GUI task function to be called after new image has been read.
def init_pdi(manip_obj):
    def func():
        if manip_obj.processed_data_item not in manip_obj.api.library.data_items:
            manip_obj.rdy_create_pdi.clear()
            create_pdi(manip_obj)
        while not manip_obj.rdy_create_pdi.wait(1):
            pass # waiting for creation of pdi
        manip_obj.processed_data_item.title = _('[LIVE] ') + 'AtomManipulator_' + manip_obj.source_title
        manip_obj.processed_data_item.xdata = copy.deepcopy(manip_obj.source_xdata)
        
        data = np.array(manip_obj.processed_data_item.data)
        manip_obj.processed_data_item.original_data_rgb = np.tile(
            ((data - data.min()) / data.ptp() * 255).astype(np.uint8)[..., None], (1, 1, 3))
        manip_obj.processed_data_item.original_data = data
        
        manip_obj.processed_data_item.data = manip_obj.processed_data_item.original_data_rgb

        if manip_obj.snapshot_counter is not None:
            sdi = manip_obj.api.library.snapshot_data_item(manip_obj.processed_data_item)
            sdi.title = _('AtomManipulator frame ' + str(manip_obj.snapshot_counter) + ' RAW_' + manip_obj.source_title)
            manip_obj.snapshot_counter += 1

        manip_obj.rdy_init_pdi.set()
    manip_obj.api.queue_task(func)
    
# GUI task function for updating the data in the processed_data_item
def update_pdi(manip_obj, new_data):
    def func():
        if manip_obj.processed_data_item not in manip_obj.api.library.data_items:
            manip_obj.rdy_init_pdi.clear()
            init_pdi()
        while not manip_obj.rdy_init_pdi.wait(1):
            pass
        manip_obj.processed_data_item.data = new_data
        manip_obj.rdy_update_pdi.set()
    manip_obj.api.queue_task(func)
    
# Supoort function: Add a listener to "graphic changed" events
def add_listener_graphic_changed(manip_obj, graphic):
    def check_site():
        shape = manip_obj.source_xdata.data_shape
        all_coords = map(lambda x: x.coords, manip_obj.sites)
        nearest_site = manip_obj.sites[ np.linalg.norm(
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
    #graphic = manip_obj.api.library.get_data_item_by_uuid(manip_obj.processed_data_item.uuid).graphics[-1]
    manip_obj.listeners.append( graphic._graphic.graphic_changed_event.listen(check_site) )
    
# Elemental identification
def elemental_identification(manip_obj):
    # Calculate intensity values
    while not manip_obj.rdy_create_pdi.wait(1) or not manip_obj.rdy_init_pdi.wait(1) or not manip_obj.rdy_update_pdi.wait(1):
        pass # waiting for tasks on pdi to be completed
    
    # Aliases
    sampling = manip_obj.structure_recognition_module.sampling
    labels = manip_obj.structure_recognition_module.nn_output['labels']
    int_radius_A = manip_obj.structure_recognition_module.elemental_id_int_radius
    Z_exponent = manip_obj.structure_recognition_module.elemental_id_exponent
    
    # Get double Gaussian blur
    if np.isnan(sampling):
        print("Elemental identfication cannot be perfomed, because there is no sampling value [Angstroem/px] available.")
        return

    sigma1 = 0.25 # in Angstroems
    sigma1 /= sampling
    data = dgb(manip_obj.processed_data_item.original_data, sigma1=sigma1, sigma2=3*sigma1, weight2=0.4)
               
    intensities = integrate_intensities(data, manip_obj.maxima_locations,
                integration_radius = int_radius_A/sampling
                )
    mean_intensity_carbon = np.nanmean( intensities[labels == 0] )
    Z_carbon = 6
        
    labels = []
    graphics = []
    for atom in manip_obj.sources:
        if atom.defined_by_user:
            #continue
            pass
        intensity = intensities[atom.site.id]
        if math.isnan(intensity):
            labels.append( "n.a." )
        else:
            Z_estimator = (intensity / mean_intensity_carbon) ** (1/Z_exponent) * Z_carbon
            try: labels.append( pt_elements[round(Z_estimator)] )
            except: labels.append( str(round(Z_estimator,1)) )
        graphics.append( atom.graphic )

    def func():
        for label, graphic in zip(labels, graphics):
            graphic.label = label
    manip_obj.api.queue_task(func)

def integrate_intensities(data, maxima_locations, integration_radius=1):
    # data ... image data (numpy.ndarray)
    # integration_radius ... (scalar)
    # maxima_locations ... (N x 2 numpy.ndarray)
    
    # Conditioning
    if type(maxima_locations) is not np.ndarray:
        maxima_locations = np.array(maxima_locations)
    if type(data) is not np.ndarray:
        data = np.array(data)
    shape = data.shape
    N = maxima_locations.shape[0]
    
    # Init array for intensity values
    values = np.full(N, np.nan)
    
    # Create mask
    integration_radius_floor = math.floor(integration_radius) # integer value; NOT making the area of integration smaller
    M = 2*integration_radius_floor + 1
    mask = np.ones((M, M))
    for i in range(-integration_radius_floor, integration_radius_floor+1):
        for j in range(-integration_radius_floor, integration_radius_floor+1):
            if math.sqrt(i**2 + j**2) > integration_radius:
                mask[i+integration_radius_floor, j+integration_radius_floor] = np.nan
    
    # Background subtraction
    #TODO check
    intensity_min = np.min(data)
    data_shifted = data - intensity_min
    
    for i in range(N):
        # Location of a maximum
        loc = np.array(maxima_locations[i,:]+0.5, dtype=int) # integer values; pixel position in the image
        
        # Check if integration range is non-negative and smaller than shape, no try-block
        lower = loc-integration_radius_floor
        upper = loc+integration_radius_floor+1
        if any(lower < 0) or any(upper > np.array(shape)):
            #TODO For now, just skip that
            continue
        else: # Do the integration
            tmp = data_shifted[lower[0]:upper[0], lower[1]:upper[1]] # Shift minimum to 0
            values[i] = np.nanmean( tmp*mask )
            
    return values

def log_message(input: string):
    return f"Atom Manipulator: {input} "#.center(80, '-')