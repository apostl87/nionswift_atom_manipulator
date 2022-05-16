import gettext
import threading
import numpy as np

# GUI classes
from .gui_tractor_beam import TractorBeamModule
from .gui_path_finding import PathFindingModule
from .gui_structure_recognition import StructureRecognitionModule
from .gui_manipulation import ManipulationModule
from .lib_widgets import ScrollArea, push_button_template
from . import lib_utils

_ = gettext.gettext

defaults = {'simulation_mode': False}

devices_dict = {False: 'scan_controller',
                True: 'usim_scan_device'}

class AtomManipulatorDelegate:
     
    def __init__(self, api):
        self.api = api
        self.panel_id = "atom-manipulator-panel"
        self.panel_name = _("Atom Manipulator")
        self.panel_positions = ["left", "right"]
        self.panel_position = "right"
        
        # Init
        self.simulation_mode = None
        self.superscan = None
        self.scan_parameters = None
        self.scan_parameters_changed = None
        self.snapshot_counter = None
        
        # Nion Swift objects that are used by this plug-in.
        self.source_xdata = None
        self.processed_data_item = None
        
        # Objects internal to the plug-in.
        self.sites = []
        self.sources = []
        self.targets = []
        self.bonds = None
        self.paths = None
        
        # Workaround: Reposition graphics.
        self.listeners = []
        self.point_regions = []
        self.line_regions = []
        self.rectangle_regions = [] 
        self.rectangle_regions_auto = []
        self.ellipse_regions = []
        
        # Threads.
        self.t1 = None
        self.t5 = None
        self.t6 = None
        
        # Events.
        self.sr_rdy = threading.Event()
        self.pf_rdy = threading.Event()
        self.tb_rdy = threading.Event()
        self.rdy_create_pdi = threading.Event()
        self.rdy_create_pdi.set()
        self.rdy_init_pdi = threading.Event()
        self.rdy_init_pdi.set()
        self.rdy_update_pdi = threading.Event()
        self.rdy_update_pdi.set()
    
    def clear_manipulator_objects(self):
        self.sites = []
        self.sources = []
        self.targets = []
        self.bonds = None
        self.paths = None
        self.listeners = []
        self.point_regions = []
        self.line_regions = []
        self.rectangle_regions = [] 
        self.rectangle_regions_auto = []
        self.ellipse_regions = []

    def create_panel_widget(self, ui, document_controller):
        self.ui = ui
        self.document_controller = document_controller
        
        # Callback functions.
        def simulation_mode_changed(checked):
            self.superscan = self.api.get_hardware_source_by_id(devices_dict[checked], "1")
            self.simulation_mode_checkbox.checked = checked
            self.simulation_mode = checked
        def create_new_data_item_clicked():
            lib_utils.create_pdi(self)

        # GUI init.
        main_col = ui.create_column_widget()
        scroll_area = ScrollArea(ui._ui)
        scroll_area.content = main_col._widget
        
        # GUI elements.
        feature_row, self.data_item_button = push_button_template(
            self.ui, _("Create new data item"), callback=create_new_data_item_clicked)
        feature_row.add_stretch()
        
        self.simulation_mode_checkbox = ui.create_check_box_widget("Simulation mode")
        self.simulation_mode_checkbox.on_checked_changed = simulation_mode_changed
        feature_row.add(self.simulation_mode_checkbox)
        
        # Modules.
        self.structure_recognition_module = StructureRecognitionModule(self.ui, self.api, self.document_controller, self)
        self.path_finding_module = PathFindingModule(self.ui, self.api, self.document_controller, self)
        self.tractor_beam_module = TractorBeamModule(self.ui, self.api, self.document_controller, self)
        self.manipulation_module = ManipulationModule(self.ui, self.api, self.document_controller, self)
        
        # Build main column.
        main_col.add_spacing(5)
        main_col.add(feature_row)
        self.structure_recognition_module.create_widgets(main_col)
        self.path_finding_module.create_widgets(main_col)
        self.tractor_beam_module.create_widgets(main_col)
        self.manipulation_module.create_widgets(main_col)
        main_col.add_stretch()
        
        # Set defaults.
        simulation_mode_changed(defaults['simulation_mode'])
        
        #return main_col
        return scroll_area 
        

class AtomManipulatorExtension(object):
    # Required for Nion Swift to recognize this as an extension class.
    extension_id = "nion.swift.extension.atom_manipulator"
    
    def __init__(self, api_broker):
        # Grab the API object.
        api = api_broker.get_api(version='~1.0', ui_version='~1.0')
        # Be sure to keep a reference or it will be closed immediately.
        self.__panel_ref = api.create_panel(AtomManipulatorDelegate(api))
  
    def close(self):
        # Close will be called when the extension is unloaded. In turn, close any references so they get closed.
        # This is not strictly necessary since the references will be deleted naturally when this object is deleted.
        self.__panel_ref.close()
        self.__panel_ref = None 
