import gettext
import threading

import logging

# Custom libraries
from . import lib_utils
from . import lib_structure_recognition
from . import lib_pathfinding
from .lib_widgets import Section, line_edit_template, check_box_template, combo_box_template
from adf_feedback import adf_feedback as adffb

_ = gettext.gettext

# Defaults on initialization.
defaults = {'snapshot': True}


class ManipulationModule(lib_utils.AtomManipulatorModule):

    def __init__(self, ui, api, document_controller, manipulator):
        super().__init__(ui, api, document_controller)
        self.manipulator = manipulator # AtomManipulatorDelegate object
        self.snapshot = None
        
        # Events.
        self.stop_auto_manipulate_event = threading.Event()

    # GUI creation method. 
    def create_widgets(self, column):
        self.section = Section(self.ui, 'Manipulation')
        column.add(self.section)
        
        # Callback functions.
        def snapshot_changed(checked):
            self.snapshot_row_checkbox.checked = checked
            self.snapshot = checked
        
        ## GUI elements.

        # Snapshot row.
        snapshot_row = self.ui.create_row_widget()
        self.snapshot_row_checkbox = self.ui.create_check_box_widget(_("Snapshot every frame (1x raw, 1x w/ insets)"))
        self.snapshot_row_checkbox.on_checked_changed = snapshot_changed
        snapshot_row.add_spacing(10)
        snapshot_row.add(self.snapshot_row_checkbox)
        
        # Button row (Start/Stop).
        automanip_button = self.ui.create_push_button_widget(_("Start automatical manipulation"))
        automanip_button.state = False
        def automanip_button_clicked():
            automanip_button.state = not automanip_button.state
            if automanip_button.state:
                logging.info(lib_utils.log_message("Starting automatical manipulation ..."))
                ## Init
                self.manipulator.tractor_beam_module.rdy.set()
                self.manipulator.pathfinding_module.rdy.clear()
                self.manipulator.structure_recognition_module.rdy.clear()
                self.stop_auto_manipulate_event.clear()
                
                ## Call modules
                # Structure recognition.
                if self.snapshot:
                    self.manipulator.snapshot_counter = 1
                else:
                    self.manipulator.snapshot_counter = None
                lib_structure_recognition.analyze_and_show(self.manipulator.structure_recognition_module, auto_manipulate=True)

                # Pathfinding.
                lib_pathfinding.find_paths(self.manipulator, auto_manipulate=True)
                
                # Tractor beam (wrapper).
                try:
                    # In this (auto-manipulation) operation mode, ADFFeedback itself listens to comm_events and will stop execution accordingly.
                    
                     # Communication events.
                    comm_events = [self.manipulator.structure_recognition_module.new_image,
                                   self.manipulator.pathfinding_module.rdy,
                                   self.manipulator.tractor_beam_module.rdy]
                    
                    self.ADFFeedback = adffb.ADFFeedbackDelegate(
                        self.api, 
                        offline_test_mode = self.manipulator.tractor_beam_module.otm_check_box.checked,
                        autodetect_changes = self.manipulator.tractor_beam_module.adc_check_box.checked,
                        external_stop_event = self.stop_auto_manipulate_event
                    )
                    
                    self.ADFFeedback.startmap(
                        self.manipulator.tractor_beam_module.frame_timeout,
                        self.manipulator.tractor_beam_module.jump_threshold,
                        self.manipulator.tractor_beam_module.drift_threshold,
                        self.manipulator.tractor_beam_module.reposition_timeout,
                        comm_events=comm_events,
                        manip_obj=self.manipulator
                    )
                        
                except Exception as exc:
                    logging.info(lib_utils.log_message("Exception when calling adf_feedback"))
                    print(exc)

                # Change button text.
                automanip_button.text = _('Stop Auto-Manipulator')
            else:
                self.stop_auto_manipulate_event.set()
                automanip_button.text = _('Start Auto-Manipulator')     
        automanip_button.on_clicked = automanip_button_clicked
        automanip_button_row = self.ui.create_row_widget()
        automanip_button_row.add(automanip_button)
        
        # Set defaults.
        snapshot_changed(defaults['snapshot'])
        
        # Assemble GUI elements.
        self.section.column.add(snapshot_row)
        self.section.column.add(automanip_button_row)