import gettext
import threading
import numpy as np

import time
import logging

# Non-standard packages
from nionswift_plugin.nionswift_structure_recognition.model import load_preset_model
from nionswift_plugin.nionswift_structure_recognition.scale import RealSpaceCalibrator

# Custom libraries
from . import lib_utils
from . import lib_structure_recognition
from .lib_utils import AtomManipulatorModule
from .lib_widgets import Section, line_edit_template, check_box_template, combo_box_template, push_button_template

_ = gettext.gettext

# Defaults on initialization.
defaults = {'visualize_atoms': True,
            'auto_detect_foreign_atoms': False,
            'elemental_identification_integration_radius_A': 0.25, # in Angstroem
            'elemental_identification_exponent': 1.8,
            'image_source': 0,   # 0: MAADF, 1: HAADF, 2: Selected data item
            'scale_calibration_mode': 1 # 0: Manual, 1: Live
            }


class StructureRecognitionModule(AtomManipulatorModule):

    def __init__(self, ui, api, document_controller, manipulator):
        super().__init__(ui, api, document_controller)
        self.manipulator = manipulator # AtomManipulatorDelegate object
        self.auto_detect_foreign_atoms = None
        self.elemental_id_int_radius = None
        self.elemental_id_exponent = None
        self.nn_output = None
        self.scale_calibration_mode = None
        self.sampling = None
        self.fov = None
        self.visualize_atoms = None
        self.live_analysis = None
        self.was_playing = None
        self.stop_live_analysis_event = threading.Event()
        self.rdy = threading.Event()
        self.new_image = threading.Event()
        self.model = load_preset_model('graphene')
    
    # GUI creation method.
    def create_widgets(self, column):
        section1 = Section(self.ui, 'Scale calibration')     # Scale calibration recognition section.
        section2 = Section(self.ui, 'Structure recognition') # Structure recognition section.
        column.add(section1)
        column.add(section2)
        
        # Callback functions.
        def scale_calibration_mode_changed(item):
            # Sanitization of input.
            if type(item) == int:
                item = self.scale_calibration_combo_box.items[item]
            # Ensure combo box displays the correct mode.
            self.scale_calibration_combo_box.current_item = item

            # Remove widgets if they are in the container.
            try:
                scale_calibration_button_container._widget.remove_all()
            except:
                pass

            # Build widget container.
            if self.scale_calibration_combo_box.current_index == 0: # Single data item for RealSpaceCalibrator.
                self.scale_calibration_mode = 0
                dummy_row, self.scale_calibration_button = push_button_template(self.ui, 'Run with selected data item')
                self.scale_calibration_button.on_clicked = scale_calibration_button_clicked
                scale_calibration_button_container.add(self.scale_calibration_button)
            elif self.scale_calibration_combo_box.current_index == 1: # Live calibration with FourierScaleCalibrator.
                self.scale_calibration_mode = 1
                self.sampling = None
                self.fov = None
                lib_utils.refresh_GUI(self.manipulator, ['sampling'])

        def scale_calibration_button_clicked():
            # Calculate target values that directly follow from Nion Swift Scan settings.
            tdi = self.manipulator.document_controller.target_data_item
            fov_1d_target_value_nm = tdi.metadata['scan']['fov_nm']
            sampling_target_value = fov_1d_target_value_nm*10/np.sqrt(tdi.data.size)
            min_sampling = sampling_target_value*.75
            max_sampling = sampling_target_value*1.25
            step_size=(max_sampling-min_sampling)/5
        
            calibrator = RealSpaceCalibrator(model=self.model,
                               template='hexagonal',
                               lattice_constant=2.46,
                               min_sampling=min_sampling,
                               max_sampling=max_sampling,
                               step_size=step_size
            )
            
            def do_this():
                t = time.time()
                logging.info(lib_utils.log_message("Calling RealSpaceCalibrator"))
                self.sampling = calibrator(tdi.data)
                t = time.time()-t
                logging.info(lib_utils.log_message(f"RealSpaceCalibrator finished after {t:.5f} seconds"))
                
                self.fov = [self.sampling*s for s in tdi.data.shape]
                lib_utils.refresh_GUI(self.manipulator, ['sampling'])
            # Run calibration, but not in main thread.
            self.manipulator.t6 = threading.Thread(target = do_this, name = 'RealSpaceCalibrator')
            self.manipulator.t6.start()

        def image_source_changed(item):
            if type(item) == int:
                item = self.image_source_combo_box.items[item]
            self.image_source_combo_box.current_item = item
            
            # Save whether live analysis is checked when switching between MAADF and HAADF.
            try:
                saved_var = self.live_analysis_check_box.checked
            except:
                saved_var = False    
            live_analysis_row._widget.remove_all()

            if self.image_source_combo_box.current_index in [0, 1]:
                self.live_analysis_check_box = self.ui.create_check_box_widget('Live analysis')
                live_analysis_row.add(self.live_analysis_check_box)
                self.live_analysis_check_box.on_checked_changed = live_analysis_changed
                self.live_analysis_check_box.checked = saved_var
            live_analysis_changed(saved_var) # Calling this ensures that the push button goes into the correct state.

        def visualize_atoms_changed(checked):
            self.visualize_atoms = checked

        def auto_detect_foreign_atoms_changed(checked):
            self.auto_detect_foreign_atoms = checked

        def elemental_id_int_radius_changed(text):
            if len(text) > 0:
                try:
                    self.elemental_id_int_radius = float(text)
                except:
                    pass
                finally:
                    self.elemental_id_int_radius_line_edit.text = f"{self.elemental_id_int_radius:.2f}"

        def elemental_id_exponent_changed(text):
            if len(text) > 0:
                try:
                    self.elemental_id_exponent = float(text)
                except: pass
                finally: self.elemental_id_exponent_line_edit.text = f"{self.elemental_id_exponent:.2f}"

        def live_analysis_changed(checked):
            self.stop_live_analysis_event.set() # always stop live analysis
            start_stop_analysis_button_next_state(state = None)   
        
        # Start procedure depending on the state of the button.
        def start_stop_analysis():            
            state = self.start_stop_analysis_button.state
            if state == 0: # Start single image analysis.
                self.was_playing = None
                self.stop_live_analysis_event.clear()
                lib_structure_recognition.analyze_and_show(self, live_analysis=False)
            elif state == 2: # Start live analysis.
                self.was_playing = self.manipulator.superscan.is_playing
                self.stop_live_analysis_event.clear()
                lib_structure_recognition.analyze_and_show(self, live_analysis=True)
            elif state == 3: # Stop live analysis.
                if self.was_playing == False:
                    self.manipulator.superscan.stop_playing()
                self.stop_live_analysis_event.set()
                self.manipulator.superscan.stop_playing()
            start_stop_analysis_button_next_state(state)

        # Change the functionality and appearance of the button.
        def start_stop_analysis_button_next_state(state = None):
            if state == None: # init
                if "live" in self.image_source_combo_box.current_item.lower() and self.live_analysis_check_box.checked:
                    state = 2 # Start live analysis.
                else:
                    state = 0 # Start single image analysis.
            elif state == 2:
                state = 3 # Stop live analysis.
            elif state == 3:
                state = 2 # Start live analysis.
            self.start_stop_analysis_button.state = state
            self.start_stop_analysis_button.text = button_state_text[state]

        ## GUI elements.

        # Scale calibration rows.
        scale_calibration_row, self.scale_calibration_combo_box = combo_box_template(
            self.ui, 'Mode', ['Before run time (Real space with NN)', 'At run time (Fourier space)']
        )
        scale_calibration_button_container = self.ui.create_row_widget()
        self.scale_calibration_combo_box.on_current_item_changed = scale_calibration_mode_changed
        scale_calibration_row.add(scale_calibration_button_container)
        scale_calibration_row.add_stretch()

        scale_calibration_display_row = self.ui.create_row_widget()
        scale_calibration_display_row.add(self.ui.create_label_widget(_('Sampling: ')))
        self.sampling_label = self.ui.create_label_widget(_('N/A'))
        scale_calibration_display_row.add(self.sampling_label)
        scale_calibration_display_row.add_stretch()
        
        # Field of view (FOV) row.
        fov_display_row = self.ui.create_row_widget()
        fov_display_row.add(self.ui.create_label_widget(_('FOV: ')))
        self.fov_label = self.ui.create_label_widget(_('N/A'))
        fov_display_row.add(self.fov_label)
        fov_display_row.add_stretch()
        
        # Image source row.
        image_source_row, self.image_source_combo_box = combo_box_template(self.ui, 'Image source',
                                                ['MAADF Live', 'HAADF Live', 'Selected data item'])       
        self.image_source_combo_box.on_current_item_changed = image_source_changed

        # Visualize atoms row.
        visualize_atoms_row, self.visualize_atoms_check_box = check_box_template(self.ui, _('Plot atom positions'))
        self.visualize_atoms_check_box.on_checked_changed = visualize_atoms_changed
        
        # Foreign atom detection row.
        auto_detect_foreign_atoms_row, self.auto_detect_foreign_atoms_check_box = check_box_template(self.ui, 'Auto-detect foreign atoms')
        self.auto_detect_foreign_atoms_check_box.on_checked_changed = auto_detect_foreign_atoms_changed

        # Elemental identification rows.
        elemental_id_row1, self.elemental_id_int_radius_line_edit = \
            line_edit_template(self.ui, "Elemental ident.: Int. radius [A]: ")
        self.elemental_id_int_radius_line_edit.on_editing_finished = elemental_id_int_radius_changed
        
        elemental_id_row2, self.elemental_id_exponent_line_edit = \
            line_edit_template(self.ui, "Elemental ident.: Z-exponent: ")
        self.elemental_id_exponent_line_edit.on_editing_finished = elemental_id_exponent_changed
        
        # Control and start/stop analysis rows.
        live_analysis_row = self.ui.create_row_widget()        
        button_state_text = [
            'Start single image analysis',
            'Processing',
            'Start live analysis',
            'Stop live analysis']
        start_stop_analysis_row, self.start_stop_analysis_button = push_button_template(self.ui, 'init')
        self.start_stop_analysis_button.state = None                    
        self.start_stop_analysis_button.on_clicked = start_stop_analysis
        
        # Number of atoms row.
        N_atoms_row = self.ui.create_row_widget()
        self.N_atoms_label = self.ui.create_label_widget('0')
        N_atoms_row.add(self.N_atoms_label)
        N_atoms_row.add(self.ui.create_label_widget(_(' atoms detected')))
        N_atoms_row.add_stretch()

        # Set defaults.
        self.auto_detect_foreign_atoms_check_box.checked = defaults['auto_detect_foreign_atoms']
        auto_detect_foreign_atoms_changed(self.auto_detect_foreign_atoms_check_box.checked)
        self.visualize_atoms_check_box.checked = defaults['visualize_atoms']
        visualize_atoms_changed(self.visualize_atoms_check_box.checked)
        elemental_id_int_radius_changed(str(defaults['elemental_identification_integration_radius_A']))
        elemental_id_exponent_changed(str(defaults['elemental_identification_exponent']))
        image_source_changed(defaults['image_source'])
        scale_calibration_mode_changed(defaults['scale_calibration_mode'])

        # Assemble GUI elements.
        section1.column.add(scale_calibration_row)
        section1.column.add(scale_calibration_display_row)
        section1.column.add(fov_display_row)

        section2.column.add(image_source_row)
        section2.column.add(visualize_atoms_row)
        section2.column.add(auto_detect_foreign_atoms_row)
        section2.column.add(elemental_id_row1)
        section2.column.add(elemental_id_row2)
        section2.column.add_spacing(5)
        section2.column.add(live_analysis_row)
        section2.column.add(start_stop_analysis_row)
        section2.column.add(N_atoms_row)