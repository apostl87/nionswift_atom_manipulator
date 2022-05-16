import threading
import logging
import gettext
from .lib_utils import AtomManipulatorModule
from .lib_widgets import Section, line_edit_template, check_box_template, combo_box_template
from adf_feedback import adf_feedback as adffb
from . import lib_structure_recognition
from . import lib_pathfinding

_ = gettext.gettext

presets = { '55 kV':
               {'frame_timeout': '45',
                'reposition_timeout': '2',
                'jump_threshold': '12',
                'drift_threshold': '4',
                },
            '60 kV':
               {'frame_timeout': '15',
                'reposition_timeout': '2',
                'jump_threshold': '15',
                'drift_threshold': '4',
                }
           }
    
defaults = {'offline_test_mode': False, 'auto_detect_changes': True}
    
class TractorBeamModule(AtomManipulatorModule):

    def __init__(self, ui, api, document_controller, manipulator=None):
        super().__init__(ui, api, document_controller)
        self.frame_timeout = None
        self.reposition_timeout = None
        self.jump_threshold = None
        self.drift_threshold = None
        self.manipulator = manipulator
        self.rdy = threading.Event()

    def create_widgets(self, column):
        self.section = Section(self.ui, 'TractorBeam')
        column.add(self.section)

        otm_row, self.otm_check_box = check_box_template(self.ui, 'Offline test mode')
        self.otm_check_box.checked = defaults['offline_test_mode']
        adc_row, self.adc_check_box = check_box_template(self.ui, 'Auto-detect changes')
        self.adc_check_box.checked = defaults['auto_detect_changes']
        
        # Row edits.
        frame_timeout_row, self.frame_timeout_line_edit = line_edit_template(self.ui, 'Tractor time [s]')
        reposition_timeout_row, self.reposition_timeout_line_edit = line_edit_template(self.ui, 'Reposition timeout [s]')
        jump_threshold_row, self.jump_threshold_line_edit = line_edit_template(self.ui, 'Auto-detect jump threshold [%]')
        drift_threshold_row, self.drift_threshold_line_edit = line_edit_template(self.ui, 'Auto-detect drift threshold [%]')
        
        def frame_timeout_editing_finished(text):
            if len(text) > 0:
                try:
                    self.frame_timeout = round(float(text))
                    #logging.info(lib_utils.log_message(f"frame_timeout: {self.frame_timeout}"))
                except: pass
                finally: self.frame_timeout_line_edit.text = f"{self.frame_timeout:d}"
        self.frame_timeout_line_edit.on_editing_finished = frame_timeout_editing_finished
                    
        def reposition_timeout_editing_finished(text):
            if len(text) > 0:
                try:
                    self.reposition_timeout = round(float(text), 1)
                    #logging.info(lib_utils.log_message(f"reposition_timeout: {self.reposition_timeout}"))
                finally:
                    self.reposition_timeout_line_edit.text = f"{self.reposition_timeout:.1f}"
        self.reposition_timeout_line_edit.on_editing_finished = reposition_timeout_editing_finished            
        
        def jump_threshold_editing_finished(text):
            if len(text) > 0:
                try:
                    self.jump_threshold = round(float(text), 2)
                    #logging.info(lib_utils.log_message(f"jump_threshold: {self.jump_threshold}"))
                finally:
                    self.jump_threshold_line_edit.text = f"{self.jump_threshold:.2f}"
        self.jump_threshold_line_edit.on_editing_finished = jump_threshold_editing_finished            
                    
        def drift_threshold_editing_finished(text):
            if len(text) > 0:
                try:
                    self.drift_threshold = round(float(text), 2)
                    #logging.info(lib_utils.log_message(f"drift_threshold: {self.drift_threshold}"))
                finally:
                    self.drift_threshold_line_edit.text = f"{self.drift_threshold:.2f}"
        self.drift_threshold_line_edit.on_editing_finished = drift_threshold_editing_finished            
        
        # Buttons.
        onbutton = self.ui.create_push_button_widget(_("Start"))
        def onbutton_clicked():
            logging.info(lib_utils.log_message("Starting TractorBeam"))
            try:
                self.ADFFeedback = adffb.ADFFeedbackDelegate(self.api, offline_test_mode = self.otm_check_box.checked, autodetect_changes = self.adc_check_box.checked)        
                self.ADFFeedback.startmap(self.frame_timeout, self.jump_threshold, self.drift_threshold, self.reposition_timeout)
            except Exception as exc:
                logging.info(lib_utils.log_message("Exception when calling adf_feedback"))
                print(exc)

        offbutton = self.ui.create_push_button_widget(_("Stop"))
        def offbutton_clicked():
            logging.info(lib_utils.log_message("Stopping TractorBeam"))
            try:
                self.ADFFeedback.stopmap()
            except:
                pass
        onbutton.on_clicked = onbutton_clicked
        offbutton.on_clicked = offbutton_clicked
        startstop_button_row = self.ui.create_row_widget()
        startstop_button_row.add(onbutton)
        startstop_button_row.add(offbutton)
        
        # Presets.
        preset_row, self.preset_combo_box = combo_box_template(
            self.ui, 'Preset', ['55 kV', '60 kV'], indent=False)
        def set_preset(name):
            frame_timeout_editing_finished(presets[name]['frame_timeout'])
            reposition_timeout_editing_finished(presets[name]['reposition_timeout'])
            jump_threshold_editing_finished(presets[name]['jump_threshold'])
            drift_threshold_editing_finished(presets[name]['drift_threshold'])
        self.preset_combo_box.on_current_item_changed = set_preset
        
        # Defaults.
        self.preset_combo_box.current_item = '60 kV' 
    
        self.section.column.add(otm_row)
        self.section.column.add(preset_row)
        self.section.column.add(adc_row)
        self.section.column.add(frame_timeout_row)
        self.section.column.add(reposition_timeout_row)
        self.section.column.add(jump_threshold_row)
        self.section.column.add(drift_threshold_row)
        self.section.column.add(startstop_button_row)