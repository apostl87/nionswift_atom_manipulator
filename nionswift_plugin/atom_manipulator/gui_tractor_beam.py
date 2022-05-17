import gettext
import threading

import logging

# Custom libraries
from . import lib_utils
from . import lib_structure_recognition
from . import lib_pathfinding
from .lib_utils import AtomManipulatorModule
from .lib_widgets import Section, line_edit_template, check_box_template, combo_box_template
from adf_feedback import adf_feedback as adffb

_ = gettext.gettext

# Parameter presets.
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

# Defaults on initialization.
defaults = {'offline_test_mode': False,
            'auto_detect_changes': True,
            'preset': '60 kV'}

 
class TractorBeamModule(AtomManipulatorModule):

    def __init__(self, ui, api, document_controller, manipulator=None):
        super().__init__(ui, api, document_controller)
        self.manipulator = manipulator # AtomManipulatorDelegate object
        self.frame_timeout = None
        self.reposition_timeout = None
        self.jump_threshold = None
        self.drift_threshold = None
        
        # Events.
        self.rdy = threading.Event()

    # GUI creation method.
    def create_widgets(self, column):
        self.section = Section(self.ui, 'TractorBeam')
        column.add(self.section)

        # Callback functions.
        def onbutton_clicked():
            logging.info(lib_utils.log_message("Starting TractorBeam"))
            try:
                # In this (manual) operation mode, ADFFeedback needs to be stopped by a direct call of the method stopmap.
                self.ADFFeedback = adffb.ADFFeedbackDelegate(
                    self.api,
                    offline_test_mode = self.otm_check_box.checked,
                    autodetect_changes = self.adc_check_box.checked
                )

                self.ADFFeedback.startmap(
                    self.frame_timeout,
                    self.jump_threshold,
                    self.drift_threshold,
                    self.reposition_timeout
                )

            except Exception as exc:
                logging.info(lib_utils.log_message("Exception when calling adf_feedback"))
                print(exc)

        def offbutton_clicked():
            logging.info(lib_utils.log_message("Stopping TractorBeam"))
            try:
                self.ADFFeedback.stopmap()
            except:
                pass

        def frame_timeout_editing_finished(text):
            if len(text) > 0:
                try:
                    self.frame_timeout = round(float(text))
                    #logging.info(lib_utils.log_message(f"frame_timeout: {self.frame_timeout}"))
                except: pass
                finally: self.frame_timeout_line_edit.text = f"{self.frame_timeout:d}"

        def reposition_timeout_editing_finished(text):
            if len(text) > 0:
                try:
                    self.reposition_timeout = round(float(text), 1)
                    #logging.info(lib_utils.log_message(f"reposition_timeout: {self.reposition_timeout}"))
                finally:
                    self.reposition_timeout_line_edit.text = f"{self.reposition_timeout:.1f}"

        def jump_threshold_editing_finished(text):
            if len(text) > 0:
                try:
                    self.jump_threshold = round(float(text), 2)
                    #logging.info(lib_utils.log_message(f"jump_threshold: {self.jump_threshold}"))
                finally:
                    self.jump_threshold_line_edit.text = f"{self.jump_threshold:.2f}"

        def drift_threshold_editing_finished(text):
            if len(text) > 0:
                try:
                    self.drift_threshold = round(float(text), 2)
                    #logging.info(lib_utils.log_message(f"drift_threshold: {self.drift_threshold}"))
                finally:
                    self.drift_threshold_line_edit.text = f"{self.drift_threshold:.2f}"

        def set_preset(name):
            frame_timeout_editing_finished(presets[name]['frame_timeout'])
            reposition_timeout_editing_finished(presets[name]['reposition_timeout'])
            jump_threshold_editing_finished(presets[name]['jump_threshold'])
            drift_threshold_editing_finished(presets[name]['drift_threshold'])

        ## GUI elements.

        # Offline test mode row.
        otm_row, self.otm_check_box = check_box_template(self.ui, 'Offline test mode')
        self.otm_check_box.checked = defaults['offline_test_mode']
        
        # Auto detect changes row.
        adc_row, self.adc_check_box = check_box_template(self.ui, 'Auto-detect changes')
        self.adc_check_box.checked = defaults['auto_detect_changes']
        
        # Frame timeout row.
        frame_timeout_row, self.frame_timeout_line_edit = line_edit_template(self.ui, 'Tractor time [s]')
        self.frame_timeout_line_edit.on_editing_finished = frame_timeout_editing_finished

        # Probe reposition timeout row.
        reposition_timeout_row, self.reposition_timeout_line_edit = line_edit_template(self.ui, 'Reposition timeout [s]')
        self.reposition_timeout_line_edit.on_editing_finished = reposition_timeout_editing_finished

        # Jump threshold row.
        jump_threshold_row, self.jump_threshold_line_edit = line_edit_template(self.ui, 'Auto-detect jump threshold [%]')
        self.jump_threshold_line_edit.on_editing_finished = jump_threshold_editing_finished

        # Drift threshold row.
        drift_threshold_row, self.drift_threshold_line_edit = line_edit_template(self.ui, 'Auto-detect drift threshold [%]')
        self.drift_threshold_line_edit.on_editing_finished = drift_threshold_editing_finished
        
        # Buttons row.
        onbutton = self.ui.create_push_button_widget(_("Start"))
        onbutton.on_clicked = onbutton_clicked
        
        offbutton = self.ui.create_push_button_widget(_("Stop"))
        offbutton.on_clicked = offbutton_clicked

        startstop_button_row = self.ui.create_row_widget()
        startstop_button_row.add(onbutton)
        startstop_button_row.add(offbutton)
        
        # Preset row.
        preset_row, self.preset_combo_box = combo_box_template(
            self.ui, 'Preset', ['55 kV', '60 kV'], indent=False)
        self.preset_combo_box.on_current_item_changed = set_preset
        
        # Set defaults.
        self.preset_combo_box.current_item = defaults['preset']

        # Assemble GUI elements.
        self.section.column.add(otm_row)
        self.section.column.add(preset_row)
        self.section.column.add(adc_row)
        self.section.column.add(frame_timeout_row)
        self.section.column.add(reposition_timeout_row)
        self.section.column.add(jump_threshold_row)
        self.section.column.add(drift_threshold_row)
        self.section.column.add(startstop_button_row)