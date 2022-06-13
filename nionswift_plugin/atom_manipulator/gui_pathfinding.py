import gettext
import threading
import numpy as np

# Custom libraries
from . import lib_pathfinding
from .lib_utils import AtomManipulatorModule
from .lib_widgets import Section, line_edit_template, check_box_template, combo_box_template, push_button_template

_ = gettext.gettext

# Defaults on initialization.
defaults = {'max_bond_length': 1.8} # in Angstroem


class PathFindingModule(AtomManipulatorModule):
    
    def __init__(self, ui, api, document_controller, manipulator):
        super().__init__(ui, api, document_controller)
        self.manipulator = manipulator # AtomManipulatorDelegate object
        self.max_bond_length = None # Internally Nion Swift calculates in nm.
        
        # Events.
        self.rdy = threading.Event()

    # GUI creation method. 
    def create_widgets(self, column):
        self.section = Section(self.ui, 'Pathfinding')
        column.add(self.section)

        # Callback functions.
        def add_foreign_atoms_clicked():
            add_remove_atoms_sites(0)   
        def remove_foreign_atoms_clicked():
            add_remove_atoms_sites(1)
        def add_target_sites_clicked():
            add_remove_atoms_sites(2)
        def remove_target_sites_clicked():
            add_remove_atoms_sites(3)

        #### GUI elements.
        ## Adding, removing, and displaying foreign atoms and target sites.
        # Toggle text of buttons.
        def toggle_button1_text(): # Foreign atoms
            while True:
                yield _('Add')
                yield _('Stop adding')
        def toggle_button2_text(): # Foreign atoms
            while True:
                yield _('Remove')
                yield _('Stop removing')
        def toggle_button3_text(): # Target sites
            while True:
                yield _('Add')
                yield _('Stop adding')
        def toggle_button4_text(): # Target sites
            while True:
                yield _('Remove')
                yield _('Stop removing')
        buttons_text = [toggle_button1_text(), toggle_button2_text(), toggle_button3_text(), toggle_button4_text()]
 
        # Button widgets.
        self.add_foreign_atoms_button = self.ui.create_push_button_widget(next(buttons_text[0]))
        self.add_foreign_atoms_button.state = False
        self.add_foreign_atoms_button.on_clicked = add_foreign_atoms_clicked
        self.remove_foreign_atoms_button = self.ui.create_push_button_widget(next(buttons_text[1]))
        self.remove_foreign_atoms_button.state = False
        self.remove_foreign_atoms_button.on_clicked = remove_foreign_atoms_clicked
        self.add_target_sites_button = self.ui.create_push_button_widget(next(buttons_text[2]))
        self.add_target_sites_button.state = False
        self.add_target_sites_button.on_clicked = add_target_sites_clicked
        self.remove_target_sites_button = self.ui.create_push_button_widget(next(buttons_text[3]))
        self.remove_target_sites_button.state = False
        self.remove_target_sites_button.on_clicked = remove_target_sites_clicked
        
        self.add_remove_buttons = []
        self.add_remove_buttons.append(self.add_foreign_atoms_button)
        self.add_remove_buttons.append(self.remove_foreign_atoms_button)
        self.add_remove_buttons.append(self.add_target_sites_button)
        self.add_remove_buttons.append(self.remove_target_sites_button)

        # Label widgets.
        self.N_foreign_atoms_label = self.ui.create_label_widget('0')
        self.N_target_sites_label = self.ui.create_label_widget('0')

        def add_remove_atoms_sites(button_idx):
            last_mode = np.where(list(map(lambda x: x.state, self.add_remove_buttons)))[0]
            if last_mode.size > 0 and button_idx != last_mode[0]:
                self.add_remove_buttons[last_mode[0]].text = next(buttons_text[last_mode[0]])
                self.add_remove_buttons[last_mode[0]].state = False
            self.add_remove_buttons[button_idx].text = next(buttons_text[button_idx])
            self.add_remove_buttons[button_idx].state = not self.add_remove_buttons[button_idx].state
            lib_pathfinding.add_or_remove_foreign_atoms_or_target_sites(self.manipulator,
                                                                        mode=button_idx,
                                                                        startstop=self.add_remove_buttons[button_idx].state)

        # Row widgets.
        foreign_atoms_row = self.ui.create_row_widget()
        foreign_atoms_row.add(self.ui.create_label_widget(_('Foreign atoms: ')))
        foreign_atoms_row.add(self.N_foreign_atoms_label)
        foreign_atoms_row.add(self.add_foreign_atoms_button)
        foreign_atoms_row.add(self.remove_foreign_atoms_button)
        foreign_atoms_row.add_stretch()
        
        target_sites_row = self.ui.create_row_widget()
        target_sites_row.add(self.ui.create_label_widget(_('Target sites: ')))
        target_sites_row.add(self.N_target_sites_label)
        target_sites_row.add(self.add_target_sites_button)
        target_sites_row.add(self.remove_target_sites_button)
        target_sites_row.add_stretch()
        
        # Other widget rows.
        max_bond_length_row, self.max_bond_length_line_edit = line_edit_template(self.ui, 'Max. bond length [A]')
        def max_bond_length_editing_finished(text):
            if len(text) > 0:
                try:
                    self.max_bond_length = float(text)
                except:
                    pass
                finally:
                    self.max_bond_length_line_edit.text = f"{self.max_bond_length:.2f}"
        self.max_bond_length_line_edit.on_editing_finished = max_bond_length_editing_finished
        
        find_paths_row, self.find_paths_button = push_button_template(self.ui, 'Find paths')
        def find_paths_clicked():
            lib_pathfinding.find_paths(self.manipulator)
        self.find_paths_button.on_clicked = find_paths_clicked
        self.move_probe_button = self.ui.create_push_button_widget('Move probe') # to the next demanded position')
        def move_probe_clicked():
            lib_pathfinding.move_probe(self.manipulator)
        self.move_probe_button.on_clicked = move_probe_clicked
        find_paths_row.add(self.move_probe_button)
        
        # Set defaults.
        max_bond_length_editing_finished(str(defaults['max_bond_length']))
        
        # Assemble GUI elements.
        self.section.column.add(foreign_atoms_row)
        self.section.column.add(target_sites_row)
        self.section.column.add(max_bond_length_row)
        self.section.column.add(find_paths_row)