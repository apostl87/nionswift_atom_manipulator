"""
Widgets library.
"""
#########################################################################################################
##  Special thanks for providing the code base for this library go to Jacob Madsen (jacobjma@github).  ##
#########################################################################################################

# Nion libraries
from nion.ui import Widgets


def line_edit_template(ui, label, default_text=None, placeholder_text=None):
    row = ui.create_row_widget()
    row.add(ui.create_label_widget(label))
    row.add_spacing(5)
    widget = ui.create_line_edit_widget()
    row.add(widget)
    row.add_spacing(5)
    widget.text = default_text
    return row, widget


def push_button_template(ui, label, callback=None):
    row = ui.create_row_widget()
    widget = ui.create_push_button_widget(label)
    widget.on_clicked = callback
    row.add(widget)
    return row, widget


def combo_box_template(ui, label, items, indent=False):
    row = ui.create_row_widget()
    if indent:
        row.add_spacing(8)
    row.add(ui.create_label_widget(label))
    row.add_spacing(5)
    widget = ui.create_combo_box_widget(items=items)
    row.add(widget)
    row.add_stretch()
    return row, widget


def check_box_template(ui, label):
    row = ui.create_row_widget()
    widget = ui.create_check_box_widget(label)
    row.add(widget)
    return row, widget


class ScrollArea:

    def __init__(self, ui):
        self.__ui = ui
        self.__scroll_area_widget = ui.create_scroll_area_widget()

    @property
    def _ui(self):
        return self.__ui

    @property
    def _widget(self):
        return self.__scroll_area_widget

    @property
    def content(self):
        return self._widget.content

    @content.setter
    def content(self, value):
        self._widget.content = value


class Section:

    def __init__(self, ui, title):
        self.__ui = ui
        self.__section_content_column = self.__ui._ui.create_column_widget()
        self.__section_widget = Widgets.SectionWidget(self.__ui._ui, title, self.__section_content_column, 'test')
        self.column = ui.create_column_widget()
        self.__section_content_column.add(self.column._widget)

    @property
    def _ui(self):
        return self.__ui

    @property
    def _widget(self):
        return self.__section_widget

    @property
    def _section_content_column(self):
        return self.__section_content_column
