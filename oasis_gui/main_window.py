"""Main window."""
from __future__ import annotations
import os
import shutil
import zipfile
from typing import Optional

from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QSplitter,
    QListWidget, QStackedWidget, QAction, QToolBar, QFileDialog,
    QMessageBox, QLabel, QStatusBar,
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon

from app_settings import AppSettings
from data.case_data import CaseData
from data.yaml_io import save_yaml, load_yaml

from editors.problem_editor import ProblemEditor
from editors.seafloor_editor import SeaFloorEditor
from editors.waves_editor import WavesEditor
from editors.bcps_editor import BCPsEditor
from editors.bodies_editor import BodiesEditor
from editors.lines_editor import LinesEditor
from editors.springs_editor import SpringsEditor
from editors.morison_editor import MorisonEditor
from editors.winches_editor import WinchesEditor
from editors.owc_editor import OWCEditor
from editors.sinking_editor import SinkingEditor
from editors.turbines_editor import TurbinesEditor

from views.schematic_view import SchematicView
from views.run_view import RunView
from views.results_view import ResultsView
from views.viewer_3d import Viewer3D

# Module list in display order
_MODULES = [
    ("Problem",   ProblemEditor),
    ("Sea Floor", SeaFloorEditor),
    ("Waves",     WavesEditor),
    ("BCPs",      BCPsEditor),
    ("Bodies",    BodiesEditor),
    ("Lines",     LinesEditor),
    ("Springs",   SpringsEditor),
    ("Morison",   MorisonEditor),
    ("Winches",   WinchesEditor),
    ("OWC",       OWCEditor),
    ("Sinking",   SinkingEditor),
    ("Turbines",  TurbinesEditor),
]


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self._settings = AppSettings()
        self._case: CaseData = CaseData()
        self._dirty = False
        self.setWindowTitle("OASIS GUI")
        self.resize(1400, 900)

        self._build_ui()
        self._build_menus()
        self._load_case_into_ui()

    # ── Build UI ──────────────────────────────────────────────────────────────
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_h = QHBoxLayout(central)
        main_h.setContentsMargins(4, 4, 4, 4)

        # Left: module list
        self._mod_list = QListWidget()
        self._mod_list.setMaximumWidth(130)
        for name, _ in _MODULES:
            self._mod_list.addItem(name)
        self._mod_list.currentRowChanged.connect(self._on_module_changed)

        # Centre: editors stacked
        self._editor_stack = QStackedWidget()
        self._editors: list = []
        for name, EditorClass in _MODULES:
            if name == "Winches":
                ed = EditorClass()
            else:
                ed = EditorClass()
            self._editors.append(ed)
            self._editor_stack.addWidget(ed)

        # Right: views stacked
        self._view_stack = QStackedWidget()
        self._schematic = SchematicView()
        self._run_view = RunView(self._settings)
        self._results_view = ResultsView()
        self._viewer_3d = Viewer3D()
        for view in [self._schematic, self._run_view, self._results_view, self._viewer_3d]:
            self._view_stack.addWidget(view)

        main_h.addWidget(self._mod_list)
        main_h.addWidget(self._editor_stack, stretch=3)
        main_h.addWidget(self._view_stack, stretch=4)

        self.setStatusBar(QStatusBar())
        self._status = self.statusBar()

    def _build_menus(self):
        tb = QToolBar("Main")
        self.addToolBar(tb)

        act_new = QAction("New", self)
        act_new.setShortcut("Ctrl+N")
        act_new.triggered.connect(self._new_case)
        tb.addAction(act_new)

        act_open = QAction("Open…", self)
        act_open.setShortcut("Ctrl+O")
        act_open.triggered.connect(self._open_case)
        tb.addAction(act_open)

        act_save = QAction("Save", self)
        act_save.setShortcut("Ctrl+S")
        act_save.triggered.connect(self._save_case)
        tb.addAction(act_save)

        act_saveas = QAction("Save As…", self)
        act_saveas.triggered.connect(self._save_case_as)
        tb.addAction(act_saveas)

        act_export = QAction("Export ZIP…", self)
        act_export.triggered.connect(self._export_zip)
        tb.addAction(act_export)

        tb.addSeparator()

        act_sch = QAction("Schematic", self)
        act_sch.triggered.connect(lambda: self._show_view(0))
        tb.addAction(act_sch)

        act_run = QAction("Run", self)
        act_run.triggered.connect(lambda: self._show_view(1))
        tb.addAction(act_run)

        act_res = QAction("Results", self)
        act_res.triggered.connect(lambda: self._show_view(2))
        tb.addAction(act_res)

        act_3d = QAction("3D View", self)
        act_3d.triggered.connect(lambda: self._show_view(3))
        tb.addAction(act_3d)

        # Menu bar
        menu = self.menuBar()
        file_menu = menu.addMenu("&File")
        file_menu.addAction(act_new)
        file_menu.addAction(act_open)
        file_menu.addAction(act_save)
        file_menu.addAction(act_saveas)
        file_menu.addSeparator()
        file_menu.addAction(act_export)
        file_menu.addSeparator()
        act_quit = QAction("Quit", self)
        act_quit.triggered.connect(self.close)
        file_menu.addAction(act_quit)

    # ── Case I/O ──────────────────────────────────────────────────────────────
    def _new_case(self):
        if not self._maybe_save():
            return
        self._case = CaseData()
        self._dirty = False
        self._load_case_into_ui()
        self._status.showMessage("New case created.")

    def _open_case(self):
        if not self._maybe_save():
            return
        last = self._settings.last_dir
        path, _ = QFileDialog.getOpenFileName(
            self, "Open YAML case", last,
            "YAML files (*.yaml *.yml);;All Files (*)"
        )
        if not path:
            return
        try:
            self._case = load_yaml(path)
            self._case.project_path = os.path.dirname(os.path.dirname(path))
            self._settings.add_recent(path)
            self._settings.last_dir = os.path.dirname(path)
            self._dirty = False
            self._load_case_into_ui()
            self._status.showMessage(f"Opened: {path}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open file:\n{e}")

    def _save_case(self):
        path = self._case.project_path
        if path:
            yaml_path = os.path.join(path, "input", "dataProblem.yaml")
            self._write_case(yaml_path)
        else:
            self._save_case_as()

    def _save_case_as(self):
        last = self._settings.last_dir
        path, _ = QFileDialog.getSaveFileName(
            self, "Save YAML case", last,
            "YAML files (*.yaml *.yml);;All Files (*)"
        )
        if not path:
            return
        if not path.lower().endswith((".yaml", ".yml")):
            path += ".yaml"
        self._case.project_path = os.path.dirname(os.path.dirname(path))
        self._write_case(path)

    def _write_case(self, yaml_path: str):
        self._collect_all_editors()
        os.makedirs(os.path.dirname(yaml_path), exist_ok=True)
        try:
            save_yaml(self._case, yaml_path)
            self._settings.add_recent(yaml_path)
            self._settings.last_dir = os.path.dirname(yaml_path)
            self._dirty = False
            self._status.showMessage(f"Saved: {yaml_path}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save:\n{e}")

    def _export_zip(self):
        """Export project input directory as a ZIP file for HPC cluster."""
        self._save_case()
        project_path = self._case.project_path
        if not project_path:
            QMessageBox.warning(self, "Export", "Save the case first.")
            return
        input_dir = os.path.join(project_path, "input")
        if not os.path.isdir(input_dir):
            QMessageBox.warning(self, "Export", f"Input directory not found:\n{input_dir}")
            return
        zip_path, _ = QFileDialog.getSaveFileName(
            self, "Export as ZIP", self._settings.last_dir,
            "ZIP archives (*.zip);;All Files (*)"
        )
        if not zip_path:
            return
        if not zip_path.lower().endswith(".zip"):
            zip_path += ".zip"
        try:
            with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
                for root, dirs, files in os.walk(input_dir):
                    for fname in files:
                        fpath = os.path.join(root, fname)
                        arcname = os.path.relpath(fpath, project_path)
                        zf.write(fpath, arcname)
            self._status.showMessage(f"Exported: {zip_path}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Export failed:\n{e}")

    # ── Helpers ───────────────────────────────────────────────────────────────
    def _maybe_save(self) -> bool:
        """Prompt to save dirty state. Returns False if user cancelled."""
        if not self._dirty:
            return True
        reply = QMessageBox.question(
            self, "Unsaved changes",
            "Save current case before continuing?",
            QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel,
        )
        if reply == QMessageBox.Save:
            self._save_case()
            return True
        if reply == QMessageBox.Discard:
            return True
        return False

    def _load_case_into_ui(self):
        for ed in self._editors:
            try:
                ed.load_from_case(self._case)
            except Exception:
                pass
        self._run_view.set_case(self._case)
        self._results_view.set_case(self._case)
        self._viewer_3d.set_case(self._case, self._case.project_path)
        self._schematic.update_case(self._case)
        if self._mod_list.count() > 0:
            self._mod_list.setCurrentRow(0)

    def _collect_all_editors(self):
        for ed in self._editors:
            try:
                ed.save_to_case(self._case)
            except Exception:
                pass

    def _on_module_changed(self, row: int):
        if 0 <= row < len(self._editors):
            self._editor_stack.setCurrentIndex(row)

    def _show_view(self, idx: int):
        self._view_stack.setCurrentIndex(idx)
        if idx == 0:
            self._collect_all_editors()
            self._schematic.update_case(self._case)
        elif idx == 2:
            self._results_view.set_case(self._case)

    def closeEvent(self, event):
        if self._maybe_save():
            event.accept()
        else:
            event.ignore()
