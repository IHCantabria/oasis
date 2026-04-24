"""Run view — launch OASIS and stream log output."""
from __future__ import annotations
import os
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QPushButton, QLineEdit, QTextEdit, QProgressBar, QLabel, QComboBox,
    QFileDialog, QMessageBox,
)
from PyQt5.QtCore import Qt, QProcess, QProcessEnvironment


def _find_oasis_exe() -> str:
    """Try to locate OASIS.exe relative to this file's directory."""
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(here, "..", "bin", "OASIS.exe"),
        os.path.join(here, "..", "bin", "oasis.exe"),
        os.path.join(here, "..", "build", "Release", "OASIS.exe"),
        os.path.join(here, "..", "build", "x64", "Release", "OASIS.exe"),
    ]
    for c in candidates:
        c = os.path.normpath(c)
        if os.path.isfile(c):
            return c
    return ""


class RunView(QWidget):
    """Panel to configure and run OASIS, streaming live log."""

    def __init__(self, app_settings, parent=None):
        super().__init__(parent)
        self._settings = app_settings
        self._process: QProcess | None = None
        self._case = None
        self._build_ui()
        # Auto-fill exe path if not yet set
        if not self._settings.oasis_exe:
            auto = _find_oasis_exe()
            if auto:
                self._settings.oasis_exe = auto
                self._exe_path.setText(auto)

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)

        # ── OASIS executable ──────────────────────────────────────────────────
        exe_box = QGroupBox("OASIS executable")
        ef = QFormLayout(exe_box)
        self._exe_path = QLineEdit(self._settings.oasis_exe)
        browse_btn = QPushButton("Browse…")
        browse_btn.clicked.connect(self._browse_exe)
        exe_row = QHBoxLayout()
        exe_row.addWidget(self._exe_path)
        exe_row.addWidget(browse_btn)
        ef.addRow("oasis.exe:", exe_row)
        layout.addWidget(exe_box)

        # ── Run options ────────────────────────────────────────────────────────
        opt_box = QGroupBox("Run options")
        of = QFormLayout(opt_box)
        self._verbosity = QComboBox()
        self._verbosity.addItems(["INFO", "DEBUG", "FAST"])
        of.addRow("Verbosity:", self._verbosity)
        self._log_file = QLineEdit("")
        of.addRow("Log file (optional):", self._log_file)
        layout.addWidget(opt_box)

        # ── Buttons ────────────────────────────────────────────────────────────
        btn_row = QHBoxLayout()
        self._btn_run = QPushButton("▶  Run OASIS")
        self._btn_run.setStyleSheet("font-weight:bold; color:darkgreen;")
        self._btn_stop = QPushButton("■  Stop")
        self._btn_stop.setEnabled(False)
        self._btn_stop.setStyleSheet("color:darkred;")
        btn_row.addWidget(self._btn_run)
        btn_row.addWidget(self._btn_stop)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        # ── Progress bar ───────────────────────────────────────────────────────
        self._progress = QProgressBar()
        self._progress.setRange(0, 100)
        self._progress.setValue(0)
        layout.addWidget(self._progress)

        # ── Log output ─────────────────────────────────────────────────────────
        layout.addWidget(QLabel("Live log output:"))
        self._log = QTextEdit()
        self._log.setReadOnly(True)
        self._log.setFontFamily("Courier New")
        layout.addWidget(self._log)

        self._btn_run.clicked.connect(self._on_run)
        self._btn_stop.clicked.connect(self._on_stop)

    def set_case(self, case) -> None:
        self._case = case

    def _browse_exe(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select OASIS executable", "",
            "Executables (*.exe);;All Files (*)"
        )
        if path:
            self._exe_path.setText(path)
            self._settings.oasis_exe = path

    def _on_run(self):
        exe = self._exe_path.text().strip()
        if not exe or not os.path.isfile(exe):
            QMessageBox.warning(
                self, "OASIS executable not found",
                "Cannot find the OASIS executable.\n\n"
                "Please use the Browse button to locate OASIS.exe."
            )
            return
        if self._case is None or not self._case.project_path:
            self._log.append("<span style='color:red'>ERROR: No project path set. Save the case first.</span>")
            return
        project_dir = self._case.project_path
        args = [project_dir, "--verbosity", self._verbosity.currentText()]
        log_file = self._log_file.text().strip()
        if log_file:
            args += ["--log", log_file]

        self._log.clear()
        self._progress.setValue(0)
        self._btn_run.setEnabled(False)
        self._btn_stop.setEnabled(True)

        self._process = QProcess(self)
        self._process.setProcessEnvironment(QProcessEnvironment.systemEnvironment())
        self._process.readyReadStandardOutput.connect(self._on_stdout)
        self._process.readyReadStandardError.connect(self._on_stderr)
        self._process.finished.connect(self._on_finished)
        self._process.start(exe, args)
        if not self._process.waitForStarted(3000):
            self._log.append("<span style='color:red'>ERROR: Failed to start process.</span>")
            self._btn_run.setEnabled(True)
            self._btn_stop.setEnabled(False)

    def _on_stop(self):
        if self._process and self._process.state() == QProcess.Running:
            self._process.kill()

    def _on_stdout(self):
        data = self._process.readAllStandardOutput().data().decode("utf-8", errors="replace")
        self._append_log(data)
        self._update_progress(data)

    def _on_stderr(self):
        data = self._process.readAllStandardError().data().decode("utf-8", errors="replace")
        self._append_log(f"<span style='color:darkorange'>{_escape(data)}</span>")

    def _append_log(self, text: str) -> None:
        self._log.moveCursor(self._log.textCursor().End)
        self._log.insertPlainText(text)
        self._log.ensureCursorVisible()

    def _update_progress(self, text: str) -> None:
        """Parse lines like '  40.0% completed' and update progress bar."""
        for line in text.splitlines():
            line = line.strip()
            if "%" in line:
                for token in line.split():
                    t = token.replace("%", "")
                    try:
                        val = float(t)
                        if 0.0 <= val <= 100.0:
                            self._progress.setValue(int(val))
                    except ValueError:
                        pass

    def _on_finished(self, exit_code: int, exit_status):
        if exit_code == 0:
            self._progress.setValue(100)
            self._log.append("\n<span style='color:green'>OASIS finished successfully.</span>")
        else:
            self._log.append(
                f"\n<span style='color:red'>OASIS exited with code {exit_code}.</span>"
            )
        self._btn_run.setEnabled(True)
        self._btn_stop.setEnabled(False)
        self._process = None


def _escape(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
