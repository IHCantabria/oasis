"""OASIS GUI entry point."""
import sys
import os
from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
from main_window import MainWindow


def _load_stylesheet(app: QApplication) -> None:
    qss_path = os.path.join(os.path.dirname(__file__), "style.qss")
    if os.path.isfile(qss_path):
        with open(qss_path, "r", encoding="utf-8") as fh:
            app.setStyleSheet(fh.read())


def _set_icon(app: QApplication) -> None:
    # Try the docs/assets logo first, then docs/theory-manual logo
    candidates = [
        os.path.join(os.path.dirname(__file__), "..", "docs", "assets", "logo.png"),
        os.path.join(os.path.dirname(__file__), "..", "docs", "theory-manual", "logo.png"),
    ]
    for path in candidates:
        path = os.path.normpath(path)
        if os.path.isfile(path):
            app.setWindowIcon(QIcon(path))
            return


def main():
    # High-DPI support
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

    app = QApplication(sys.argv)
    app.setApplicationName("OASIS GUI")
    app.setOrganizationName("OASIS")
    app.setStyle("Fusion")

    _load_stylesheet(app)
    _set_icon(app)

    win = MainWindow()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
