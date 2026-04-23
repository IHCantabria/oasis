"""App settings — persist user preferences with QSettings."""
from __future__ import annotations
from PyQt5.QtCore import QSettings


class AppSettings:
    _ORG = "OASIS"
    _APP = "GUI"

    def __init__(self):
        self._qs = QSettings(self._ORG, self._APP)

    @property
    def oasis_exe(self) -> str:
        return self._qs.value("oasis_exe", "")

    @oasis_exe.setter
    def oasis_exe(self, val: str) -> None:
        self._qs.setValue("oasis_exe", val)

    @property
    def recent_files(self) -> list:
        val = self._qs.value("recent_files", [])
        return list(val) if val else []

    @recent_files.setter
    def recent_files(self, val: list) -> None:
        self._qs.setValue("recent_files", val)

    def add_recent(self, path: str, max_items: int = 10) -> None:
        recent = self.recent_files
        if path in recent:
            recent.remove(path)
        recent.insert(0, path)
        self.recent_files = recent[:max_items]

    @property
    def last_dir(self) -> str:
        return self._qs.value("last_dir", "")

    @last_dir.setter
    def last_dir(self, val: str) -> None:
        self._qs.setValue("last_dir", val)
