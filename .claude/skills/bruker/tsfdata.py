"""Minimal ctypes reader of Bruker TSF files (analysis.tsf / analysis.tsf_bin:
timsTOF / timsOmni spectra without ion mobility separation) through Bruker's
timsdata library (TDF-SDK).

The library is found in this order: the environment variable TIMSDATA_LIB
(a path to libtimsdata.so / timsdata.dll), the PyPI package
opentims-bruker-bridge (pip install opentims-bruker-bridge; it ships the Linux
and Windows builds), libtimsdata.so next to this file."""
import ctypes as C
import os
import sqlite3

import numpy as np


def find_lib():
    path = os.environ.get("TIMSDATA_LIB")
    if path:
        return path
    try:
        import opentims_bruker_bridge
        return opentims_bruker_bridge.get_appropriate_so_path()
    except ImportError:
        pass
    local = os.path.join(os.path.dirname(os.path.realpath(__file__)), "libtimsdata.so")
    if os.path.exists(local):
        return local
    raise RuntimeError("Bruker's timsdata library not found: pip install opentims-bruker-bridge "
                       "or set TIMSDATA_LIB")


class TsfData:
    def __init__(self, d_path, recalibrated=False, lib=None):
        self.dll = C.cdll.LoadLibrary(os.path.realpath(lib or find_lib()))
        dll = self.dll
        dll.tsf_open.argtypes = [C.c_char_p, C.c_uint32]
        dll.tsf_open.restype = C.c_uint64
        dll.tsf_close.argtypes = [C.c_uint64]
        dll.tsf_get_last_error_string.argtypes = [C.c_char_p, C.c_uint32]
        dll.tsf_get_last_error_string.restype = C.c_uint32
        dll.tsf_read_line_spectrum_v2.argtypes = [C.c_uint64, C.c_int64, C.POINTER(C.c_double),
                                                  C.POINTER(C.c_float), C.c_int32]
        dll.tsf_read_line_spectrum_v2.restype = C.c_int32
        dll.tsf_read_profile_spectrum_v2.argtypes = [C.c_uint64, C.c_int64, C.POINTER(C.c_uint32),
                                                     C.c_int32]
        dll.tsf_read_profile_spectrum_v2.restype = C.c_int32
        dll.tsf_index_to_mz.argtypes = [C.c_uint64, C.c_int64, C.POINTER(C.c_double),
                                        C.POINTER(C.c_double), C.c_uint32]
        dll.tsf_index_to_mz.restype = C.c_uint32
        self.handle = dll.tsf_open(d_path.encode(), 1 if recalibrated else 0)
        if self.handle == 0:
            raise RuntimeError(self.error())
        self.conn = sqlite3.connect(os.path.join(d_path, "analysis.tsf"))

    def error(self):
        buf = C.create_string_buffer(512)
        self.dll.tsf_get_last_error_string(buf, 512)
        return buf.value.decode()

    def close(self):
        self.dll.tsf_close(self.handle)
        self.conn.close()

    def index_to_mz(self, frame, idx):
        idx = np.ascontiguousarray(idx, dtype=np.float64)
        mz = np.empty_like(idx)
        ok = self.dll.tsf_index_to_mz(self.handle, frame,
                                      idx.ctypes.data_as(C.POINTER(C.c_double)),
                                      mz.ctypes.data_as(C.POINTER(C.c_double)), len(idx))
        if ok == 0:
            raise RuntimeError(self.error())
        return mz

    def line_spectrum(self, frame):
        """Bruker's centroided (line) spectrum: m/z, intensity."""
        n = 4096
        while True:
            idx = np.empty(n, dtype=np.float64)
            inten = np.empty(n, dtype=np.float32)
            got = self.dll.tsf_read_line_spectrum_v2(self.handle, frame,
                                                     idx.ctypes.data_as(C.POINTER(C.c_double)),
                                                     inten.ctypes.data_as(C.POINTER(C.c_float)), n)
            if got < 0:
                raise RuntimeError(self.error())
            if got <= n:
                return self.index_to_mz(frame, idx[:got]), inten[:got].astype(np.float64)
            n = got

    def profile_spectrum(self, frame):
        """The profile spectrum on the TOF index grid: intensity per index (uint32)."""
        n = 1 << 21
        while True:
            prof = np.empty(n, dtype=np.uint32)
            got = self.dll.tsf_read_profile_spectrum_v2(self.handle, frame,
                                                        prof.ctypes.data_as(C.POINTER(C.c_uint32)), n)
            if got < 0:
                raise RuntimeError(self.error())
            if got <= n:
                return prof[:got]
            n = got
