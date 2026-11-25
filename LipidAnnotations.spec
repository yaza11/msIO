# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller build of the msIO GUI (msIO/gui/main.py), used locally and by .github/workflows/build-exe.yml:

    pyinstaller --noconfirm --clean LipidAnnotations.spec

The build goes to dist/msIO_GUI/. Run it with the Python environment that has msIO's dependencies installed
(including LipidCalculator and IsoSpecPy), the libraries are collected from there.
"""
import glob
import os
import sysconfig

from PyInstaller.utils.hooks import collect_all

NAME = 'msIO_GUI'

# main.py looks for the guide "about" next to its own __file__, which for the entry script is the bundle root
datas = [(os.path.join(SPECPATH, 'msIO', 'gui', 'about'), '.')]
binaries = []
hiddenimports = ['scipy._external.array_api_compat.numpy.fft']

for package in ['matchms', 'LipidCalculator', 'IsoSpecPy']:
    package_datas, package_binaries, package_hiddenimports = collect_all(package)
    datas += package_datas
    binaries += package_binaries
    hiddenimports += package_hiddenimports

# IsoSpecPy loads its C library (IsoSpecCppPy*.dll / .pyd / .so) with cffi at runtime, PyInstaller can't see it.
# IsoSpecPy searches its own folder and the folder above it, so it goes to the bundle root. Depending on the
# IsoSpecPy version it is in site-packages, site-packages/bin or site-packages/IsoSpecPy.
_site_packages = sysconfig.get_paths()['purelib']
_isospec_libs = sorted({
    path
    for pattern in ['IsoSpecCppPy*', os.path.join('*', 'IsoSpecCppPy*')]
    for path in glob.glob(os.path.join(_site_packages, pattern))
    if path.endswith(('.dll', '.pyd', '.so', '.dylib'))
})
if not _isospec_libs:
    raise SystemExit(f'IsoSpecCppPy library not found in {_site_packages}, is IsoSpecPy installed?')
binaries += [(path, '.') for path in _isospec_libs]


a = Analysis(
    [os.path.join(SPECPATH, 'msIO', 'gui', 'main.py')],
    pathex=[SPECPATH],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    # windowed exe has no stdout / stderr, redirect them to a log file before anything else is imported
    runtime_hooks=[os.path.join(SPECPATH, 'msIO', 'gui', '_gui_log.py')],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name=NAME,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name=NAME,
)
