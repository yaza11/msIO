# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_data_files
from PyInstaller.utils.hooks import collect_all
import os

# paths relative to this spec file, so the build works on any machine
# main.py looks for "about" next to its own __file__, which for the entry script is the bundle root
datas = [(os.path.join(SPECPATH, 'msIO', 'gui', 'about'), '.')]
binaries = []
# older IsoSpecPy builds ship a separate DLL in site-packages/bin, newer ones a .pyd that PyInstaller finds itself
_isospec_dll = os.path.join(SPECPATH, '.venv', 'Lib', 'site-packages', 'bin', 'IsoSpecCppPy.dll')
if os.path.isfile(_isospec_dll):
    binaries.append((_isospec_dll, 'bin'))
hiddenimports = []
datas += collect_data_files('LipidCalculator')
tmp_ret = collect_all('matchms')
datas += tmp_ret[0]; binaries += tmp_ret[1]; hiddenimports += tmp_ret[2]
tmp_ret = collect_all('IsoSpecPy')
datas += tmp_ret[0]; binaries += tmp_ret[1]; hiddenimports += tmp_ret[2]


a = Analysis(
    ['msIO/gui/main.py'],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
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
    name='LipidAnnotations',
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
    name='LipidAnnotations',
)
