from distutils.core import setup, Extension


rwmodels_module = Extension('RWModelsCPP',sources = ['RWModels.cpp'])
setup(name = 'RWModelsCPP',version='0.0',description = 'This is a rw phys net module',ext_modules = [rwmodels_module])
