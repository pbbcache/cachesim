#!/usr/bin/env python

from distutils.core import setup, Extension

schedctl_module = Extension('_schedctl',sources = ['schedctl_wrap.cxx', 'schedctl.cpp'])

setup (name = 'schedctl',
       version = '0.1',
       author      = "Juan Carlos Saez",
       description = """Module to use schedctl subsystem from Python""",
       ext_modules = [schedctl_module],
       py_modules = ["schedctl"],
       )
