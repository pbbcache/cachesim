#!/bin/bash
swig -c++ -python schedctl.i   
python3 setup.py build_ext --inplace 
cp _schedctl.cpython-311-x86_64-linux-gnu.so schedctl.py ../../simulator/common


