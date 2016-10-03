#!/bin/bash


## Installation comamnds
pip3 uninstall pydair -y
pip3 install --no-deps -e .

pip uninstall pydair -y
pip install --no-deps -e .



## Clean up
#rm -rf PyDAIR.egg-info
#rm -rf build
rm ./PyDAIR/test/data/test_output_*
rm ./PyDAIR/test/data/results/*
rm ./PyDAIR/*.pyc
rm ./PyDAIR/bin/*.pyc
rm ./PyDAIR/app/*.pyc
rm ./PyDAIR/io/*.pyc
rm ./PyDAIR/plot/*.pyc
rm ./PyDAIR/seq/*.pyc
rm ./PyDAIR/stats/*.pyc
rm ./PyDAIR/templates/*.pyc
rm ./PyDAIR/test/*.pyc
rm ./PyDAIR/utils/*.pyc



## PyPI registration commands
# python setup.py sdist
# python setup.py bdist_wheel --universal
# twine upload dist/*

