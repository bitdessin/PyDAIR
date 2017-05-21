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
rm ./PyDAIR/.DS_Store
rm ./PyDAIR/bin/*.pyc
rm ./PyDAIR/bin/.DS_Store
rm ./PyDAIR/app/*.pyc
rm ./PyDAIR/app/.DS_Store
rm ./PyDAIR/io/*.pyc
rm ./PyDAIR/io/.DS_Store
rm ./PyDAIR/plot/*.pyc
rm ./PyDAIR/plot/.DS_Store
rm ./PyDAIR/seq/*.pyc
rm ./PyDAIR/seq/.DS_Store
rm ./PyDAIR/stats/*.pyc
rm ./PyDAIR/stats/.DS_Store
rm ./PyDAIR/templates/*.pyc
rm ./PyDAIR/templates/.DS_Store
rm ./PyDAIR/test/*.pyc
rm ./PyDAIR/test/.DS_Store
rm ./PyDAIR/utils/*.pyc
rm ./PyDAIR/utils/.DS_Store
rm ./PyDAIR/sim/*.pyc
rm ./PyDAIR/sim/.DS_Store

rm -rf ./PyDAIR/__pycache__
rm -rf ./PyDAIR/app/__pycache__
rm -rf ./PyDAIR/bin/__pycache__
rm -rf ./PyDAIR/io/__pycache__
rm -rf ./PyDAIR/plot/__pycache__
rm -rf ./PyDAIR/seq/__pycache__
rm -rf ./PyDAIR/sim/__pycache__
rm -rf ./PyDAIR/stats/__pycache__
rm -rf ./PyDAIR/test/__pycache__
rm -rf ./PyDAIR/utils/__pycache__


## PyPI registration commands
# python setup.py sdist
# python setup.py bdist_wheel --universal
# twine upload dist/*

