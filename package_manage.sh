#!/bin/bash


## Installation comamnds
pip3 uninstall pydair -y
pip3 install --no-deps .

pip uninstall pydair -y
pip install --no-deps .


## Create API documents
<< '#__C__'
rm -rf docs/api
sphinx-apidoc -f -F -o docs/api/ PyDAIR
cd docs/api/
cat << EOF >> conf.py

def skip(app, what, name, obj, skip, options):
    if name == "__init__":
        return False
    return skip

def setup(app):
    app.connect("autodoc-skip-member", skip)

EOF

make html
cd -
# Change 'html_theme' in docs/api/conf.py to 'sphinx_rtd_theme'.
#__C__

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
rm -rf ./PyDAIR/PyDAIR.egg-info



## PyPI registration commands
# python setup.py sdist
# python setup.py bdist_wheel --universal
# twine upload dist/*

