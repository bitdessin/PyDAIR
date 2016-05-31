from setuptools import setup, find_packages


try:
    import pypandoc
    read_md = lambda f: pypandoc.convert(f, 'rst')
except ImportError:
    print('Warning: pypandoc module not found.')
    read_md = lambda f: open(f, 'r').read()


setup(
    name        = 'PyDAIR',
    version     = '0.1.0',
    description = 'Python library for diversity analysis of immune repertoire.',
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
    ],
    keywords     = 'blast, bioinformatics',
    author       = 'Jianqiang Sun, Xi Fu',
    author_email = 'wukong@bi.a.u-tokyo.ac.jp',
    url          = 'https://github.com/jqsunac/PyDAIR',
    license      = 'GNU',
    packages     = find_packages(exclude=['examples', 'tests']),
    scripts      = ['PyDAIR/bin/pydair-parseseq',
                    'PyDAIR/bin/pydair-analysis'],
    include_package_data = True,
    zip_safe = True,
    long_description = read_md('README.md'),
    install_requires = ['matplotlib', 'pandas', 'biopython'],
)

