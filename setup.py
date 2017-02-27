from setuptools import setup
#from distutils.core import setup # no longer recommended

import os

setup(
    name='sba',
    version='1.6.7',
    packages=['sba',"sba.test"],
    scripts=['bin/eucsbademo.py'],

    # dependencies
    #include_package_data=True,
    install_requires=[
        #"quaternions (>=0.0)" no longer required
        "numpy>=1.9.2",
    ],
    dependency_links=[
        'http://www.ics.forth.gr/~lourakis/sba', # lourakis C library
        'https://bitbucket.org/devangel77b/libsbaprojs' # projections
    ],
    #zip_safe=False,

    # metadata for upload to PyPI
    author='Dennis Evangelista',
    author_email='devangel77b@gmail.com',
    description="wrapper for Lourakis' sparse bundle adjustment C library",
    license='GNU GPLv3',
    keywords='SBA, sparse bundle adjustment, calibration, camera, camera calibration, photogrammetry',
    url='hg+https://bitbucket.org/devangel77b/python-sba',

    long_description=open(os.path.join(os.path.dirname(__file__),'README.txt'),'r').read(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Scientific/Engineering',
        'Topic :: Multimedia :: Video',
        'Topic :: Multimedia :: Graphics',
        'Topic :: Multimedia :: Graphics :: 3D Modeling',
        'Topic :: Multimedia :: Graphics :: Capture :: Digital Camera',
        'Operating System :: POSIX :: Linux',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research'
        ],
)
