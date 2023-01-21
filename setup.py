import sys

from setuptools import setup, Extension

sys.path.append('src')
import sattrack

description, long_description = sattrack.__doc__.split('\n', 1)

ext_modules = [Extension(
    'sattrack.sgp4',
    include_dirs=['extensions/sgp4/include'],
    sources=[
        'extensions/sgp4/src/SGP4.cpp',
        'extensions/sgp4/src/sgp4module.cpp'
    ]
)]

setup(
    name='sattrack',
    version='0.1.0',
    author="Quinton Barnes",
    author_email="devqbizzle68@gmail.com",
    description=description,
    long_description=long_description,
    long_description_content_type='text',
    license='MIT',
    url='https://github.com/qbizzle68/sattrack',
    install_requires=['pyevspace>=0.0.12.4', 'requests'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy'
    ],
    packages=['sattrack', 'sattrack.spacetime', 'sattrack.util'],
    package_dir={'': 'src'},
    ext_modules=ext_modules,
)
