from setuptools import setup, Extension, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

ext_modules = [Extension(
    'sattrack.sgp4',
    include_dirs = ['extensions/sgp4/include'],
    sources = [
        'extensions/sgp4/src/SGP4.cpp',
        'extensions/sgp4/src/sgp4module.cpp'
    ]
)]

setup(
    name = 'sattrack',
    version = '0.0.1',
    author = "Quinton Barnes",
    author_email = "devqbizzle68@gmail.com",
    description = "A satellites tracking framework.",
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    license= 'MIT',
    url = 'https://github.com/qbizzle68/sattrack',
    install_requires=['pyevspace>=0.0.8', 'requests'],
    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy'
    ],
    packages = ['sattrack', 'sattrack.rotation', 'sattrack.spacetime', 'sattrack.structures', 'sattrack.util'],
    package_dir = {'': 'src'},
    ext_modules = ext_modules,
)