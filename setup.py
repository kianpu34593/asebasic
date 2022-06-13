"""
BASIC
A Python package for Bulk Adsorption Surface energy calculation with automatIc Convergence
"""
import sys
from setuptools import setup

with open('VERSION', 'r') as f:
    VERSION = f.read()

with open('README.md', 'r') as f:
    LONG_DESCRIPTION = f.read()

# # from https://github.com/pytest-dev/pytest-runner#conditional-requirement
# needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
# pytest_runner = ['pytest-runner'] if needs_pytest else []

# try:
#     with open("README.md", "r") as handle:
#         long_description = handle.read()
# except:
#     long_description = "\n".join(short_description[2:])


setup(
    # Self-descriptive entries which should always be present
    name='asebasic',
    author='Jiankun Pu',
    author_email='jiankunp@andrew.cmu.edu',
    description='A python package for bulk, adsorption and surface energy calculation for GPAW in Atomic Simulation Environment.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    yrl="https://github.com/kianpu34593/asebasic"
    version=VERSION,
    #cmdclass=versioneer.get_cmdclass(),
    license='MIT',

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    package_dir={'':'src'}
    packages=["basic"]
    python_requires=">=3.6"
    install_requies=[
        #"gpaw",
        "ase",
        "pymatgen",
        "autocat"
    ]
    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    #include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    #setup_requires=[] + pytest_runner,

    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # url='http://www.my_package.com',  # Website
    # install_requires=[],              # Required packages, pulls from pip if needed; do not use for Conda deployment
    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    # python_requires=">=3.5",          # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
