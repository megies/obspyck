import inspect
import os
import re
from setuptools import setup

INSTALL_REQUIRES = [
    'obspy>=1.1.0',
    # pyqt can not be declared as a dependency cleanly it seems, see
    # http://stackoverflow.com/questions/4628519/
    # 'PyQt5',
    'numpy',
    'scipy',
    'matplotlib',
    'requests',
    ]
ENTRY_POINTS = {
    'console_scripts': [
        'obspyck = obspyck.obspyck:main',
        ]
    }
PACKAGE_DATA = {
    'obspyck': ['example.cfg', 'obspyck.gif', 'obspyck_16x16.gif',
                'obspyck_24x24.gif', 'obspyck_32x32.gif', 'obspyck_48x48.gif']}

SETUP_DIRECTORY = os.path.dirname(os.path.abspath(inspect.getfile(
    inspect.currentframe())))


# get the package version from from the main __init__ file.
version_regex_pattern = r"__version__ += +(['\"])([^\1]+)\1"
for line in open(os.path.join(SETUP_DIRECTORY, 'obspyck',
                              '__init__.py')):
    if '__version__' in line:
        version = re.match(version_regex_pattern, line).group(2)


def find_packages():
    """
    Simple function to find all modules under the current folder.
    """
    modules = []
    for dirpath, _, filenames in os.walk(
            os.path.join(SETUP_DIRECTORY, "obspyck")):
        if "__init__.py" in filenames:
            modules.append(os.path.relpath(dirpath, SETUP_DIRECTORY))
    return [_i.replace(os.sep, ".") for _i in modules]

setup(
    name="obspyck",
    version=version,
    description="A GUI application for seismogram analysis",
    author="Tobias Megies",
    author_email="megies@geophysik.uni-muenchen.de",
    url="https://github.com/megies/obspyck",
    download_url="https://github.com/megies/obspyck.git",
    install_requires=INSTALL_REQUIRES,
    keywords=["obspy", "github", "seismology", "earthquake", "seismogram"],
    packages=find_packages(),
    package_data=PACKAGE_DATA,
    entry_points=ENTRY_POINTS,
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Topic :: Software Development :: Libraries :: Python Modules",
        ],
    long_description="ObsPyck is a GUI application that is intended to cover "
                     "the tasks in a standard analysis workflow for seismic "
                     "events in seismological observatory practice.",
    )
