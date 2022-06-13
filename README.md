# nionswift_atom_manipulator

Single-atom manipulation tool for Nion Swift

-----
**Prerequisite**
--
Anaconda environment with Nion Swift installed (see https://nionswift.readthedocs.io/en/stable/installation.html)

-----
**Installation**
--

1. Activate the Anaconda environment with a Nion Swift installation. Install the required Python packages, which are listed below (see also requirements.txt). These will not be installed in step 3.
    - numpy
    - scipy
    - matplotlib
    - scikit-image
    - *periodictable* (optional, for elemental identification)
    - double_gaussian_blur (see https://github.com/arpostl/double_gaussian_blur)
    - fourier_scale_calibration (see https://github.com/jacobjma/fourier-scale-calibration)
    - nionswift_structure_recognition (see https://github.com/jacobjma/nionswift-structure-recognition)
    - nionswift_univie_tractorbeam (see [T.B.A.])
    - *nionswift-usim* (fork, see https://github.com/jacobjma/nionswift-usim) (optional, for simulation mode)

2. In a terminal (Linux, MacOS) or command prompt (Windows) window, navigate to the root folder of this package.

3. Install this package by executing the provided setup.py:
```
$ python3 ./setup.py install
```
or
```
$ pip3 install .
```

-----
**Infographics**
--
![Task overview and description](./infographics/tasks-and-descriptions.png)
--
![Operating modes](./infographics/operating-modes.png)

-----
**Screenshot**
--
![Plug-in screenshot](./infographics/plugin-screenshot.png)

-----
**Acknowledgements**
--

Cordial thanks for coding and support go to
- Toma Susi (https://github.com/TomaSusi)
- Jacob Madsen (https://github.com/jacobjma)
- Andreas Mittelberger (https://github.com/Brow71189)
- Christoph Hofer (https://github.com/christophhofer40)
- Nion company (https://github.com/nion-software)
