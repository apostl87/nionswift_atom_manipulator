import setuptools
import os

setuptools.setup(
    name="nionswift_atom_manipulator",
    version="1.1",
    author="Andreas Postl",
    author_email="andreas.postl@univie.ac.at",
    description= "A Nion Swift plug-in for single atom manipulation using an STEM",
    packages=["nionswift_plugin.atom_manipulator", "nionswift_plugin.atom_manipulator.classes"],
    install_requires=[], # see requirements.txt
    license='GPLv3',
    classifiers=[
        "Development Status :: 1 - Alpha",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3.8",
    ],
    python_requires='~=3.8',
    zip_safe=False,
    dependency_links=[]
)
