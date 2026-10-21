# WRF-Chem-Training-2026
WRF-Chem Training 2026

Installation, Configuration and Introduction to Atmospheric Chemistry Modelling with WRF-Chem v4.5.2

Lancaster University | WRF-Chem Training Course | 2026
Author of the course: Dr Andrea Mazzeo a.mazzeo@lancaster.ac.uk

Overview

This repository provides training materials, installation resources, configuration scripts and supporting documentation for a practical course on the Weather Research and Forecasting model coupled with Chemistry (WRF-Chem).

The course is designed to introduce postgraduate students and early-career researchers to the fundamental concepts and practical procedures required to install, configure, execute and analyse simulations using WRF-Chem on a High-Performance Computing (HPC) system of the Lancaster University HEC.

The training focuses on the practical implementation of WRF-Chem, using the HEC HPC facility at Lancaster University as the reference computing environment.

A particular emphasis is placed on understanding the complete modelling workflow, from the installation of software dependencies and geographical domain configuration to the preparation of meteorological input data, model execution and analysis of atmospheric chemistry outputs.

The repository provides a prepared installation package containing WRF-Chem and WPS source archives, supporting libraries and installation scripts. This package was developed to simplify the installation process while allowing students to compile the software independently within their own HPC environment.


The repository provides course documentation and access to a downloadable student installation package.

The installation package is distributed through GitHub Releases rather than as an ordinary repository file.

Student installation package

Package: WRF-Chem_Course_Student_Package.tgz

Release tag: WRFchem

Approximate download size: 105 MB

The package includes:

WRF-Chem_Course_Student_Package.tgz
│
├── set_env.sh
│
├── scripts/
│   ├── install_netcdf.sh
│   ├── configure_wrfchem.sh
│   ├── compile_wrfchem.sbatch
│   ├── configure_wps.sh
│   ├── compile_wps.sbatch
│   ├── check_installation.sh
│   ├── prepare_from_working_install.sh
│   └── package_wps_source.sh
│
└── package/
    ├── WRFV4.5.2_Course_Source.tgz
    ├── WPS-4.5_Course_Source.tgz
    └── GRIB2_Course_Libraries.tgz

The source archives contain the prepared software needed for the training installation.

Software Attribution, Licensing and Disclaimer

This repository provides educational materials and installation resources for training in atmospheric modelling using WRF-Chem.

WRF-Chem and the WRF Preprocessing System (WPS) are developed and maintained by their respective scientific development communities. The original software, associated documentation and applicable copyright notices remain attributable to their original developers and contributors.

The WRF-Chem source included in the training package is based on WRF version 4.5.2 and contains limited modifications to support compilation in the Lancaster University Wayland HPC environment. These modifications are provided for educational purposes and do not constitute an official WRF-Chem release.

The original WRF public-domain notice is retained, and any third-party components remain subject to their own applicable licences and attribution requirements.

Any ERA2GRIB-related material included in this repository is intended for educational use, with attribution and redistribution governed by the terms applicable to its original source.

This repository is independently maintained for teaching and research training purposes and is not officially affiliated with, endorsed by or maintained by NCAR, UCAR or the original developers of the third-party software.

Users should consult the original software documentation and licence notices before redistributing or modifying individual components.
