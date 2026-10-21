## WRF-Chem-Training-2026
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

The source archives contain the prepared software needed for the training installation.



## Software Attribution, Licensing and Disclaimer

This repository provides educational materials and installation resources for training in atmospheric modelling using WRF-Chem.

WRF-Chem and the WRF Preprocessing System (WPS) are developed and maintained by their respective scientific development communities. The original software, associated documentation and applicable copyright notices remain attributable to their original developers and contributors.

The WRF-Chem source included in the training package is based on WRF version 4.5.2 and contains limited modifications to support compilation in the Lancaster University Wayland HPC environment. These modifications are provided for educational purposes and do not constitute an official WRF-Chem release.

The original WRF public-domain notice is retained, and any third-party components remain subject to their own applicable licences and attribution requirements.

Any ERA2GRIB-related material included in this repository is intended for educational use, with attribution and redistribution governed by the terms applicable to its original source.

This repository is independently maintained for teaching and research training purposes and is not officially affiliated with, endorsed by or maintained by NCAR, UCAR or the original developers of the third-party software.

Users should consult the original software documentation and licence notices before redistributing or modifying individual components.


## ERA5 Meteorological Data Processing — Modified NCAR era5_to_int

This repository includes a modified implementation of the NCAR `era5_to_int` software, developed to simplify the preparation of ERA5 meteorological forcing data for WRF and WRF-Chem simulations.

The original software was developed by the National Center for Atmospheric Research (NCAR) and is available at:

https://github.com/NCAR/era5_to_int

For the WRF-Chem Training 2026 course, the original workflow has been extended to integrate automated meteorological data retrieval and conversion within a single execution procedure.

The modifications include:

- Automated downloading of the required ERA5 meteorological input files.
- Integration of the download and conversion stages into a single workflow.
- Additional file handling and download checks.
- Adaptation of the execution procedure for the Lancaster University Wayland HPC environment.

The modified software is used to allow students to prepare the WPS intermediate files required by `metgrid.exe`.

**Licensing and attribution:** The modified implementation is derived from NCAR's `era5_to_int`, distributed under the Apache License, Version 2.0. Original copyright and licence notices are retained. Modifications introduced for the WRF-Chem Training 2026 course are independently maintained and are not an official NCAR release.

ERA5 data access and use remain subject to the applicable data-provider terms and conditions.
