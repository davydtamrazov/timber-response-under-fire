# Structural Response of Timber Structures under Fire

## Overview

The topic of timber structures was explored as a part of the [final project](//github.com/davydtamrazov/timber-response-under-fire/blob/main/aux/CEE282_FinalProject.pdf) for the Nonlinear Structural Analysis class (CEE282) at Stanford University. The code posted in this repository contains similar ideas to the original code, but fully rearranged and generalised to include the possibility of analysing any structural arrangement under fire.


#### Example of the analysis:

![](https://github.com/davydtamrazov/timber-response-under-fire/blob/main/aux/frame_example.gif)


## Methodology
A general framework of the implementation of the fire loading analysis of a timber structure is as follows:

1. Define timber section dimensions.
2. Define fire exposure time and surface temperature.
3. Define ambient timber material multi-linear properties.
4. Discretize section and split into fibers with unique temperature-adjusted material properties assigned to each.
5. Define geometry of the arrangement: nodes, elements, fixities.
6. Define point loads.
7. *Optional*: Define node recorders outputting results for a specific node in a .txt file.
8. Perform analysis with OpenSeesPy framework.

#### Section definition

Timber section is defined in [timber_section.py](https://github.com/davydtamrazov/timber-response-under-fire/blob/main/src/timber_section.py) file. The process is illustrated on the diagram below.

![](https://github.com/davydtamrazov/timber-response-under-fire/blob/main/aux/section_splitting_process.png)
