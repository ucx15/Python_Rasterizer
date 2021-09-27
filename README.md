# Python_Rasterizer

## Description 
This is an implementation of a Software Rasterizer written in Cython.

## Dependencies
* PyGame
* Cython 
* `Python 3.x` is recommended.

## Notice
This is the Cython implementation you're looking for

## Building
* Clone this Repository
* Simply run `main.py` file
* You can either use pyximport or build the `renderlib.pyx` file using cython compiler
* Changes in the scene can be made by editing `scenefile.py` or `scenefile2.py`, whichever is imported in `main.py`

## Features
* Perspective and Orthographic Projection
* Rotation, Translation and Scale
* Real-time Preview using PyGame
* Diffuse and Specular Shading
* Directional and Point Lighting
* Real-time Profiling
* OBJ import support (only triangluar geometry is supported for now)

## ChangeLog
###Sept 28,2021
    * Initial Commit (see Features)
