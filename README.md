# Python_Rasterizer

## Description 
This is an implementation of a Software Rasterizer written in Python.

## Dependencies
* PyGame
* `Python 3.x` is recommended.

## Notice
This is a pure Python 3.x code, and will be very slow on CPython interpreter.

For speed boost, kindly use PyPy version.

I've also implemented it in Cython, check that branch `Cython_Ver`. Cython code is much more verbose compared Pure Python due to Types, but in return you get blazing fast speed of C (Native C will almost always be faster)

## Building
* Clone this Repository
* Simply run `main.py` file
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
### Sept 28,2021
    * Initial Commit (see Features)

## Editing Scenes
_Kindly Read `RenderLib.py` file and if you understand a little bit of it, Scenes are defined by using a `Scene` object which contains all the data of scene._

I know this is not enough to explain but I'll be writing doc-strings in each every class and functions.. And in future, I maybe create a wiki or documentation for using this project 

## Contributions 
If you like this Project, and want to improve something in it, or just want to reimplement this from your perspective (any class/function/algorithm or whole project), feel free to do.. It will be very appreciated!
We're going to take over the VFX industry with this 🙂
