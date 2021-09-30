from RenderLib.RenderLibL import *

w,h = 1600,720

Obj = Mesh.LoadMesh("OBJ/teapot.obj")

Obj.Rotation = [30,0,0]
Obj.Scale = Vec(.5, .5, .5)
Obj.mat = Material(color=Color(0,0.2,1),spec=1)


camr = Camera(Vec(), 50)
pxlt = Light(Vec(0,1,0), 550)

RenderScene = Scene([Obj],[pxlt],camr,w,h)

Anim = [0,-2,0]
#RenderScene.HDRClip=True
#RenderScene.DrwMode= "WIRE"
