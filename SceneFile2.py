from RenderLib.RenderLib import *

w,h = 720,720

Obj = Mesh.LoadMesh("OBJ/torus.obj")
Obj2 = Mesh.LoadMesh("OBJ/sphere.obj")

Obj.Rotation = [-10,0,0]
Obj.Scale = Vec(.5, .5, .5)
Obj2.Scale = Vec(.25, .25, .25)

Obj.mat = Material(color=Color(0,.2,1),spec=1)
Obj2.mat = Material(color=Color(1,.2,0),spec=1)

camr = Camera(Vec(), 50)
pxlt = Light(Vec(0,2,0), 600,typ="POINT")


RenderScene = Scene([Obj,Obj2],[pxlt],camr,w,h)

Anim = [1,1,0]

#RenderScene.HDRClip=True
#RenderScene.DrwMode= "ALL"
