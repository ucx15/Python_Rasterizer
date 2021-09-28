from RenderLib.RenderLib import *

w,h = 720,720

Obj = Mesh.LoadMesh("OBJ/torus.obj")
Obj2 = Mesh.LoadMesh("OBJ/sphere.obj")

Obj.Rotation = [-10,0,0]
Obj.Translate = [0,.3,6]
Obj2.Translate = [0,.3,6]
Obj.Scale = [.5]*3
Obj2.Scale = [.25]*3

Obj.mat = Material(color=Color(0,.2,1),spec=1)
Obj2.mat = Material(color=Color(1,.2,0),spec=1)

camr = Camera(Vec(), 50)

dxlt = Light(Vec(),power=.8,dirx=Vec(0,-.51,1))
pxlt = Light(Vec(0,2,0), 600,typ="POINT")


RenderScene = Scene([Obj,Obj2],[pxlt],camr,w,h)

Anim = [1,1,0]

#RenderScene.HDRClip=True
#RenderScene.DrwMode= "ALL"
