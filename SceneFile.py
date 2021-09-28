from RenderLib.RenderLib import *

w,h = 1600,720

Obj = Mesh.LoadMesh("OBJ/teapot.obj")

Obj.Rotation = [30,0,0]
Obj.Scale = [1/1.3]*3
Obj.Translate = [0,0,6]
Obj.mat = Material(color=Color(0,0.2,1),spec=1)

camr = Camera(Vec(), 60)


dxlt = Light(Vec(),.8,dirx=Vec(0,-.5,1))
pxlt = Light(Vec(0,1,0), 250,typ="POINT")


RenderScene = Scene([Obj],[pxlt],camr,w,h)

Anim = [0,-2,0]
#RenderScene.HDRClip=True
#RenderScene.DrwMode= "WIREVERT"
