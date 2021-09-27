from RenderLib.RenderLibC import *

w,h = 1200,720

Obj = Mesh.LoadMesh("OBJ/sphere2.obj")

Obj.mat.color = Color(1,0,0)
Obj.Rotation = [0,0,0]
Obj.Scale = [.5]*3
Obj.PreTranslate = [0,0,0]
Obj.Translate = [0,0,6]

Obj.mat = Material(color=Color(1,0,.1),spec=1)

camr = Camera(Vec(), 50)

dxlt = Light(Vec(),power=.8,dirx=Vec(0,-.51,1))
pxlt = Light(Vec(0,2,0), 280,typ="POINT")


RenderScene = Scene([Obj],[pxlt],camr,w,h)

Anim = [0,1,0]
#RenderScene.HDRClip=True
#RenderScene.DrwMode= "ALL"