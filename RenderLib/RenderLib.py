from time import time
from math import sin,cos,tan,pi,radians

#vectors
class Vec:
	def __init__(self, i=0,j=0,k=0):
		self.i = i
		self.j = j
		self.k = k

	def __repr__(self):
		return f"{self.i} {self.j} {self.k}"
	
	def __neg__(self):
		return Vec(-self.i,-self.j,-self.k)
	
	def __add__(self, v):
		return Vec((self.i + v.i), (self.j + v.j), (self.k + v.k))
	def __sub__(self, v):
		return Vec((self.i - v.i), (self.j - v.j), (self.k - v.k))	
	
	def __mul__(self, scl):
		if type(scl) ==Vec:
			return Vec((self.i *scl.i), (self.j *scl.j), (self.k *scl.k))	
		return Vec((self.i *scl), (self.j *scl), (self.k *scl))	
	
	def __truediv__(self,scl):
		invscl = 1/scl
		return Vec((self.i*invscl), (self.j*invscl), (self.k*invscl))
	
	def __eq__(self, v):
		if (self.i == v.i) and (self.j==v.j) and (self.k==v.k):
			return True
		else:
			return False
	
	def mag(self):
		return (self.dot(self))**.5
	def magsq(self):
		return self.dot(self)
	def normalize(self):
		return self/(self.mag())
	def dot(self, v):
		return ((self.i * v.i)+(self.j*v.j)+(self.k*v.k))	
	def cross(self, v):
		return Vec((self.j*v.k - self.k*v.j), (self.k*v.i - self.i*v.k), (self.i*v.j - self.j*v.i))


#___Color___
class Color:
	def __init__(self, r=0,g=0,b=0):
		self.r = r
		self.g = g
		self.b = b
	
	def __repr__(self):
		return f"{self.r} {self.g} {self.b}"
	def __add__(self, c):
		return Color((self.r + c.r), (self.g + c.g), (self.b + c.b))
	def __iadd__(self, c):
		return Color((self.r + c.r), (self.g + c.g), (self.b + c.b))	
	def __mul__(self, scl):
		if type(scl) == Color:
			return Color((self.r *scl.r), (self.g *scl.g), (self.b *scl.b))
		else:
			return Color((self.r *scl) , (self.g *scl), (self.b *scl))	

	def __truediv__(self, scl):
		return self*(1/scl)
	
	def __eq__(self, c):
		if ((self.r == c.r) and (self.g==c.g) and (self.b==c.b)):
			return True
		return False
	def quantize(self,bpc=8,clip=False):	
		qv = 2**bpc -1
		gam = .45
		if clip:
			if self.r>=1 or self.g>=1 or self.b>=1:
				return [0,0,255]
				
		r,g,b= ((min(int((self.r**gam)*qv),qv)),
				(min(int((self.g**gam)*qv),qv)),
				(min(int((self.b**gam)*qv),qv)) )
		return [r,g,b]

#__Material__
class Material:
	def __init__(self,color=Color(),rough=0,spec=0):
		self.color = color
		self.rough = rough
		self.spec = spec



#__Geometry__
class Triangle:
	def __init__(self,a,b,c):
		self.a, self.b, self.c = a,b,c
		self.color = Color(1,1,1)
		self.mat = Material()

	def normal(self):
		return ((self.b-self.a).cross(self.c-self.a))
	def center(self):
		return Vec(
		(self.a.i+self.b.i+self.c.i)*.33,
		(self.a.j+self.b.j+self.c.j)*.33,
		(self.a.k+self.b.k+self.c.k)*.33)


class Mesh:
	def __init__(self,fcs):
		self.name = ""
		self.triangles = fcs
		self.mat = Material(Color(1,1,1))
		
		self.Rotation = [0,0,0]
		self.Scale = [1,1,1]
		self.Translate = [0,0,5]
	
	@staticmethod
	def LoadMesh(path,norm=1):
		verts = []
		tris = []
		ctr = 0
		with open(path, "r") as file:
			Lines = file.readlines()
			file.close()
	
		for Line in Lines:
			if Line.startswith("v "):
				a,b,c = list(map(float,(Line[2::]).split()))
				verts.append(Vec(a,b,c))
	
	
		for i,Line in enumerate(Lines):
			
			if Line.startswith("f"):
				a,b,c = Line[2::].split()
	
				if "//" in Line:			
					a,_ = a.split("//")
					b,_ = b.split("//")
					c,_ = c.split("//")
	
				if "/" in Line:			
					a = a.split("/")[0]
					b = b.split("/")[0]
					c = c.split("/")[0]
		
				a,b,c = int(a),int(b),int(c)
	
				try:
					if norm == -1:
						tris.append(Triangle(verts[a-1],verts[c-1],verts[b-1]))
					else:
						tris.append(Triangle(verts[a-1],verts[b-1],verts[c-1]))
						
				except IndexError:
					ctr += i
	
		if ctr: print("faulty",ctr)			
		return Mesh(tris)	
	

#__Camera_Scene_&_Lights__

class Camera:
	def __init__(self, loc, fov):
		self.loc = loc
		self.fov = fov
		self.typ = 1


class Scene:
	def __init__(self,objs,lt,cam,w,h):
		self.Objects = objs
		self.Lights = lt
		self.camera = cam
		self.W, self.H = w,h
		
		self.Size = (self.W,self.H)
		self.Scale = 100
		
		self.NearPoint = .1
		self.FarPoint = 1e4
		
		self.Lighting = True
		self.DrwMode = "RENDER"
		
		self.VertColor = (0,120,255)
		self.EdgeColor = (250,240,240)
		
		self.HDRClip = False


class Light:
	def __init__(self, loc, power, dirx=Vec(0,0,1), typ="DIRX"):
		self.loc = loc
		self.power = power
		self.typ = typ
		self.dirx = dirx
		self.color = Color(1,1,1)


#___FUNCTIONS___
def PersMatrixGen(fov, AsR, zN, zF):
	fov /= 2
	f = round(1/tan(radians(fov)),8)
	q = round(zF/(zF-zN),8)
	return (AsR*f,f,q,zN*q)	


def RotMatrixGen(AnTup):
	a1,a2,a3 = map(radians,AnTup)
	return (Vec((cos(a2)*cos(a3)),(cos(a2)*sin(a3)),(-sin(a2))),
		Vec((sin(a1)*sin(a2)*cos(a3) - cos(a1)*sin(a3)),(sin(a1)*sin(a2)*sin(a3)+cos(a1)*cos(a3)),(sin(a1)*cos(a2))),
		Vec((cos(a1)*sin(a2)*cos(a3) - sin(a1)*sin(a3)),(cos(a1)*sin(a2)*sin(a3)-sin(a1)*cos(a3)),(cos(a1)*cos(a2))))


def Transform(O,Anim,NrmOp,cmra):
	RotATup = O.Rotation
	RotATup[0] += Anim[0]
	RotATup[1] += Anim[1]
	RotATup[2] += Anim[2]
	
	sx,sy,sz = O.Scale				
	tx,ty,tz = O.Translate		
	r1,r2,r3 = RotMatrixGen(RotATup)
	LocTrisDic = {}
	mt = O.mat
	
	for T in O.triangles:	
		#_Transformation
		a,b,c = T.a,T.b,T.c
		FinalTris = Triangle(
			(Vec(r1.dot(a)*sx+tx,r2.dot(a)*sy+ty,r3.dot(a)*sz+tz) ),
			(Vec(r1.dot(b)*sx+tx,r2.dot(b)*sy+ty,r3.dot(b)*sz+tz) ),
			(Vec(r1.dot(c)*sx+tx,r2.dot(c)*sy+ty,r3.dot(c)*sz+tz) ))

		#_Apply_Colors
		FinalTris.mat = mt

		#_BackFace_Culling
		rVec = (cmra.loc - FinalTris.a)
		if NrmOp:
			if rVec.dot(FinalTris.normal()) > 0:
				LocTrisDic[FinalTris] = rVec.mag()
		else:
			LocTrisDic[FinalTris] = rVec.mag()

	return LocTrisDic	



def Project(PM,tr,cx,cy,Scl,Typ):
	p1,p2,p3 = tr.a, tr.b, tr.c
	
	if Typ:
		Px1,Py1,Pz1 = (PM[0]*p1.i,
				PM[1]*p1.j,
				PM[2]*p1.k - PM[3])
		Px2,Py2,Pz2 = (PM[0]*p2.i,
				PM[1]*p2.j,
				PM[2]*p2.k - PM[3])
		Px3,Py3,Pz3 = (PM[0]*p3.i,
				PM[1]*p3.j,
				PM[2]*p3.k - PM[3])

		if p1.k:
			z1 = 1/p1.k
			Px1,Py1,Pz1 = Px1*z1, Py1*z1, Pz1*z1
		if p2.k:
			z2 = 1/p2.k
			Px2,Py2,Pz2 = Px2*z2, Py2*z2, Pz2*z2
		if p3.k:
			z3 = 1/p3.k
			Px3,Py3,Pz3 = Px3*z3, Py3*z3, Pz3*z3
		
		return (
		Vec( (Px1+1)*cx,(-Py1+1)*cy,Pz1 ),
		Vec( (Px2+1)*cx,(-Py2+1)*cy,Pz1 ),		
		Vec( (Px3+1)*cx,(-Py3+1)*cy,Pz1 ))

	else:
		return (
			Vec(p1.i*Scl + cx, cy-p1.j*Scl),
			Vec(p2.i*Scl + cx, cy-p2.j*Scl),
			Vec(p3.i*Scl + cx, cy-p3.j*Scl))



def Shade(T,L,C):
	TNor = (T.normal()).normalize()
	RetCol = None
	LCor = L.color
	
	#dif
	LAtten = 1
	if L.typ == "DIRX":
		DiffConst = ( (TNor).dot((-L.dirx).normalize()) )

	elif L.typ == "POINT":
		LVec = L.loc - T.center()
		LAtten = 1/(4*pi*LVec.magsq())		
		DiffConst = (TNor).dot(LVec.normalize())

	if DiffConst > 0:
		RetCol = T.mat.color*(DiffConst*min(L.power,L.power*LAtten))
		if LCor != Color(1,1,1):
			RetCol = RetCol*LCor

		#spec
		if T.mat.spec:
			if L.typ == "DIRX":
				lv = -L.dirx
				H = ((C.loc-T.a)+lv).normalize()
			elif L.typ =="POINT":
				lv = L.loc-T.a
				H = ((C.loc-T.a)+lv).normalize()
			

			SpecConst = (H.dot(TNor))**150
			if SpecConst >0:
				RetCol += L.color*SpecConst
	return RetCol



def Render(PG,UserScene,RotAn):
	PG.init()
	myfont = PG.font.SysFont("monospace", 15)
	
	Display = PG.display
	Display.set_caption("URaster")
	Surface = Display.set_mode(UserScene.Size)
	clock = PG.time.Clock()
	
	Cam = UserScene.camera
	Sze = UserScene.Size
	OX,OY = Sze[0]/2, Sze[1]/2
	Scl = UserScene.Scale
	Lts = UserScene.Lights

	VertSize = 2
	EdgeSize = 1
	NOpt = True
	
	DVt = DEd = DFc = False
	DrwMode = UserScene.DrwMode
	if DrwMode == "RENDER":
		DFc = True
	elif DrwMode == "ALL":
		DVt = DEd = DFc = True
	elif DrwMode == "ALL-NV":
		DEd = DFc = True
	elif DrwMode == "ALL-NE":
		DVt = DFc = True
	elif DrwMode == "WIRE":
		DEd = True
		NOpt = False
	elif DrwMode == "WIREVERT":
		DEd = DVt = True
		NOpt = False
	elif DrwMode == "VERT":
		DVt = True
		NOpt = False
	
	VColor,EColor = (UserScene.VertColor,
			UserScene.EdgeColor)

		
	PersMatrix = PersMatrixGen(Cam.fov,
				(UserScene.H/UserScene.W),
				UserScene.NearPoint,
				UserScene.FarPoint)


	while True:
		t1 = time()
		clock.tick(60)
		Surface.fill((0,0,0))
			
		UnsortedTris = {}
		
		tTrans1 = time()
		#__GEOMETRY_CALCULATION___	
		for Object in UserScene.Objects:
			UnsortedTris.update(Transform(Object,RotAn,NOpt,Cam))
		tTrans2 = time()



		#__FACE_SORTING___
		SortedTris = sorted(UnsortedTris, key=UnsortedTris.get)[::-1]

				
		tPrj = RastT = tLgt = 0
		#___RENDERING__
		for TRIS in SortedTris:
			#_Projection_
			tPrj1 = time()
			a,b,c = Project(PersMatrix,TRIS,OX,OY,Scl,Cam.typ)
			tPrj2 = time()
			tPrj += (tPrj2-tPrj1)
			
			#_Lighting
			if DFc and Lts:
				for light in Lts:
					
					SurfCol = Shade(TRIS,light,Cam)
					tlgt2 = time()
					tLgt += tlgt2-tPrj2
					
					if SurfCol:
						R,G,B = SurfCol.quantize(clip=UserScene.HDRClip)		
						PG.draw.polygon(Surface, (R,G,B), ((a.i,a.j),(b.i,b.j),(c.i,c.j)))
	
	
			if DEd:
				PG.draw.lines(Surface, EColor, 1, ((a.i,a.j),(b.i,b.j),(c.i,c.j)), EdgeSize)

			if DVt:		
				PG.draw.circle(Surface, VColor, (a.i,a.j), VertSize)
				PG.draw.circle(Surface, VColor, (b.i,b.j), VertSize)
				PG.draw.circle(Surface, VColor, (c.i,c.j), VertSize)
			tRen2 = time()
			RastT += (tRen2 - tPrj2)
		t2 = time()


		#__PROFILING__
		tt = (t2-t1)
		TransT = round((tTrans2-tTrans1)*100,3)
		tPrj = round(tPrj*100,3)
		tLgt = round(tLgt*100,3)
		RastT = round(RastT*100-tLgt,3)

		fps = round(1/tt,2)
		ttris = len(SortedTris)
		
		fpsLab = myfont.render(f"{ttris}Tri  {fps}fps  {round((tt*100),2)}:Elapsed", 1, (255,255,0))
		timerLab = myfont.render(f"{TransT}:Transform  {tPrj}:Projection  {tLgt}:Lighting  {RastT}:Rendering", 1, (255,255,0))	

		Surface.blit(fpsLab, (20, Sze[1]-220))
		Surface.blit(timerLab, (20, Sze[1]-200))

		Display.update()

#		break


