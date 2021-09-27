from time import time
from math import sin,cos,tan,pi,radians

#vectors
cdef class Vec:
	cdef public double i,j,k
	def __init__(self, double i=0, double j=0, double k=0):
		self.i = i
		self.j = j
		self.k = k

	def __repr__(self):
		return f"{self.i} {self.j} {self.k}"
	
	def __neg__(self):
		return Vec(-self.i,-self.j,-self.k)
	
	def __add__(self, Vec v):
		return Vec((self.i + v.i), (self.j + v.j), (self.k + v.k))
	def __sub__(self, Vec v):
		return Vec((self.i - v.i), (self.j - v.j), (self.k - v.k))	
	
	def __mul__(self, scl):
		if type(scl) ==Vec:
			return Vec((self.i *scl.i), (self.j *scl.j), (self.k *scl.k))	
		return Vec((self.i *scl), (self.j *scl), (self.k *scl))	
	
	def __truediv__(self, double scl):
		invscl = 1/scl
		return Vec((self.i*invscl), (self.j*invscl), (self.k*invscl))
	
	def __eq__(self, Vec v):
		if (self.i == v.i) and (self.j==v.j) and (self.k==v.k):
			return True
		return False
	
	cpdef double mag(self):
		return (self.dot(self))**.5
	cpdef double magsq(self):
		return self.dot(self)
	
	cpdef Vec normalize(self):
		return self/(self.mag())
	
	cpdef double dot(self, Vec v):
		return ((self.i * v.i)+(self.j*v.j)+(self.k*v.k))	
	cpdef Vec cross(self, Vec v):
		return Vec((self.j*v.k - self.k*v.j), (self.k*v.i - self.i*v.k), (self.i*v.j - self.j*v.i))


#___Color___
cdef class Color:
	cdef public double r,g,b

	def __init__(self,double r=0, double g=0, double b=0):
		self.r = r
		self.g = g
		self.b = b
	
	def __repr__(self):
		return f"{self.r} {self.g} {self.b}"
	def __add__(self, Color c):
		return Color((self.r + c.r), (self.g + c.g), (self.b + c.b))
	def __iadd__(self, Color c):
		return self+c	
	def __mul__(self, scl):
		if type(scl) == Color:
			return Color((self.r *scl.r), (self.g *scl.g), (self.b *scl.b))
		else:
			return Color((self.r *scl) , (self.g *scl), (self.b *scl))	

	def __truediv__(self, double scl):
		return self/scl
	
	def __eq__(self, Color c):
		if ((self.r == c.r) and (self.g==c.g) and (self.b==c.b)):
			return True
		return False
	
	cpdef list quantize(self,int bpc=8, int clip=0):	
		cdef int qv = 2**bpc -1
		cdef double gam = .45
		cdef double rh,gh,bh
		rh,gh,bh = self.r, self.g, self.b

		if clip:
			if rh>=1 or gh>=1 or bh>=1:
				return [0,0,255]
		
		rh = (rh**gam)*qv
		gh = (gh**gam)*qv
		bh = (bh**gam)*qv

		return [min(rh,qv),
				min(gh,qv),
				min(bh,qv) ]


#__Material__
cdef class Material:
	cdef public Color color
	cdef public double rough,spec

	def __init__(self,Color color=Color(), double rough=0, double spec=0):
		self.color = color
		self.rough = rough
		self.spec = spec



#__Geometry__
cdef class Triangle:
	cdef public Vec a,b,c,BA,CA
	cdef public Material mat

	def __init__(self, Vec a, Vec b, Vec c):
		self.a, self.b, self.c = a,b,c
		self.BA,self.CA = (b-a),(c-a)
		self.mat = Material(color=Color())

	cpdef Vec normal(self):
		return (self.BA.cross(self.CA))
	cpdef Vec center(self):
		return (self.a+self.b+self.c)*.33


cdef class Mesh:
	cdef public list triangles,Rotation,Scale,PreTranslate,Translate
	cdef public Material mat
	
	def __init__(self,list fcs):
		self.triangles = fcs
		self.mat = Material(Color(1,1,1))
		
		self.Rotation = [0,0,0]
		self.Scale = [1,1,1]		
		self.PreTranslate = [0,0,0]
		self.Translate = [0,0,5]
	
	@staticmethod
	def LoadMesh(path,norm=1):
		cdef list verts,tris
		verts,tris = [],[]
		cdef int ctr = 0

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

cdef class Camera:
	cdef public double fov
	cdef public int typ
	cdef public Vec loc

	def __init__(self, loc, fov):
		self.loc = loc
		self.fov = fov
		self.typ = 1


cdef class Scene:
	cdef public list Objects,Lights,Size,VertColor,EdgeColor
	cdef public Camera camera
	cdef public int W,H,HDRClip,Lighting
	cdef public double Scale,NearPoint,FarPoint
	cdef public str DrwMode
	
	def __init__(self,list objs, list lt, Camera cam, int w, int h):
		self.Objects = objs
		self.Lights = lt
		self.camera = cam
		self.W, self.H = w,h
		
		self.Size = [self.W,self.H]
		self.Scale = 100
		
		self.NearPoint = .1
		self.FarPoint = 1e4
		
		self.Lighting = 1
		self.DrwMode = "RENDER"
		
		self.VertColor = [0,120,255]
		self.EdgeColor = [250,240,240]
		
		self.HDRClip = 0


cdef class Light:
	cdef public Vec loc,dirx
	cdef public Color color
	cdef public str typ
	cdef public double power

	def __init__(self,Vec loc, double power,Vec dirx=Vec(0,0,1),str typ="DIRX"):
		self.loc = loc
		self.power = power
		self.typ = typ
		self.dirx = dirx
		self.color = Color(1,1,1)


#___FUNCTIONS___
cdef list PersMatrixGen(fov, AsR, zN, zF):
	fov /= 2
	f = 1/tan(radians(fov))
	q = zF/(zF-zN)
	return [AsR*f,f,q,zN*q]	


cdef list RotMatrixGen(AnTup):
	a1,a2,a3 = map(radians,AnTup)
	return [Vec((cos(a2)*cos(a3)),(cos(a2)*sin(a3)),(-sin(a2))),
			Vec((sin(a1)*sin(a2)*cos(a3) - cos(a1)*sin(a3)),(sin(a1)*sin(a2)*sin(a3)+cos(a1)*cos(a3)),(sin(a1)*cos(a2))),
			Vec((cos(a1)*sin(a2)*cos(a3) - sin(a1)*sin(a3)),(cos(a1)*sin(a2)*sin(a3)-sin(a1)*cos(a3)),(cos(a1)*cos(a2)))]


cdef dict Transform(Mesh O, int NrmOp, Camera cmra):

	cdef Triangle T,FinalTris
	cdef Vec a,b,c, r1,r2,r3,rVec
	cdef double sx,sy,sz, rtx,rty,rtz

	sx,sy,sz = O.Scale
	rtx,rty,rtz = O.PreTranslate				
	tx,ty,tz = O.Translate		
	r1,r2,r3 = RotMatrixGen(O.Rotation)
	cdef Vec PreTr = Vec(rtx,rty,rtz)

	cdef dict LocTrisDic = {}
	cdef Material mt = O.mat
	
	for T in O.triangles:	
		#_Transformation
		a,b,c = T.a+PreTr, T.b+PreTr, T.c+PreTr
		FinalTris = Triangle(
					Vec(r1.dot(a)*sx+tx,r2.dot(a)*sy+ty,r3.dot(a)*sz+tz),
					Vec(r1.dot(b)*sx+tx,r2.dot(b)*sy+ty,r3.dot(b)*sz+tz),
					Vec(r1.dot(c)*sx+tx,r2.dot(c)*sy+ty,r3.dot(c)*sz+tz))
					
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



cdef list Project(list PM, Triangle tr, int cx, int cy, double Scl, int Typ):
	cdef double p1x,p2x,p3x
	cdef double p1y,p2y,p3y
	cdef double p1z,p2z,p3z
	cdef double z1,z2,z3
	cdef double Px1,Px2,Px3,Py1,Py2,Py3,Pz1,Pz2,Pz3

	p1x,p1y,p1z = tr.a.i, tr.a.j, tr.a.k
	p2x,p2y,p2z = tr.b.i, tr.b.j, tr.b.k	
	p3x,p3y,p3z = tr.c.i, tr.c.j, tr.c.k	
	
	if Typ:
		Px1,Py1,Pz1 = (PM[0]*p1x,
						PM[1]*p1y,
						PM[2]*p1z - PM[3])
		Px2,Py2,Pz2 = (PM[0]*p2x,
						PM[1]*p2y,
						PM[2]*p2z - PM[3])
		Px3,Py3,Pz3 = (PM[0]*p3x,
						PM[1]*p3y,
						PM[2]*p3z - PM[3])

		if p1z:
			z1 = 1/p1z
			Px1,Py1,Pz1 = Px1*z1, Py1*z1, Pz1*z1
		if p2z:
			z2 = 1/p2z
			Px2,Py2,Pz2 = Px2*z2, Py2*z2, Pz2*z2
		if p3z:
			z3 = 1/p3z
			Px3,Py3,Pz3 = Px3*z3, Py3*z3, Pz3*z3

		return [
		Vec( (Px1+1)*cx,(1-Py1)*cy,Pz1 ),
		Vec( (Px2+1)*cx,(1-Py2)*cy,Pz2 ),		
		Vec( (Px3+1)*cx,(1-Py3)*cy,Pz3 )]

	else:
		return [
				Vec(p1x*Scl + cx, cy-p1y*Scl),
				Vec(p2x*Scl + cx, cy-p2y*Scl),
				Vec(p3x*Scl + cx, cy-p3y*Scl)]



cdef Color Shade(Triangle T, Light L, Camera C):
	cdef Vec TNor,Lv,H
	cdef Color RetCol,LCor
	cdef double LAtten,DiffConst,SpecConst

	LCor = L.color
	LV = -L.dirx
	cdef str Lty = L.typ
	
	TNor = (T.normal()).normalize()
	RetCol = Color(0,0,0)

	LAtten = 1
	
	#dif
	if Lty == "DIRX":
		DiffConst = ( (TNor).dot(LV.normalize()) )

	else:
		LVec = L.loc - T.center()
		LAtten = 1/(4*pi*LVec.magsq())		
		DiffConst = (TNor).dot(LVec.normalize())

	if DiffConst > 0:
		RetCol = T.mat.color*(DiffConst*min(L.power,L.power*LAtten))
		if LCor != Color(1,1,1):
			RetCol = RetCol*LCor

		#spec
		if T.mat.spec:
			if Lty == "DIRX":
				H = ((C.loc-T.a)+LV).normalize()
			else:
				Lv = L.loc-T.a
				H = ((C.loc-T.a)+LV).normalize()
			

			SpecConst = (H.dot(TNor))**150
			if SpecConst >0:
				RetCol += L.color*SpecConst
	return RetCol



def Render(PG, Scene UserScene, list RotAn):
	PG.init()
	myfont = PG.font.SysFont("monospace", 15)
	
	Display = PG.display
	Display.set_caption("URaster")
	Surface = Display.set_mode(UserScene.Size)
	clock = PG.time.Clock()
	
	cdef Camera Cam = UserScene.camera
	cdef list Sze = UserScene.Size
	cdef int OX,OY 
	OX,OY = (Sze[0]/2), (Sze[1]/2)
	cdef double Scl = UserScene.Scale
	cdef list Lts = UserScene.Lights

	cdef int VertSize = 2
	cdef int EdgeSize = 1
	cdef int NOpt,DVt,DEd,DFc
	
	NOpt = 1
	DVt = DEd = DFc = 0
	DrwMode = UserScene.DrwMode
	if DrwMode == "RENDER":
		DFc = 1
	elif DrwMode == "ALL":
		DVt = DEd = DFc = 1
	elif DrwMode == "ALL-NV":
		DEd = DFc = 1
	elif DrwMode == "ALL-NE":
		DVt = DFc = 1
	elif DrwMode == "WIRE":
		DEd = 1
		NOpt = 0
	elif DrwMode == "WIREVERT":
		DEd = DVt = 1
		NOpt = 0
	elif DrwMode == "VERT":
		DVt = 1
		NOpt = 0
	
	cdef list VColor,EColor
	VColor,EColor = UserScene.VertColor, UserScene.EdgeColor

		
	cdef list PersMatrix 
	PersMatrix = PersMatrixGen(Cam.fov,
							(<double>UserScene.H/<double>UserScene.W),
							UserScene.NearPoint,
							UserScene.FarPoint)

	cdef dict UnsortedTris
	cdef list SortedTris
	cdef Triangle TRIS
	cdef Vec a,b,c
	cdef Color SurfCol
	cdef int R,G,B
	cdef Mesh Object

	while 1:
		t1 = time()
		clock.tick(60)
		Surface.fill((0,0,0))
			
		UnsortedTris = {}
		
		tTrans1 = time()
		#__GEOMETRY_CALCULATION___	
		for Object in UserScene.Objects:
			Object.Rotation = list(map(sum, zip(Object.Rotation, RotAn)))		
			UnsortedTris.update(Transform(Object,NOpt,Cam))
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


		#__PERFORMANCE__
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
