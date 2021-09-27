import pyximport; pyximport.install()

import pygame
from SceneFile2 import RenderScene,Anim
from RenderLib.RenderLibC import Render


if __name__  == "__main__":
	Render(pygame,RenderScene,Anim)	