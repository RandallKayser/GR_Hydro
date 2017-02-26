import matplotlib as mp
import matplotlib.pyplot as pp
import numpy as np
from itertools import islice

f = open("output/first_test.txt", 'r')

i = 0
frame = 0
time_now = 0
xcellnum = 25
ycellnum = 25
rhovals_at_time_t = [0] * xcellnum * ycellnum

def makeframe(i):
	pp.clf()
	string = "vis/frame" + (4 - len(str(i))) * "0" + str(i) + ".png"
	pp.pcolormesh(rhovals_reshaped, cmap="RdBu", vmin=1., vmax=100.)
	pp.title("Density at Time = " + str(time_now))
	pp.ylabel("Position")
	pp.xlabel("Position")

	pp.savefig(string)

for line in f:
	if(i % (xcellnum * ycellnum) == 0):
		print i / (xcellnum * ycellnum)

	rhovals_at_time_t[i % (xcellnum * ycellnum)] = float(line.split(", ")[4])
	
	if i % (xcellnum * ycellnum) == 0:
		time_now = float(line.split(", ")[3])
	
	if i % (xcellnum * ycellnum) == xcellnum * ycellnum - 1:
		rhovals_reshaped = np.reshape(rhovals_at_time_t, [xcellnum, ycellnum])
		makeframe(frame)
		frame += 1
	
	i += 1
