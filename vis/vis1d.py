import matplotlib as mp
import matplotlib.pyplot as pp
import numpy as np

f = open("output/first_test.txt", 'r')

myarray = np.loadtxt(f, delimiter=", ", usecols=(0, 3, 4, 5, 6, 7))

myarray = np.transpose(myarray)
xcellnum = 250
time = myarray[1]
rhovals = myarray[2]
pvals = myarray[3]
u0vals = myarray[4]
u1vals = myarray[5]
vvals = [u1vals[i] / u0vals[i] for i in xrange(len(u0vals))]
frameskip = max(len(time) / xcellnum / 250, 1)
def makeframe(i):
	pp.clf()
	string = "vis/frame" + (4 - len(str(i))) * "0" + str(i) + ".png"
	x = np.linspace(0., 1., num=xcellnum)
	pp.subplot(411)
	pp.plot(x, rhovals[xcellnum * (i * frameskip): xcellnum * (i * frameskip +1)], "r-")
	pp.axis([0., 1., 0., 15.])
	pp.title("Time = " + str(time[xcellnum * i]))
	pp.ylabel("Density")

	pp.subplot(412)
	pp.plot(x, pvals[xcellnum * i * frameskip: xcellnum * (i * frameskip + 1)], "r-")
	pp.axis([0., 1., 0., 18.])
	pp.ylabel("Pressure")

	pp.subplot(413)
	pp.plot(x, u0vals[xcellnum * i * frameskip: xcellnum * (i * frameskip + 1)], "r-")
	pp.axis([0., 1., 0., 2.])
	pp.ylabel("lorentz factor")

	pp.subplot(414)
	pp.plot(x, vvals[xcellnum * i * frameskip: xcellnum * (i * frameskip + 1)], "r-")
	pp.ylabel("v")
	pp.axis([0., 1., -1., 1.])
	pp.xlabel("Position");

	pp.savefig(string)

for i in xrange(len(time) / xcellnum / frameskip):
	makeframe(i)
