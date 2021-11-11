"""
31Oct2021 -cutkosky
Compare the linear algebra nullspace Grasp matrix
with what Yoshikawa & Nagai (1991 TRO) suggest. The latter approach
is what we implemented for controlling our two-fingered hand
(Griffin 2003 PhD thesis).

As usual, it's best to select and run a few lines at a time, inspecting
intermediate results to make sure they make sense to you.
"""
import numpy as np
import matplotlib.pyplot as plt

"""
Consider a planar block of width = 2 units, held on the left and
right by two fingers with frictional point contacts. To keep things
simple, let the finger frames be aligned with the body frame with X
axes all pointing to the right and Y axes pointing upward.
"""

"""
PART 1: Linear Mapping from fbody to finger forces.
Set up a planar grasp matrix in the usual way to map finger
forces [fx1,fy1,fx2,fy2]' to a body wrench [fxb, fyb, mzb]' and an
internal force fint (very similar to FlatBot). Finger 1 is on the left.
What should the internal force be? A 'reasonable' solution
could be to make it the average of the inward-pointing contact forces: 
"""
Gtransp = np.array([[1,0,1,0],[0,1,0,1],[0,-1,0,1],[0.5,0,-0.5,0]])

"""
Suppose we have a downward force in the Y direction of 1N and a coefficient 
of friction mu = 0.5  An inward force of 1N on each finger (+1N in
in the X direction on the left, -1N on the right) and an upward force of 0.5N
at each contact satisfies the friction constraint.
"""
ffingers = np.array([1,0.5,-1,0.5])
fbody = np.dot(Gtransp,ffingers)
#fbody = [fbx,fby,mbz,fint] = [0, 1, 0, 1]

"""
Now suppose the external force can be in any direction. What happens if 
we use the above definition of Gtransp and its inverse to solve for 
the finger forces?
"""
fext = 1       #You can experiment with different external force magnitudes
fint = 1       #You can experiment with different fint magnitudes
fbody[2] = 0   #External moment (if any)
fbody[3] = fint

#We need inverse grasp matrix to go from body to finger forces
Ginv = np.linalg.inv(Gtransp)

#Consider a range of directions for the external force
numsteps = 36
thetas = np.linspace(0,2*np.pi,numsteps)
plotpts = np.zeros((numsteps,2))
markers = ['.'] * numsteps

for i in range(0,numsteps):
    fbody[0] = fext*np.cos(thetas[i])
    fbody[1] = fext*np.sin(thetas[i])
    ffingers = np.dot(Ginv,fbody)
    plotpts[i,0] = ffingers[0]
    plotpts[i,1] = ffingers[1]
#Highlight any cases where friction constraint is not satisfied
    if (np.abs(ffingers[1])/ffingers[0] > 0.5):
        markers[i] = 'x'  

#Plot the forces on the left finger as a function of theta
plt.figure(1)
fig1 = plt.gcf()

plt.plot(thetas,plotpts[:,0],color='b')
plt.plot(thetas,plotpts[:,1],color='g')
for i in range(0,numsteps):
    if (markers[i] == 'x'):
        plt.scatter(thetas[i],plotpts[i,1],color='k',marker=markers[i])

plt.xlabel('angle (0 to 2pi)')
plt.legend(['fx1', 'fy1'])
plt.title('left finger forces (linear mapping); x marks friction fails')






"""
PART 2: Nonlinear Mapping from fbody to finger forces
what Yoshikawa and Nagai suggest, and what we used on our two-finger hand (Grffin)
is a different definition of internal grasp force. Let the internal force
be based on min{fx1,-fx2}. Which one is less will depend on the external force.
"""
#Use same fint, fext as before

#Internal force based only on the left finger
Gleft = np.array([[1,0,1,0],[0,1,0,1],[0,-1,0,1],[1,0,0,0]])
Gleftinv = np.linalg.inv(Gleft)
#Internal force based only on the right finger
Gright = np.array([[1,0,1,0],[0,1,0,1],[0,-1,0,1],[0,0,-1,0]])
Grightinv = np.linalg.inv(Gright)

#reset markers list
markers = ['.'] * numsteps

for i in range(0,numsteps):
    fbody[0] = fext*np.cos(thetas[i])
    fbody[1] = fext*np.sin(thetas[i])
    if (fbody[0] > 0):
        ffingers = np.dot(Grightinv,fbody)
    else:
        ffingers = np.dot(Gleftinv,fbody)
    
    plotpts[i,0] = ffingers[0]
    plotpts[i,1] = ffingers[1]
#Highlight any cases where friction constraint is not satisfied
    if (np.abs(ffingers[1])/ffingers[0] > 0.5):
        markers[i] = 'x' 


#Plot the forces on the left finger as a function of theta
plt.figure(2)
fig2 = plt.gcf()
plt.plot(thetas,plotpts[:,0],color='b')
plt.plot(thetas,plotpts[:,1],color='g')
for i in range(0,numsteps):
    if (markers[i] == 'x'):
        plt.scatter(thetas[i],plotpts[i,1],color='k',marker=markers[i])

plt.xlabel('angle (0 to 2pi)')
plt.legend(['fx1', 'fy1'])
plt.title('left finger forces (nonlinear mapping); x marks friction fails')