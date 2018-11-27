from modsim import *

state = State(cx=0,cy=0,cvx=1,cvy=0,b1x=1,b1y=0,b1vx=0,b1vy=0)
system = System()
params = Params()

def slope_func(state, t, system):
    ids = linrange(0,1,1,endpoint=True)
    pos_vectors = []
    vel_vectors = []
    for b in ids:
        pos_vectors.append(Vector(state[4*b],state[4*b+1]))
        vel_vectors.append(Vector(state[4*b+2],state[4*b+3]))

slope_func(state,1,system)
