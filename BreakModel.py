from modsim import *

state = State(cx=0,cy=0,cvx=1,cvy=0,b1x=1,b1y=0,b1vx=0,b1vy=0)
system = System(f_force_mag=0.0166698, m=0.1701, r=0.028575, balls=1)
params = Params()

def f_force(system, v_vector):
    unpack(system)
    f_force_vector = v_vector.hat() * -f_force_mag
    return f_force_vector

def dist(v1, v2):
    dist = v1.dist(v2)
    return dist

def newv1(p1, p2, v1, v2):
    v1diff = v1 - v2
    x1diff = p1 - p2
    v1prime = v1 - x1diff * v1diff.dot(x1diff) / (x1diff.mag**2)
    return v1prime

def newv2(p1, p2, v1, v2):
    v2diff = v2 - v1
    x2diff = p2 - p1
    v2prime = v2 - x2diff * v2diff.dot(x2diff) / (x2diff.mag**2)
    return v2prime

def slope_func(state, t, system):
    unpack(system)
    ids = linrange(0,balls,1,endpoint=True)
    pos_vectors = []
    vel_vectors = []
    dvs = []
    collisions = []
    for b in ids:
        pos_vectors.append(Vector(state[4*b],state[4*b+1]))
        vel_vectors.append(Vector(state[4*b+2],state[4*b+3]))
    checkdist = linrange(0,balls,1)
    for b1 in checkdist:
        if b1 == 0:
            collisions.append(0)
        rest = linrange(b1+1,balls,1,endpoint=True)
        for b2 in rest:
            if b1 == 0:
                collisions.append(0)
            if dist(pos_vectors[b1],pos_vectors[b2]) <= 2*r:
                collisions[b1] = collisions[b1] + 1
                collisions[b2] = collisions[b2] + 1


    for b1 in checkdist:
        rest = linrange(b1+1,balls,1,endpoint=True)
        for b2 in rest:
            if dist(pos_vectors[b1],pos_vectors[b2]) <= 2*r:
                if collisions[b1] > 1 or collisions[b2] > 1:
                    if collisions[b1] > 2 or collisions[b2] > 2:
                        print('well fuck')
                    # do some shit here for when threee balllllls touch eachother
                v1prime = newv1(pos_vectors[b1],pos_vectors[b2],vel_vectors[b1],vel_vectors[b2])
                v2prime = newv2(pos_vectors[b1],pos_vectors[b2],vel_vectors[b1],vel_vectors[b2])
                vel_vectors[b1] = v1prime
                vel_vectors[b2] = v2prime

    for b in ids:
        dvb = vel_vectors[b] * f_force(system, vel_vectors[b])/m
        dvs.append(dvb)
        state[4*b+2] = dvb.x
        state[4*b+3] = dvb.y



    print('delta vs')
    print(dvs)
    print('pos')
    print(pos_vectors)
    print('vel')
    print(vel_vectors)

slope_func(state,1,system)
slope_func(state,2,system)
