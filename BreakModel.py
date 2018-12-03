from modsim import *

init = State(cx=0,cy=0,cvx=1,cvy=0,b1x=1,b1y=.02,b1vx=0,b1vy=0)
state = init
system = System(init=init, f_force_mag=0.0166698, m=0.1701, r=0.028575, balls=1, dt=.05, t_0=0, t_end=20, minvel=.01)

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

def update_func(state, t, system):
    unpack(system)
    ids = linrange(0,balls,1,endpoint=True)
    pos_vectors = []
    vel_vectors = []
    ps = []
    vs = []
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
                    print('three ball collision')
                print(str(b1) + ' and ' + str(b2) + ' collided')
                v1prime = newv1(pos_vectors[b1],pos_vectors[b2],vel_vectors[b1],vel_vectors[b2])
                v2prime = newv2(pos_vectors[b1],pos_vectors[b2],vel_vectors[b1],vel_vectors[b2])
                vel_vectors[b1] = v1prime
                vel_vectors[b2] = v2prime

    for b in ids:
        pb = pos_vectors[b] + vel_vectors[b] * dt
        vb = vel_vectors[b] + (f_force(system, vel_vectors[b])/m) * dt
        if vb.mag <= minvel:
            vb = Vector(0,0)
        vs.append(vb)
        ps.append(pb)
        print('ps' + str(b) + '--' + str(t+dt))
        print(ps)
        print('vs' + str(b) + '--' + str(t+dt))
        print(vs)

    results = State(cx=ps[0].x,cy=ps[0].y,cvx=vs[0].x,cvy=vs[0].y,b1x=ps[1].x,b1y=ps[1].y,b1vx=vs[1].x,b1vy=vs[1].y)
    return results

def run_sim(system, update_func):
    unpack(system)

    frame = TimeFrame(columns=init.index)
    frame.loc[t_0] = init
    ts = linrange(t_0,t_end,dt)

    for t in ts:
        frame.row[t+dt] = update_func(frame.row[t], t, system)

    return frame

results = run_sim(system, update_func)

plot(results.cx, results.cy)
plot(results.b1x, results.b1y)
savefig('graphs\cue.png')
