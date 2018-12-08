from modsim import *

init = State(cx=0.5,cy=0,cvx=3.9,cvy=0, b1x=1.5,b1y=.02,b1vx=0,b1vy=0, b2x=1.55,b2y=-0.01,b2vx=0,b2vy=0)
state = init
system = System(init = init, fs=0.0166698, fk=0.333396, m=0.1701, r=0.028575, balls=2, dist_step=2, t_0=0, t_end=10, minvel=.01, h=1, w=2, t_cut=1)

def f_force(system, v_vector):
    unpack(system)
    if v_vector.mag > 0.5:
        f_force_vector = v_vector.hat() * -fk
        return f_force_vector
    f_force_vector = v_vector.hat() * -fs
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
    tlc = Vector(0,h/2)
    blc = Vector(0,-h/2)
    trc = Vector(w,h/2)
    brc = Vector(w,-h/2)

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

    for b in ids:
        bld = pos_vectors[b].x
        brd = 2 - pos_vectors[b].x
        bbd = 0.5 + pos_vectors[b].y
        btd = 0.5 - pos_vectors[b].y
        if bld <= r or brd <= r:
            vel_vectors[b] = Vector(-vel_vectors[b].x, vel_vectors[b].y)
        if bbd <= r or btd <=r:
            vel_vectors[b] = Vector(vel_vectors[b].x, -vel_vectors[b].y)

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
                print('1position-1: ' + str(pos_vectors[b1]))
                print('1position-2: ' + str(pos_vectors[b2]))
                dist_error = r - pos_vectors[b1].dist(pos_vectors[b2])
                print('dist error ' + str(dist_error))
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

    for b in ids:
        if pos_vectors[b].dist(tlc) <= 2*r or pos_vectors[b].dist(blc) <= 2*r or pos_vectors[b].dist(trc) <= 2*r or pos_vectors[b].dist(brc) <= 2*r:
            print('ball sunk')
            system.score += 1
            ps[b] = Vector(2+b,2)
            vs[b] = Vector(0,0)

    results = State(cx=ps[0].x,cy=ps[0].y,cvx=vs[0].x,cvy=vs[0].y, b1x=ps[1].x,b1y=ps[1].y,b1vx=vs[1].x,b1vy=vs[1].y, b2x=ps[2].x,b2y=ps[2].y,b2vx=vs[2].x,b2vy=vs[2].y)

    return results

def run_sim(system, update_func):
    system = System(system, dt = system.dist_step / (100* sqrt(state.cvx**2 + state.cvy**2)), score=0)
    unpack(system)

    frame = TimeFrame(columns=init.index)
    frame.loc[t_0] = init
    ts = linrange(t_0, t_end, dt)

    for t in ts:
        frame.row[t+dt] = update_func(frame.row[t], t, system)

    return frame, score

def sweep_sim(init,state,system,update_func, sweep_y, sweep_v):
    unpack(system)
    columnlabels = ['y','v','score']
    outcome = DataFrame(index=sweep_y, columns=sweep_v)
    for ystart in sweep_y:
        print('ystart')
        print(ystart)
        init.cy = ystart

        for vstart in sweep_v:
            init.cvx = vstart
            state = init
            system.init = init
            print('vstart')
            print(vstart)
            results, score = run_sim(system, update_func)
            plt.figure(10*vstart+ystart)
            plot(results.cx, results.cy)
            plot(results.b1x, results.b1y)
            plot(results.b2x, results.b2y)
            savefig('graphs/balls----y'+str(ystart)+'----v'+str(vstart) + '--.png')
            outcome.set_value(ystart, vstart, score)

    return outcome

sweep_y = linspace(-.4,.4,3, endpoint=True)
sweep_v = linspace(1,20,3, endpoint=True)

print(sweep_y)
print(sweep_v)

results = sweep_sim(init,state,system,update_func,sweep_y,sweep_v)

print(results)
