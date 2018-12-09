from modsim import *
p0=1.5
h=0.0581536058641
rp=0.033575
init = State(cx=0.5,cy=0,cvx=1,cvy=0,
             b1x=1.5,b1y=0,b1vx=0,b1vy=0,
             b2x=p0+h,b2y=rp,b2vx=0,b2vy=0,
             b3x=p0+h,b3y=-rp,b3vx=0,b3vy=0,
             b4x=p0+2*h,b4y=0,b4vx=0,b4vy=0,
             b5x=p0+2*h,b5y=2*rp,b5vx=0,b5vy=0,
             b6x=p0+2*h,b6y=-2*rp,b6vx=0,b6vy=0,
             b7x=p0+3*h,b7y=rp,b7vx=0,b7vy=0,
             b8x=p0+3*h,b8y=3*rp,b8vx=0,b8vy=0,
             b9x=p0+3*h,b9y=-rp,b9vx=0,b9vy=0,
             b10x=p0+3*h,b10y=-3*rp,b10vx=0,b10vy=0,
             b11x=p0+4*h,b11y=0,b11vx=0,b11vy=0,
             b12x=p0+4*h,b12y=2*rp,b12vx=0,b12vy=0,
             b13x=p0+4*h,b13y=4*rp,b13vx=0,b13vy=0,
             b14x=p0+4*h,b14y=-2*rp,b14vx=0,b14vy=0,
             b15x=p0+4*h,b15y=-4*rp,b15vx=0,b15vy=0)
state = init
system = System(init = init, fs=0.0166698, fk=0.333396, m=0.1701, r=0.028575, balls = 15, dist_step=0.5,  t_0=0, t_end=7, minvel=.01, h=1, w=2, t_cut=1)

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
                dist_error = 100 * (2*r - pos_vectors[b1].dist(pos_vectors[b2])) / (2*r)
                print(str(b1) + ' and ' + str(b2) + ' collided --- ' + 'dist error ' + str(dist_error))
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
            if b == 0:
                system.score = -15
            if b != 0:
                system.score += 1
            ps[b] = Vector(2+b,2)
            vs[b] = Vector(0,0)

    results = State(cx=ps[0].x,cy=ps[0].y,cvx=vs[0].x,cvy=vs[0].y,
                    b1x=ps[1].x,b1y=ps[1].y,b1vx=vs[1].x,b1vy=vs[1].y,
                    b2x=ps[2].x,b2y=ps[2].y,b2vx=vs[2].x,b2vy=vs[2].y,
                    b3x=ps[3].x,b3y=ps[3].y,b3vx=vs[3].x,b3vy=vs[3].y,
                    b4x=ps[4].x,b4y=ps[4].y,b4vx=vs[4].x,b4vy=vs[4].y,
                    b5x=ps[5].x,b5y=ps[5].y,b5vx=vs[5].x,b5vy=vs[5].y,
                    b6x=ps[6].x,b6y=ps[6].y,b6vx=vs[6].x,b6vy=vs[6].y,
                    b7x=ps[7].x,b7y=ps[7].y,b7vx=vs[7].x,b7vy=vs[7].y,
                    b8x=ps[8].x,b8y=ps[8].y,b8vx=vs[8].x,b8vy=vs[8].y,
                    b9x=ps[9].x,b9y=ps[9].y,b9vx=vs[9].x,b9vy=vs[9].y,
                    b10x=ps[10].x,b10y=ps[10].y,b10vx=vs[10].x,b10vy=vs[10].y,
                    b11x=ps[11].x,b11y=ps[11].y,b11vx=vs[11].x,b11vy=vs[11].y,
                    b12x=ps[12].x,b12y=ps[12].y,b12vx=vs[12].x,b12vy=vs[12].y,
                    b13x=ps[13].x,b13y=ps[13].y,b13vx=vs[13].x,b13vy=vs[13].y,
                    b14x=ps[14].x,b14y=ps[14].y,b14vx=vs[14].x,b14vy=vs[14].y,
                    b15x=ps[15].x,b15y=ps[15].y,b15vx=vs[15].x,b15vy=vs[15].y)

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
            plot(results.b3x, results.b3y)
            plot(results.b4x, results.b4y)
            plot(results.b5x, results.b5y)
            plot(results.b6x, results.b6y)
            plot(results.b7x, results.b7y)
            plot(results.b8x, results.b8y)
            plot(results.b9x, results.b9y)
            plot(results.b10x, results.b10y)
            plot(results.b11x, results.b11y)
            plot(results.b12x, results.b12y)
            plot(results.b13x, results.b13y)
            plot(results.b14x, results.b14y)
            plot(results.b15x, results.b15y)
            savefig('graphs/balls----y'+str(ystart)+'----v'+str(vstart) + '--.png')
            outcome.set_value(ystart, vstart, score)
            print(outcome)

    return outcome

sweep_y = linspace(0,0.05715,8, endpoint=True)
sweep_v = linspace(1,10,8, endpoint=True)

print(sweep_y)
print(sweep_v)

results = sweep_sim(init,state,system,update_func,sweep_y,sweep_v)

print(results)
