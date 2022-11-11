include("jsat_mtk.jl")

#Make the Rigid Body
rb = make(RigidBody(1000*I(3),[pi/4,0,0],zeros(3)))

#Make the Reaction Wheels
r1 = ReactionWheel(0.25, 1, [ 1, 0, 0], 20)
r2 = ReactionWheel(0.25, 1, [ 0, 1, 0], 20)
r3 = ReactionWheel(0.25, 1, [ 0, 0, 1], 20)
rw = make([r1,r2,r3])

#Make the controller
@named c = controller(samplerate = 0.1,kd=0,kp=1)

#Fake the Thrusters
@named th = step(times = [Inf,Inf,Inf], durations = ones(3),values = ones(3))

eqs = [    
    c.u ~ rw.u
    th.u ~ rb.Te
    rw.T ~ rb.Ti
    rw.H ~ rb.Hi
    rb.θ ~ c.θ
    rb.ω ~ c.ω
]
sys = compose(ODESystem(eqs,t;name=:sys),rb,rw,c,th)
#sys_ = structural_simplify(sys)
#sol = test(sys_)
