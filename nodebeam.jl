include("ode_solvers.jl")

struct Node
    mass::Float64
    pos::Vector{Float64}
    vel::Vector{Float64}
    fixed::Bool
    Node(m, p) = new(m, p, [0.0, 0, 0], false)
    Node(m, p, f::Bool) = new(m, p, [0.0, 0, 0], f)
end

mutable struct Sim
    N_nodes::Int64
    N_fixed::Int64
    masses::Vector{Float64}
    fixed::Vector{Bool}
    phase
    beams
    extforce::Function
    history
    function Sim(nodes, beams, extforce)
        nn = length(nodes)
        ms = [n.mass for n in nodes]
        fs = [n.fixed for n in nodes]
        nf = sum(fs)
        return new(
            nn,
            nf,
            ms,
            fs,
            [[n.pos for n in nodes], [n.vel for n in nodes]],
            beams,
            extforce,
        )
    end
end

function stiff_matrix(s::Sim)
    m = zeros(s.N_nodes, s.N_nodes)
    for beam in s.beams
        i = beam[1]
        j = beam[2]
        k = beam[3]
        if !(s.fixed[i] && s.fixed[j])
            m[i, i] += k
            m[i, j] -= k
            m[j, i] -= k
            m[j, j] += k
        else
            if s.fixed[i]
                m[j, j] += k
            elseif s.fixed[j]
                m[i, i] += k
            end
        end
    end
    return m
end

function damp_matrix(s::Sim)
    m = zeros(s.N_nodes, s.N_nodes)
    for beam in s.beams
        i = beam[1]
        j = beam[2]
        k = beam[5]
        if !(s.fixed[i] && s.fixed[j])
            m[i, i] += k
            m[i, j] -= k
            m[j, i] -= k
            m[j, j] += k
        else
            if s.fixed[i]
                m[j, j] += k
            elseif s.fixed[j]
                m[i, i] += k
            end
        end
    end
    return m
end

mass_matrix(s::Sim) = diagm([m for m in s.masses])

normsq(v) = sum(v .* v)
norm(v) = sqrt(normsq(v))
dir(v) = v / norm(v)
dot(u, v) = sum(u .* v)

function calc_phase_diff(t, phase, sim::Sim)
    acc = 0.0(phase[2])

    for (i, j, ks, l, kd) in sim.beams
        rij = phase[1][j] - phase[1][i]
        vij = phase[2][j] - phase[2][i]
        dir_rij = dir(rij)

        F = ((norm(rij) - l) * ks + dot(vij, dir_rij) * kd) * dir_rij

        if sim.fixed[i] == false
            acc[i] .+= F / sim.masses[i]
        end
        if sim.fixed[j] == false
            acc[j] .-= F / sim.masses[j]
        end
    end

    ext = sim.extforce(t, sim)
    for i in 1:sim.N_nodes
        if sim.fixed[i] == false
            acc[i] .+= ext
        end
    end

    [phase[2], acc]
end

function run(sim::Sim, tspan, dt::Float64; ode_solver::Function, substeps::Int64)
    ts = tspan[1]:dt:tspan[2]
    N_steps = length(ts)
    ys = [0.0sim.phase for i in 1:N_steps]
    ys[1] .= 1sim.phase

    ẏ(t, y) = calc_phase_diff(t, y, sim)

    for i in 2:N_steps
        ode_step!(ts[i], dt, sim.phase, ẏ; ode_solver=ode_solver, steps=substeps)
        ys[i] .= 1sim.phase
    end

    sim.history = [ts, ys]
    return sim
end

