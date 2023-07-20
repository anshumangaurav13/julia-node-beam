include("nodebeam.jl")

function main(ode, steps=500)
    l = 10.0
    m = 5
    k = 50.0
    damp = 10.0

    n1 = Node(m, [0, 0, 0], true)
    n2 = Node(m, [0, l, 0], true)
    n3 = Node(m, [l, 0, 0])
    n4 = Node(m, [2l, 0, 0])
    n5 = Node(m, [l, l, 0])
    n6 = Node(m, [2l, l, 0])
    n7 = Node(m, [(2.1 + √3 / 2)l, 0.5l, 0])
    nodes = [n1, n2, n3, n4, n5, n6, n7]

    b1 = (1, 3, k, l, damp)
    b2 = (2, 5, k, l, damp)
    b3 = (5, 3, k, l, damp)
    b4 = (3, 4, k, l, damp)
    b5 = (5, 6, k, l, damp)
    b6 = (4, 6, k, l, damp)
    b7 = (7, 6, k, l, damp)
    b8 = (4, 7, k, l, damp)
    b9 = (2, 3, 10k, √2l, damp)
    b10 = (4, 5, 10k, √2l, damp)
    b11 = (1, 5, 10k, √2l, damp)
    b12 = (3, 6, 10k, √2l, damp)
    beams = [b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12]

    extforce(t, sim) = [0, -9.81, 0]

    sim = Sim(nodes, beams, extforce)

    run(sim, [0, 10.0], 1 / 30.0; ode_solver=ode, substeps=steps)
end

function doublepend(ode, steps=500)
    l = 200.0
    m = 50.0
    k = 55000.0
    damp = 5.0

    n1 = Node(m, [0, 0, 0], true)
    n2 = Node(2m, [l, 0, 0])
    n3 = Node(0.5m, [2l, 0, 0])
    nodes = [n1, n2, n3]

    b1 = (1, 2, k, l, damp)
    b2 = (2, 3, k, l, damp)
    beams = [b1, b2]

    extforce(t, sim) = [0, -981, 0]

    sim = Sim(nodes, beams, extforce)

    run(sim, [0, 30.0], 1 / 60.0; ode_solver=ode, substeps=steps)
end

using Plots
function save_frames(sim::Sim, height)
    phases = sim.history[2]

    poss = [phase[1] for phase in phases]
    plot()
    return @animate for i in 1:length(phases)
        p = plot(
            aspect_ratio=:equal,
            size=(height, height)
        )
        for beam in sim.beams
            n1 = beam[1]
            n2 = beam[2]
            xpos = [poss[i][n1][1], poss[i][n2][1]]
            ypos = [poss[i][n1][2], poss[i][n2][2]]
            plot!(p, xpos, ypos,
                label="$(beam[3]), $(beam[4]), $(beam[5])",
                linewidth=map(beam[3], 0.5, 4)
            )
        end
        for j in 1:sim.N_nodes
            possj = [possn[j] for possn in poss]
            xpos = [p[1] for p in possj]
            ypos = [p[2] for p in possj]
            scatter!(p, [xpos[i]], [ypos[i]],
                markersize=map(sim.masses[j], 5, 10),
                label="$(sim.masses[j])"
            )
            plot!(p, xpos[1:i], ypos[1:i], label="")
        end
    end
end

map(x, a, b) = (0.5(x - 10) / (1 + abs(x - 10)) + 0.5) * (b - a) + a