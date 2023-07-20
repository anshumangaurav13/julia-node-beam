rk1(t, dt, y, ẏ) = y + dt * ẏ(t, y)

function rk1!(t, dt, y, ẏ)
    y .+= dt * ẏ(t, y)

    return
end

euler(t, dt, y, ẏ) = rk1(t, dt, y, ẏ)
euler!(t, dt, y, ẏ) = rk1!(t, dt, y, ẏ)

euler_cromer(t, dt, y, ẏ) = y + dt * ẏ(t, [y[1] + dt * y[2], y[2]])

function euler_cromer!(t, dt, y, ẏ)
    y .+= dt * ẏ(t, [y[1] + dt * y[2], y[2]])

    return
end

rk21(t, dt, y, ẏ) = y + dt * ẏ(t + 0.5dt, y + 0.5dt * ẏ(t, y))

function rk21!(t, dt, y, ẏ)
    y .+= dt * ẏ(t + 0.5dt, y + 0.5dt * ẏ(t, y))

    return
end

function rk3(t, dt, y, ẏ)
    k1 = dt * ẏ(t, y)
    k2 = dt * ẏ(t + 0.5dt, y + 0.5k1)
    k_ = dt * ẏ(t + dt, y + k1)
    k3 = dt * ẏ(t + dt, y + k_)
    y + (k1 + 4k2 + k3) / 6
end

function rk3!(t, dt, y, ẏ)
    k1 = dt * ẏ(t, y)
    k2 = dt * ẏ(t + 0.5dt, y + 0.5k1)
    k_ = dt * ẏ(t + dt, y + k1)
    k3 = dt * ẏ(t + dt, y + k_)
    y .+= (k1 + 4k2 + k3) / 6

    return
end

function rk4(t, dt, y, ẏ)
    k1 = dt * ẏ(t, y)
    k2 = dt * ẏ(t + 0.5dt, y + 0.5k1)
    k3 = dt * ẏ(t + 0.5dt, y + 0.5k2)
    k4 = dt * ẏ(t + dt, y + k3)

    y + (k1 + 2k2 + 2k3 + k4) * (1 / 6)
end

function rk4!(t, dt, y, ẏ)
    k1 = dt * ẏ(t, y)
    k2 = dt * ẏ(t + 0.5dt, y + 0.5k1)
    k3 = dt * ẏ(t + 0.5dt, y + 0.5k2)
    k4 = dt * ẏ(t + dt, y + k3)
    y .+= (k1 + 2k2 + 2k3 + k4) * (1 / 6)

    return
end

function verlet(t, dt, y, ẏ)
    (x, v) = (y[1], y[2])

    x .+= 0.5dt * v
    v .+= dt * ẏ(t, [x, v])[2]
    x .+= 0.5dt * v

    [x, v]
end

function verlet!(t, dt, y, ẏ)
    y[1] .+= 0.5dt * y[2]
    y[2] .+= dt * ẏ(t, y)[2]
    y[1] .+= 0.5dt * y[2]

    return
end

const θ = 1 / (2 - ∛2)
function forest_ruth(t, dt, y, ẏ)
    (x, v) = (y[1], y[2])

    x .+= θ * 0.5dt * v
    v .+= θ * dt * ẏ(t, [x, v])[2]
    x .+= (1 - θ) * 0.5dt * v
    v .+= (1 - 2θ) * dt * ẏ(t, [x, v])[2]
    x .+= (1 - θ) * 0.5dt * v
    v .+= θ * dt * ẏ(t, [x, v])[2]
    x .+= θ * 0.5dt * v

    [x, v]
end

function forest_ruth!(t, dt, y, ẏ)
    y[1] .+= θ * 0.5dt * y[2]
    y[2] .+= θ * dt * ẏ(t, y)[2]
    y[1] .+= (1 - θ) * 0.5dt * y[2]
    y[2] .+= (1 - 2θ) * dt * ẏ(t, y)[2]
    y[1] .+= (1 - θ) * 0.5dt * y[2]
    y[2] .+= θ * dt * ẏ(t, y)[2]
    y[1] .+= θ * 0.5dt * y[2]

    return
end

const ξ = +0.17861789584480910
const λ = -0.21234183106260540
const χ = -0.06626458266981849
function PEFRL(t, dt, y, ẏ)
    (x, v) = (y[1], y[2])

    x .+= ξ * dt * v
    v .+= (1 - 2λ) * 0.5dt * ẏ(t, [x, v])[2]
    x .+= χ * dt * v
    v .+= λ * dt * ẏ(t, [x, v])[2]
    x .+= (1 - 2(χ + ξ)) * dt * v
    v .+= λ * dt * ẏ(t, [x, v])[2]
    x .+= χ * dt * v
    v .+= (1 - 2λ) * 0.5dt * ẏ(t, [x, v])[2]
    x .+= ξ * dt * v

    [x, v]
end

function PEFRL!(t, dt, y, ẏ)
    y[1] .+= ξ * dt * y[2]
    y[2] .+= (1 - 2λ) * 0.5dt * ẏ(t, y)[2]
    y[1] .+= χ * dt * y[2]
    y[2] .+= λ * dt * ẏ(t, y)[2]
    y[1] .+= (1 - 2(χ + ξ)) * dt * y[2]
    y[2] .+= λ * dt * ẏ(t, y)[2]
    y[1] .+= χ * dt * y[2]
    y[2] .+= (1 - 2λ) * 0.5dt * ẏ(t, y)[2]
    y[1] .+= ξ * dt * y[2]

    return
end

function ode_step(t, DT, y, ẏ; ode_solver::Function, steps::Int64)
    if steps == 1
        return ode_solver(t, DT, y, ẏ)
    end
    dt = DT / steps
    y_tmp = y
    for tt in t:dt:t+DT-dt
        y_tmp = ode_solver(tt, dt, y_tmp, ẏ)
    end
    y_tmp
end

function ode_step!(t, DT, y, ẏ; ode_solver::Function, steps::Int64)
    if steps == 1
        ode_solver(t, DT, y, ẏ)
        return
    end
    dt = DT / steps
    for tt in t:dt:t+DT-dt
        ode_solver(tt, dt, y, ẏ)
    end
    
    return
end