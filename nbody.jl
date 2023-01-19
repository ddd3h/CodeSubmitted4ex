using Distributions
using LinearAlgebra
using PyPlot

const G::Float64 = 1.0
const N::Int64 = 2000
const m::Float64 = 20/N
const dt::Float64 = 0.01
const max_steps::Int64 = 1000
const max = 2
const min = - max

function f(r1::Vector{Float64},r2::Vector{Float64})
    return - G * m .* (r1 - r2) ./(sqrt((norm(r1 - r2))^2 + 0.01)^3)
end


function initialize()
    # Positions
    r0 = [rand(Normal(0,1),2) for i in 1:N]

    # Velocities
    v0 = [rand(Normal(0,1),2) for i in 1:N]
    a = mean(v0)
    for i in 1:N
        v0[i] .-= a
    end
    
    return r0, v0
end

function initialize2()
    # Positions
    r0 = [rand(Uniform(-2,2),2) for i in 1:N]

    # Velocities
    v0 = [zeros(2) for i in 1:N]
    
    return r0, v0
end

function initialize3()
    # Positions
    a = Int(floor(N/2))
    r0 = [rand(Normal(0,0.3),2) for i in 1:a]
    
    for i in 1:a
        r0[i][1] = r0[i][1] - 1
    end
    b = [rand(Normal(0,0.3),2) for i in a+1:N]
    for i in 1:(N-(a+1)+1)
        b[i][1] = b[i][1] + 1
    end
    r0 = cat(r0,b,dims=1)
    
    # Velocities
    c = [rand(Uniform(0,2)) for i in 1:N]
    v0 = [zeros(2) for i in 1:N]
    for i in 1:N
        v0[i][1] = - c[i]*r0[i][2]
        v0[i][2] = c[i]*r0[i][1]
    end
    
    return r0, v0
end

function update(r0, v0)
    k1 = copy(v0)
    l1::Vector{Vector{Float64}} = [zeros(2) for i in 1:N]
    for i in 1:N
        t::Vector{Float64} = zeros(2)
        for j in 1:N
            if i != j
                t .+= f(r0[i], r0[j])
            end
        end
        l1[i] = t
    end
    
    k2 = v0 .+ dt/2 .* l1
    l2::Vector{Vector{Float64}} = [zeros(2) for i in 1:N]

    for i in 1:N
        t::Vector{Float64} = zeros(2)
        for j in 1:N
            if i != j
                t .+= f(r0[i] + dt/2 * k1[i], r0[j] + dt/2 * k1[j])
            end
        end
        l2[i] = t
    end

    k3 = v0 .+ dt/2 .* l2
    l3::Vector{Vector{Float64}} = [zeros(2) for i in 1:N]

    for i in 1:N
        t::Vector{Float64} = zeros(2)
        for j in 1:N
            if i != j
                t .+= f(r0[i] + dt/2 * k2[i], r0[j] + dt/2 * k2[j])
            end
        end
        l3[i] = t
    end

    k4 = v0 .+ dt .* l3
    l4::Vector{Vector{Float64}} = [zeros(2) for i in 1:N]

    for i in 1:N
        t::Vector{Float64} = zeros(2)
        for j in 1:N
            if i != j
                t .+= f(r0[i] + dt * k3[i], r0[j] + dt * k3[j])
            end
        end
        l4[i] = t
    end

    k = (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4) ./ 6
    l = (l1 .+ 2 .* l2 .+ 2 .* l3 .+ l4) ./ 6

    r = r0 .+ dt .* k
    v = v0 .+ dt .* l

    return r, v
end


function pplot(r,step)
    x = [r[i][1] for i in 1:N]
    y = [r[i][2] for i in 1:N]
    scatter(x, y, s=1, alpha=0.8)
    xlim(min, max)
    ylim(min, max)
    s = string(step, pad=6)
    savefig("img/step"*s*".png")
    close("all")
end


function main()
    open("e_8.txt","w") do fp
        r::Vector{Vector{Float64}}, v::Vector{Vector{Float64}} = initialize3()
        println(fp, sum([sum(v[i].^2) for i in 1:length(v)]))
        # pplot(r, 0)
        for step in 1:max_steps
            r, v = update(r, v)
            # pplot(r, step)
            println(fp, sum([sum(v[i].^2) for i in 1:length(v)]))
        end
    end
end

main()