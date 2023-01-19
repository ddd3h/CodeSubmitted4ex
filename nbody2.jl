using Distributions
using LinearAlgebra
using PyPlot

const G = 1
const N = 1000
const M = 20/N
const dt = 0.01
const max_steps = 100
const max = 2
const min = - max

function f(m,r1,r2)
    return - G * m .* (r1 - r2) ./(sqrt((norm(r1 - r2))^2 + 0.01)^3)
end


function initialize()
    # Positions
    r0 = [rand(Normal(0,1),2) for i in 1:N]

    # Velocities
    v0 = [rand(Normal(0,1),2) for i in 1:N]
    for i in 1:N
        v0[i] .-= mean(v0[i])
    end
    
    return r0, v0
end


function calc(r)
    xmax = maximum([abs(r[i][1]) for i in 1:N])
    ymax = maximum([abs(r[i][2]) for i in 1:N])
    xymax = maximum([xmax, ymax])
    xymin = - xymax
    s = xymax - xymin
    xmax = xymax
    ymax = xymax
    xmin = xymin
    ymin = xymin
    treelist = []
    d = [i for i in 1:N]
    push!(treelist,[mean([r[i][1] for i in 1:N]), mean([r[i][2] for i in 1:N]), N, 0, d, xmin, xmax, ymin, ymax])
    _treelist = segmentation(xmin,xmax,ymin,ymax,r,d,0,treelist)
    cat(treelist, _treelist, dims=1)
    sort!(treelist, by= x -> x[4])

    l = [zeros(2) for i in 1:N]
    for i in 1:N
        caculated::Vector{Int64} = []
        t = zeros(2)
        for j in 1:N
            if i != j
                if !(j in caculated)
                    w = r[i] - treelist[j][1:2]
                    s = treelist[j][7] - treelist[j][6]
                    theta = s/norm(w)
                    if theta < 10
                        cat(caculated, treelist[j][5], dims=1)
                        t .+= f(treelist[j][3]*M, r[i], treelist[j][1:2])
                    else
                        continue
                    end
                end
            end
        end
        l[i] = t
    end
    return l
end

function segmentation(xmin,xmax,ymin,ymax,r,d,depth,treelist)
    n = length(r)
    s = xmax - xmin
    NSEWcount::Vector{Int64} = [0,0,0,0]
    NSEWparticlelist::Vector{Vector{Int64}} = [[],[],[],[]]
    NSEWposition_x = [0.0 for i in 1:4]
    NSEWposition_y = [0.0 for i in 1:4]
    xmaxlist = [xmax, xmin + s/2, xmin + s/2, xmax]
    xminlist = [xmax - s/2, xmin, xmin, xmax - s/2]
    ymaxlist = [ymax, ymax, ymin + s/2, ymin + s/2]
    yminlist = [ymax - s/2, ymax - s/2, ymin, ymin]

    for i in d
        if (xminlist[1] <= r[i][1] <= xmaxlist[1]) && (yminlist[1] <= r[i][2] <= ymaxlist[1])
            areanum = 1
        elseif (xminlist[2] <= r[i][1] <= xmaxlist[2]) && (yminlist[2] <= r[i][2] <= ymaxlist[2])
            areanum = 2
        elseif (xminlist[3] <= r[i][1] <= xmaxlist[3]) && (yminlist[3] <= r[i][2] <= ymaxlist[3])
            areanum = 3
        elseif (xminlist[4] <= r[i][1] <= xmaxlist[4]) && (yminlist[4] <= r[i][2] <= ymaxlist[4])
            areanum = 4
        end
        
        push!(NSEWparticlelist[areanum], i)
        NSEWcount[areanum] += 1
        NSEWposition_x[areanum] += r[i][1]
        NSEWposition_y[areanum] += r[i][2]
        
    end

    for i in 1:4
        if NSEWcount[i] == 1
            push!(treelist,[NSEWposition_x[i]/NSEWcount[i], NSEWposition_y[i]/NSEWcount[i], NSEWcount[i], depth + 1, NSEWparticlelist[i], xminlist[i], xmaxlist[i], yminlist[i], ymaxlist[i]])
        elseif NSEWcount[i] > 1
            push!(treelist,[NSEWposition_x[i]/NSEWcount[i], NSEWposition_y[i]/NSEWcount[i], NSEWcount[i], depth + 1, NSEWparticlelist[i], xminlist[i], xmaxlist[i], yminlist[i], ymaxlist[i]])
            segmentation(xminlist[i], xmaxlist[i], yminlist[i], ymaxlist[i], r, NSEWparticlelist[i],depth+1,treelist)
        end
    end
    return treelist
end

function initialize2()
    # Positions
    r0 = [rand(Uniform(-2,2),2) for i in 1:N]

    # Velocities
    v0 = [zeros(2) for i in 1:N]
    for i in 1:N
        v0[i] .-= mean(v0[i])
    end
    
    return r0, v0
end

function update(r0, v0)
    k1 = copy(v0)
    l1 = calc(r0)
    
    k2 = [v0[i] + dt/2 * l1[i] for i in 1:N]
    l2 = calc(r0 .+ dt/2 * k1)

    k3 = [v0[i] + dt/2 * l2[i] for i in 1:N]
    l3 = calc(r0 .+ dt/2 * k2)

    k4 = [v0[i] + dt * l3[i] for i in 1:N]
    l4 = calc(r0 .+ dt * k3)

    k = [(k1[i] + 2 .* k2[i] + 2 .* k3[i] + k4[i]) ./ 6 for i in 1:N]
    l = [(l1[i] + 2 .* l2[i] + 2 .* l3[i] + l4[i]) ./ 6 for i in 1:N]

    r = [r0[i] + dt * k[i] for i in 1:N]
    v = [v0[i] + dt * l[i] for i in 1:N]

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
    r, v = initialize2()
    # pplot(r, 0)
    for step in 1:max_steps
        println(step)
        r, v = update(r, v)
        # pplot(r, step)
    end
end

main()

