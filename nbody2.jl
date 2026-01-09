using Random
using Distributions
using LinearAlgebra
using StaticArrays
using PyPlot
using Printf

# ----------------------------
# 物理・シミュレーション設定
# ----------------------------
const G        = 1.0
const N        = 1000
const M        = 20.0 / N
const dt       = 0.01
const max_steps = 100

# 描画範囲
const plot_max = 2.0
const plot_min = -plot_max

# ソフトニング（発散防止）
const eps2 = 0.01^2

# Barnes–Hut 開閉パラメータ（小さいほど正確・遅い / 大きいほど速い・粗い）
const θ = 0.7

# 2次元ベクトル（確保がほぼ消えて速い）
const Vec2 = SVector{2,Float64}

# ----------------------------
# 初期化
# ----------------------------
# 初期位置分布:
#   デフォルトでは一様分布 Uniform(-2, 2) を用いる（正方形領域内に一様に配置）。
#   以前は Normal(0, 1) を用いており、より銀河形成などに近い
#   ガウス分布を使いたい場合は、たとえば
#       initialize(position_dist = Normal(0, 1))
#   のように明示的に分布を指定してください。
#
# 初期速度:
#   デフォルトでは速度ゼロ (random_velocities=false) を用いる。
#   以前のコードのように Normal(0, 1) から生成し全粒子の平均速度を引く場合は
#       initialize(random_velocities = true)
#   のように指定してください。
function initialize(; seed::Int=0, position_dist::Distribution=Uniform(-2, 2), random_velocities::Bool=false)
    seed != 0 && Random.seed!(seed)
    r0 = Vector{Vec2}(undef, N)
    v0 = Vector{Vec2}(undef, N)

    # 位置は指定された分布から初期化
    @inbounds for i in 1:N
        r0[i] = Vec2(rand(position_dist), rand(position_dist))
    end

    if random_velocities
        # 旧コードと同様: Normal(0,1) から生成し，全粒子の平均速度を引く
        vdist = Normal(0.0, 1.0)
        tmp_v = Vector{Vec2}(undef, N)
        v_sum = Vec2(0.0, 0.0)
        @inbounds for i in 1:N
            vx = rand(vdist)
            vy = rand(vdist)
            vi = Vec2(vx, vy)
            tmp_v[i] = vi
            v_sum += vi
        end
        v_mean = v_sum / N
        @inbounds for i in 1:N
            v0[i] = tmp_v[i] - v_mean
        end
    else
        # デフォルト: 速度をゼロで初期化
        @inbounds for i in 1:N
            v0[i] = Vec2(0.0, 0.0)
        end
    end

    return r0, v0
end

# ----------------------------
# Barnes–Hut 用 Quadtree ノード
# ----------------------------
mutable struct Node
    cx::Float64
    cy::Float64
    h::Float64                 # half-width (正方形領域)
    mass::Float64
    com::Vec2                  # center of mass
    idx::Int                   # 葉で粒子1個なら粒子番号、そうでなければ0
    child::NTuple{4,Union{Node,Nothing}}  # 1:NW 2:NE 3:SW 4:SE
end

@inline new_node(cx, cy, h) = Node(cx, cy, h, 0.0, Vec2(0.0,0.0), 0, (nothing,nothing,nothing,nothing))

@inline function quadrant(n::Node, p::Vec2)
    east  = p[1] >= n.cx
    north = p[2] >= n.cy
    return north ? (east ? 2 : 1) : (east ? 4 : 3)
end

@inline function child_center(n::Node, q::Int)
    hh = n.h/2
    if q == 1
        return (n.cx - hh, n.cy + hh)   # NW
    elseif q == 2
        return (n.cx + hh, n.cy + hh)   # NE
    elseif q == 3
        return (n.cx - hh, n.cy - hh)   # SW
    else
        return (n.cx + hh, n.cy - hh)   # SE
    end
end

@inline function update_mass_com!(n::Node, m::Float64, p::Vec2)
    newmass = n.mass + m
    # 逐次更新（確保なし）
    n.com = (n.com * n.mass + p * m) / newmass
    n.mass = newmass
end

# 最小セルサイズ（これ以下には分割しない）= 無限再帰を防ぐ
# この値は粒子が完全に同じ位置にある場合や数値誤差で区別できない場合に、
# 無限に細分化されるのを防ぐために設定される。
# 1e-10 はソフトニングパラメータ eps2 = 0.01^2 = 1e-4 よりも十分小さく、
# 物理的に意味のある距離スケールよりも小さい値として選択されている。
const MIN_CELL_SIZE = 1e-10

function insert!(n::Node, i::Int, pos::Vector{Vec2})
    p = pos[i]
    update_mass_com!(n, M, p)

    # 葉が空なら格納
    if n.idx == 0 && all(x -> x === nothing, n.child)
        n.idx = i
        return
    end

    # セルサイズが最小値以下なら分割しない（同位置粒子の無限再帰を防ぐ）
    if n.h < MIN_CELL_SIZE
        # この場合、既存粒子と新粒子が同じセルに共存することになる。
        # idx は既存粒子のインデックスを保持するが、質量・重心は両方の粒子を反映している。
        # Barnes-Hut の力計算では、このノード全体の質量・重心が使われるため、
        # 極めて近接した粒子同士の相互作用は正しく近似される。
        # （個別の粒子として扱われないが、質量中心としては両方が考慮される）
        return
    end

    # 既に粒子が入っていたら分割して押し出す
    if n.idx != 0
        old = n.idx
        n.idx = 0

        qold = quadrant(n, pos[old])
        (ccx, ccy) = child_center(n, qold)
        ch = n.child[qold]
        if ch === nothing
            ch = new_node(ccx, ccy, n.h/2)
            n.child = Base.setindex(n.child, ch, qold)
        end
        insert!(ch, old, pos)
    end

    # 今の粒子も子へ
    q = quadrant(n, p)
    (ccx, ccy) = child_center(n, q)
    ch = n.child[q]
    if ch === nothing
        ch = new_node(ccx, ccy, n.h/2)
        n.child = Base.setindex(n.child, ch, q)
    end
    insert!(ch, i, pos)
end

function build_tree(pos::Vector{Vec2})
    xmin = +Inf; xmax = -Inf
    ymin = +Inf; ymax = -Inf
    @inbounds for p in pos
        x = p[1]; y = p[2]
        xmin = min(xmin, x); xmax = max(xmax, x)
        ymin = min(ymin, y); ymax = max(ymax, y)
    end

    cx = (xmin + xmax)/2
    cy = (ymin + ymax)/2
    h  = max(xmax - xmin, ymax - ymin)/2 + 1e-9   # 0割り回避

    root = new_node(cx, cy, h)
    @inbounds for i in 1:length(pos)
        insert!(root, i, pos)
    end
    return root
end

@inline function pair_force(m::Float64, r::Vec2, c::Vec2)
    dr = r - c
    d2 = dr[1]^2 + dr[2]^2 + eps2
    invd = inv(sqrt(d2))
    invd3 = invd^3
    return -G * m * dr * invd3
end

function force_from_node(i::Int, r::Vec2, n::Node)
    # 自分自身だけの葉なら無視
    if n.idx == i && all(x -> x === nothing, n.child)
        return Vec2(0.0, 0.0)
    end

    # BH判定: s/d < θ ならまとめる（遠いほどまとめる）
    dr = r - n.com
    d2 = dr[1]^2 + dr[2]^2 + eps2
    d  = sqrt(d2)
    # s は正方形領域の辺の長さ（half-width h の2倍）
    s  = 2n.h

    # 子が無い＝葉、または十分遠い＝近似OK
    if (s/d) < θ || all(x -> x === nothing, n.child)
        return pair_force(n.mass, r, n.com)
    end

    # 近いので子へ降りる
    acc = Vec2(0.0, 0.0)
    @inbounds for k in 1:4
        ch = n.child[k]
        if ch !== nothing
            acc += force_from_node(i, r, ch)
        end
    end
    return acc
end

function calc_accel_bh(pos::Vector{Vec2})
    root = build_tree(pos)
    a = Vector{Vec2}(undef, length(pos))
    @inbounds for i in 1:length(pos)
        a[i] = force_from_node(i, pos[i], root)
    end
    return a
end

# 比較用：素朴 O(N^2)
function calc_accel_naive(pos::Vector{Vec2})
    a = [Vec2(0.0, 0.0) for _ in 1:length(pos)]
    @inbounds for i in 1:length(pos)
        ai = Vec2(0.0, 0.0)
        ri = pos[i]
        for j in 1:length(pos)
            if i != j
                ai += pair_force(M, ri, pos[j])
            end
        end
        a[i] = ai
    end
    return a
end

# ----------------------------
# RK4 更新（確保を減らすためワーク領域を持つ）
# ----------------------------
mutable struct RK4Workspace
    # k: dr/dt = v なので、k1..k4 は速度
    k1::Vector{Vec2}; k2::Vector{Vec2}; k3::Vector{Vec2}; k4::Vector{Vec2}
    # l: dv/dt = a なので、l1..l4 は加速度
    l1::Vector{Vec2}; l2::Vector{Vec2}; l3::Vector{Vec2}; l4::Vector{Vec2}
    # 中間状態
    rtmp::Vector{Vec2}
    vtmp::Vector{Vec2}
end

function RK4Workspace(n::Int)
    z = [Vec2(0.0,0.0) for _ in 1:n]
    RK4Workspace(copy(z),copy(z),copy(z),copy(z),
                 copy(z),copy(z),copy(z),copy(z),
                 copy(z),copy(z))
end

function update_rk4!(r::Vector{Vec2}, v::Vector{Vec2}, ws::RK4Workspace; method::Symbol=:bh)
    accel = (method == :bh) ? calc_accel_bh : calc_accel_naive
    n = length(r)

    # k1 = v, l1 = a(r)
    ws.k1 .= v
    ws.l1 .= accel(r)

    # rtmp = r + dt/2*k1, vtmp = v + dt/2*l1
    @inbounds for i in 1:n
        ws.rtmp[i] = r[i] + (dt/2)*ws.k1[i]
        ws.vtmp[i] = v[i] + (dt/2)*ws.l1[i]
    end
    ws.k2 .= ws.vtmp
    ws.l2 .= accel(ws.rtmp)

    @inbounds for i in 1:n
        ws.rtmp[i] = r[i] + (dt/2)*ws.k2[i]
        ws.vtmp[i] = v[i] + (dt/2)*ws.l2[i]
    end
    ws.k3 .= ws.vtmp
    ws.l3 .= accel(ws.rtmp)

    @inbounds for i in 1:n
        ws.rtmp[i] = r[i] + dt*ws.k3[i]
        ws.vtmp[i] = v[i] + dt*ws.l3[i]
    end
    ws.k4 .= ws.vtmp
    ws.l4 .= accel(ws.rtmp)

    # r, v を更新
    @inbounds for i in 1:n
        k = (ws.k1[i] + 2ws.k2[i] + 2ws.k3[i] + ws.k4[i]) / 6
        l = (ws.l1[i] + 2ws.l2[i] + 2ws.l3[i] + ws.l4[i]) / 6
        r[i] = r[i] + dt*k
        v[i] = v[i] + dt*l
    end
    return nothing
end

# ----------------------------
# 描画
# ----------------------------
function save_snapshot(r::Vector{Vec2}, step::Int; outdir::String="img")
    isdir(outdir) || mkpath(outdir)
    x = Vector{Float64}(undef, length(r))
    y = Vector{Float64}(undef, length(r))
    @inbounds for i in 1:length(r)
        x[i] = r[i][1]
        y[i] = r[i][2]
    end
    scatter(x, y, s=1, alpha=0.8)
    xlim(plot_min, plot_max)
    ylim(plot_min, plot_max)
    s = @sprintf("%06d", step)
    savefig(joinpath(outdir, "step"*s*".png"))
    close("all")
end

# ----------------------------
# main
# ----------------------------
function main(; method::Symbol=:bh, do_plot::Bool=false, seed::Int=0)
    r, v = initialize(seed=seed)
    ws = RK4Workspace(length(r))

    for step in 1:max_steps
        println(step)
        update_rk4!(r, v, ws; method=method)
        do_plot && save_snapshot(r, step)
    end
end

# 実行（:bh or :naive）
main(method=:bh, do_plot=false, seed=1)
