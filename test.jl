include("dependencies.jl")
include("functions.jl")
ITensors.disable_warn_order()

let 
    inds,s = genPlusState(4)
    s = appSqProjN(inds,2,s)
    sd = prime(dag(s))
    U, S, V = svd(s*sd, inds...)
    @show S
end

let 
    inds,s = genProdState(2)
    H1 = op("H",inds[1])
    o2 = op("X",inds[2])
    @show o1 * (o1 * o2 * s)
end

let 
    inds,s = genProdState(4)
    # proj = genSqProjN(inds,2)
    # so = s * proj
    prime!(s)
    @show plev(inds[1])
end

# single square plaquette
let 
    inds,s = genPlusState(4)
    H1 = op("H",inds[1])
    H2 = op("H",inds[2])
    H3 = op("H",inds[3])
    H4 = op("H",inds[4])

    ncycle = 7
    suc_rate = Vector{Float64}(undef,ncycle)
    for i = 1:ncycle
        sg = appSqProjN(inds,2,s)
        suc_rate[i] = real((sg*dag(sg))[1])
        s = appSqCProjN(inds,2,s)
        s = appSqHadamard(inds,s)
        # s = s * H1
        # s = s * H2
        # s = s * H3
        # s = s * H4
        s = normalizeState(s)
    end
    s = appSqCProjN(inds,2,s)
    tot_suc_rate = findTotSucRate(suc_rate)
    @show tot_suc_rate
end

# two square vertex sharing
let 
    run(`clear`)
    inds,s = genPlusState(6)

    ncycle = 100
    suc_rate = Vector{Float64}(undef,ncycle)
    for i = 1:ncycle
        sg = appSqProjN(inds[1:4],2,s)
        sg = appSqProjN(inds[3:6],2,sg)
        suc_rate[i] = real((sg*dag(sg))[1])
        s = appSqCProjN(inds[1:4],2,s)
        s = appSqCProjN(inds[3:6],2,s)
        s = appHadamard(inds,s)
        s = normalizeState(s)
    end
    s = appSqCProjN(inds[1:4],2,s)
    s = appSqCProjN(inds[3:6],2,s)
    s = normalizeState(s)
    tot_suc_rate = findTotSucRate(suc_rate)
    @show tot_suc_rate[end-5:end]
    # @show s
    U, S, V = svd(s, inds[1],inds[2],inds[3],inds[4])
    @show diag(S)[1:4]
end

# two square vertex sharing one bad plaquette
let 
    run(`clear`)
    inds,s = genPlusState(6)

    ncycle = 100
    suc_rate = Vector{Float64}(undef,ncycle)
    for i = 1:ncycle
        sg = appSqProjN(inds[1:4],2,s)
        sg = appSqProjN(inds[3:6],2,sg)
        suc_rate[i] = real((sg*dag(sg))[1])
        s = appSqCProjN(inds[1:4],2,s)
        s = appSqProjN(inds[3:6],2,s)
        s = appHadamard(inds[1:4],s)
        s = normalizeState(s)
    end
    s = appSqCProjN(inds[1:4],2,s)
    s = appSqProjN(inds[3:6],2,s)
    s = normalizeState(s)
    tot_suc_rate = findTotSucRate(suc_rate)
    @show tot_suc_rate[end-5:end]
    # @show s
    U, S, V = svd(s, inds[1],inds[2],inds[3],inds[4])
    @show diag(S)[1:4]
end

# three square vertex sharing
let 
    run(`clear`)
    inds,s = genPlusState(8)

    ncycle = 100
    suc_rate = Vector{Float64}(undef,ncycle)
    for i = 1:ncycle
        sg = appSqProjN(inds[1:4],2,s)
        sg = appSqProjN(inds[3:6],2,sg)
        sg = appSqProjN(inds[5:8],2,sg)
        suc_rate[i] = real((sg*dag(sg))[1])
        s = appSqCProjN(inds[1:4],2,s)
        s = appSqCProjN(inds[3:6],2,s)
        s = appSqCProjN(inds[5:8],2,s)
        s = appHadamard(inds,s)
        s = normalizeState(s)
    end
    s = appSqCProjN(inds[1:4],2,s)
    s = appSqCProjN(inds[3:6],2,s)
    s = appSqCProjN(inds[5:8],2,s)
    s = normalizeState(s)
    tot_suc_rate = findTotSucRate(suc_rate)
    @show tot_suc_rate[end-5:end]
    U, S, V = svd(s, inds[1],inds[2],inds[3],inds[4])
    # U, S, V = svd(s, inds[1],inds[2])
    @show diag(S)[1:6]
    # @show s
end

# four square vertex sharing
let 
    run(`clear`)
    inds,s = genPlusState(10)

    ncycle = 101
    suc_rate = Vector{Float64}(undef,ncycle)
    for i = 1:ncycle
        sg = appSqProjN(inds[1:4],2,s)
        sg = appSqProjN(inds[3:6],2,sg)
        sg = appSqProjN(inds[5:8],2,sg)
        sg = appSqProjN(inds[7:10],2,sg)
        suc_rate[i] = real((sg*dag(sg))[1])
        s = appSqCProjN(inds[1:4],2,s)
        s = appSqCProjN(inds[3:6],2,s)
        s = appSqCProjN(inds[5:8],2,s)
        s = appSqCProjN(inds[7:10],2,s)
        s = appHadamard(inds,s)
        s = normalizeState(s)
    end
    s = appSqCProjN(inds[1:4],2,s)
    s = appSqCProjN(inds[3:6],2,s)
    s = appSqCProjN(inds[5:8],2,s)
    s = appSqCProjN(inds[7:10],2,s)
    s = normalizeState(s)
    tot_suc_rate = findTotSucRate(suc_rate)
    @show tot_suc_rate[end-5:end]
    U, S, V = svd(s, inds[1],inds[2],inds[3],inds[4])
    @show diag(S)[1:6]
end

# two square edge sharing
let
    run(`clear`)
    inds,s = genPlusState(7)
    ncycle = 100
    suc_rate = Vector{Float64}(undef,ncycle)
    for i = 1:ncycle
        sg = appSqProjN(inds[1:4],2,s)
        sg = appSqProjN(inds[4:7],2,sg)
        suc_rate[i] = real((sg*dag(sg))[1])
        s = appSqCProjN(inds[1:4],2,s)
        s = appSqCProjN(inds[4:7],2,s)
        s = appHadamard(inds,s)
        s = normalizeState(s)
    end
    s = appSqCProjN(inds[1:4],2,s)
    s = appSqCProjN(inds[4:7],2,s)
    s = normalizeState(s)
    tot_suc_rate = findTotSucRate(suc_rate)
    @show tot_suc_rate[end-5:end]
    # @show s
    U, S, V = svd(s, inds[1],inds[2],inds[3],inds[4])
    @show diag(S)[1:6]
end


# single hex plaquette
let 
    run(`clear`)
    inds,s = genPlusState(6)
    ncycle = 101
    suc_rate = Vector{Float64}(undef,ncycle)
    for i = 1:ncycle
        sg = appHexProjN(inds,2,s)
        suc_rate[i] = real((sg*dag(sg))[1])
        s = appHexCProjN(inds,2,s)
        s = appHexHadamard(inds,s)
        s = normalizeState(s)
    end
    s = appHexCProjN(inds,2,s)
    tot_suc_rate = findTotSucRate(suc_rate)
    @show tot_suc_rate[end-5:end]
    # @show s
end


# two hex plaquette shared vertex
let 
    run(`clear`)
    inds,s = genPlusState(11)
    ncycle = 101
    suc_rate = Vector{Float64}(undef,ncycle)
    for i = 1:ncycle
        sg = appHexProjN(inds[1:6],2,s)
        sg = appHexProjN(inds[6:11],2,sg)
        suc_rate[i] = real((sg*dag(sg))[1])
        s = appHexCProjN(inds[1:6],2,s)
        s = appHexCProjN(inds[6:11],2,s)
        s = appHadamard(inds,s)
        s = normalizeState(s)
    end
    s = appHexCProjN(inds[1:6],2,s)
    s = appHexCProjN(inds[6:11],2,s)
    tot_suc_rate = findTotSucRate(suc_rate)
    @show tot_suc_rate[end-5:end]
    # @show appHadamard(inds,s)
end

# three hex plaquette shared vertex
let 
    run(`clear`)
    inds,s = genPlusState(15)
    ncycle = 100
    suc_rate = Vector{Float64}(undef,ncycle)
    inds_hex1 = inds[1:6]
    inds_hex2 = inds[6:11]
    inds_hex3 = [inds[5],inds[11],inds[12],inds[13],inds[14],inds[15]]
    for i = 1:ncycle
        sg = appHexProjN(inds_hex1,1,s)
        sg = appHexProjN(inds_hex2,1,sg)
        sg = appHexProjN(inds_hex3,1,sg)
        suc_rate[i] = real((sg*dag(sg))[1])
        s = appHexCProjN(inds_hex1,1,s)
        s = appHexCProjN(inds_hex2,1,s)
        s = appHexCProjN(inds_hex3,1,s)
        s = appHadamard(inds,s)
        # s = appHadamard([inds[5],inds[6],inds[11]],s)
        s = normalizeState(s)
    end
    s = appHexCProjN(inds_hex1,1,s)
    s = appHexCProjN(inds_hex2,1,s)
    s = appHexCProjN(inds_hex3,1,s)
    tot_suc_rate = findTotSucRate(suc_rate)
    @show tot_suc_rate[end-5:end]
    # @show appHadamard(inds,s)
end