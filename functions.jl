function prepInds(nsites)
    inds_l = Vector{Index{Int64}}(undef,nsites)
    for i = 1:nsites
        inds_l[i] = Index(2,"Site,Qubit")
    end
    return inds_l
end

function genProdState(nsites::Int64)
    inds_l = prepInds(nsites)
    state = ITensor(ComplexF64,inds_l)
    state[ones(Int64,nsites)...] = 1
    return inds_l,state
end

function genProdState(inds_l::Vector{Index{Int64}})
    nsites = length(inds_l)
    state = ITensor(ComplexF64,inds_l)
    state[ones(Int64,nsites)...] = 1
    return state
end

function genPlusState(inds_l::Vector{Index{Int64}})
    nsites = length(inds_l)
    state = ITensor(ComplexF64,inds_l)
    for iv in eachindval(inds_l...)
        state[iv...] = 1/(sqrt(2)^nsites)
    end
    return state
end

function genPlusState(nsites::Int64)
    inds_l = prepInds(nsites)
    state = ITensor(ComplexF64,inds_l)
    for iv in eachindval(inds_l...)
        state[iv...] = 1/(sqrt(2)^nsites)
    end
    return inds_l,state
end

function normalizeState(s::ITensor)
    return s / sqrt((s*dag(s))[1])
end

function genTrgProjN(inds_l::Vector{Index{Int64}},nex::Int64)
    inds_o = prime(dag(inds_l))
    proj = ITensor(ComplexF64,vcat(inds_l,inds_o))
    for iv in Iterators.product(1:2,1:2,1:2)
        if sum(iv) - 3 == nex
            proj[iv...,iv...] = 1
        end
    end
    return proj
end

function appTrgProjN(inds_l::Vector{Index{Int64}},nex::Int64,s::ITensor)
    proj = genTrgProjN(inds_l,nex)
    so = noprime(proj * s)
    return so
end

function genTrgCProjN(inds_l::Vector{Index{Int64}},nex::Int64)
    inds_o = prime(dag(inds_l))
    proj = ITensor(ComplexF64,vcat(inds_l,inds_o))
    for iv in Iterators.product(1:2,1:2,1:2)
        if sum(iv) - 3 != nex
            proj[iv...,iv...] = 1
        end
    end
    return proj
end

function appTrgCProjN(inds_l::Vector{Index{Int64}},nex::Int64,s::ITensor)
    proj = genTrgCProjN(inds_l,nex)
    so = noprime(proj * s)
    return so
end

function appTrgHadamard(inds::Vector{Index{Int64}},s::ITensor)
    so = s * op("H",inds[1])
    so = so * op("H",inds[2])
    so = so * op("H",inds[3])
    return so
end

function genSqProjN(inds_l::Vector{Index{Int64}},nex::Int64)
    inds_o = prime(dag(inds_l))
    proj = ITensor(ComplexF64,vcat(inds_l,inds_o))
    for iv in Iterators.product(1:2,1:2,1:2,1:2)
        if sum(iv) - 4 == nex
            proj[iv...,iv...] = 1
        end
    end
    return proj
end

function appSqProjN(inds_l::Vector{Index{Int64}},nex::Int64,s::ITensor)
    proj = genSqProjN(inds_l,nex)
    so = noprime(proj * s)
    return so
end

function genSqCProjN(inds_l::Vector{Index{Int64}},nex::Int64)
    inds_o = prime(dag(inds_l))
    proj = ITensor(ComplexF64,vcat(inds_l,inds_o))
    for iv in Iterators.product(1:2,1:2,1:2,1:2)
        if sum(iv) - 4 != nex
            proj[iv...,iv...] = 1
        end
    end
    return proj
end

function appSqCProjN(inds_l::Vector{Index{Int64}},nex::Int64,s::ITensor)
    proj = genSqCProjN(inds_l,nex)
    so = noprime(proj * s)
    return so
end

function appSqHadamard(inds::Vector{Index{Int64}},s::ITensor)
    so = s * op("H",inds[1])
    so = so * op("H",inds[2])
    so = so * op("H",inds[3])
    so = so * op("H",inds[4])
    return so
end

function appHadamard(inds::Vector{Index{Int64}},s::ITensor)
    so = copy(s)
    for ind in inds
        so = so * op("H",ind)
    end
    return so
end

function genHexProjN(inds_l::Vector{Index{Int64}},nex::Int64)
    inds_o = prime(dag(inds_l))
    proj = ITensor(ComplexF64,vcat(inds_l,inds_o))
    for iv in Iterators.product(1:2,1:2,1:2,1:2,1:2,1:2)
        if sum(iv) - 6 == nex
            proj[iv...,iv...] = 1
        end
    end
    return proj
end

function appHexProjN(inds_l::Vector{Index{Int64}},nex::Int64,s::ITensor)
    proj = genHexProjN(inds_l,nex)
    so = noprime(proj * s)
    return so
end

function genHexCProjN(inds_l::Vector{Index{Int64}},nex::Int64)
    inds_o = prime(dag(inds_l))
    proj = ITensor(ComplexF64,vcat(inds_l,inds_o))
    for iv in Iterators.product(1:2,1:2,1:2,1:2,1:2,1:2)
        if sum(iv) - 6 != nex
            proj[iv...,iv...] = 1
        end
    end
    return proj
end

function appHexCProjN(inds_l::Vector{Index{Int64}},nex::Int64,s::ITensor)
    proj = genHexCProjN(inds_l,nex)
    so = noprime(proj * s)
    return so
end

function appHexHadamard(inds::Vector{Index{Int64}},s::ITensor)
    so = s * op("H",inds[1])
    so = so * op("H",inds[2])
    so = so * op("H",inds[3])
    so = so * op("H",inds[4])
    so = so * op("H",inds[5])
    so = so * op("H",inds[6])
    return so
end

function findTotSucRate(suc_rate::Vector{Float64})
    ncycle = length(suc_rate)
    tot_suc_rate = Vector{Float64}(undef,ncycle)
    tot_suc_rate[1] = suc_rate[1]
    for i = 2:ncycle
        tot_suc_rate[i] = tot_suc_rate[i-1] + (1-tot_suc_rate[i-1])*suc_rate[i]
    end
    return tot_suc_rate
end

function findStateSum(s::ITensor)
    sum = 0
    for iv in eachindval(s.inds)
        sum = sum + real(s[iv...])
    end
    return sum
end

function findStateAbsSum(s::ITensor)
    sum = 0
    for iv in eachindval(s.inds)
        sum = sum + abs(s[iv...])
    end
    return sum
end

function get_sites(x::Int, y::Int, dx::Int, dy::Int)
    # calculate the site numbers of the edges of the plaquette
    site1 = (x-1)*dy + y
    site2 = (x-1)*(dy+1) + y
    site3 = x*(dy+1) + y
    site4 = x*dy + y
    
    # return the site numbers
    return site1, site2, site3, site4
end
