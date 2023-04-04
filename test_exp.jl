include("dependencies.jl")
include("functions.jl")
include("exp.jl")
ITensors.disable_warn_order()

let
    run(`clear`)
    dx = 2
    dy = 3
    ncycle = 100
    nex = 2
    nsim = 1000
    ncycle_arr = collect(1:ncycle)
    exit_cycle_arr = runFixedSizeExp(simThreeKagome,ncycle,nex,nsim,true)
    
    save_idx = 16
    save1DData(ncycle_arr,exit_cycle_arr,string("data/230403/230403_",save_idx,".csv"))
end

let 
    run(`clear`)
    dx = 1
    dy = 2
    inds,s = genPlusState((dx+1)*(dy+1))
    inds_arr = [inds[1],inds[2],inds[dy+2],inds[dy+3]]
    for x = 0:(dx-1)
        for y = 0:(dy-1)
            if (x == 0) & (y == 0)
                continue
            end
            @show x,y
            inds_arr = hcat(inds_arr,[inds[x*(dy+1)+y+1],inds[x*(dy+1)+y+2],inds[(x+1)*(dy+1)+y+1],inds[(x+1)*(dy+1)+y+2]])
        end
    end
    # inds_arr = transpose(inds_arr)
    # @show inds
    @show inds_arr
    # @show inds[5],inds[6],inds[8],inds[9]
end