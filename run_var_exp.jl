include("dependencies.jl")
include("functions.jl")
include("exp.jl")
ITensors.disable_warn_order()

const dx = parse(Int64,ARGS[1])
const dy = parse(Int64,ARGS[2])
const ncycle = parse(Int64,ARGS[3])
const nex = parse(Int64,ARGS[4])
const nsim = parse(Int64,ARGS[5])
const file_name = ARGS[6]

ncycle_arr = collect(1:ncycle)
exit_cycle_arr = expSqVertex(ncycle,nex,dx,dy,nsim,true)

save1DData(ncycle_arr,exit_cycle_arr,string(file_name,".csv"))
