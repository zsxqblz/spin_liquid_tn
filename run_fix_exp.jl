include("dependencies.jl")
include("functions.jl")
include("exp.jl")
ITensors.disable_warn_order()

const ncycle = parse(Int64,ARGS[1])
const nex = parse(Int64,ARGS[2])
const nsim = parse(Int64,ARGS[3])
const file_name = ARGS[4]

ncycle_arr = collect(1:ncycle)
exit_cycle_arr = runFixedSizeExp(simFourKagome,ncycle,nex,nsim,true)

save1DData(ncycle_arr,exit_cycle_arr,string(file_name,".csv"))
