function simOneSqVertex(ncycle::Int64,nex::Int64)
    inds,s = genPlusState(4)
    for i = 1:ncycle
        prob1 = rand()
        s1 = appSqCProjN(inds[1:4],nex,s)
        if prob1 < real((s1*dag(s1))[1]) 
            s = appHadamard(inds[1:4],s1)
        else
            return true, i
        end
        s = normalizeState(s)
    end
    return false, -1
end

function simTwoSqVertex(ncycle::Int64,nex::Int64)
    inds,s = genPlusState(6)
    for i = 1:ncycle
        prob1 = rand()
        prob2 = rand()
        s1 = appSqCProjN(inds[1:4],nex,s)
        if prob1 < real((s1*dag(s1))[1]) 
            s2 = appSqCProjN(inds[3:6],nex,normalizeState(s1))
            if prob2 < real((s2*dag(s2))[1])
                s = normalizeState(s2)
                s = appHadamard(inds,s)
            else
                s = appSqProjN(inds[3:6],nex,s1)
                s = normalizeState(s)
                s = appHadamard(inds[1:4],s)
            end
        else
            s1 = appSqProjN(inds[1:4],nex,s)
            s2 = appSqCProjN(inds[3:6],nex,normalizeState(s1))
            if prob2 < real((s2*dag(s2))[1])
                s = normalizeState(s2)
                s = appHadamard(inds[3:6],s)
            else
                return true, i
            end
        end
    end
    return false, -1
end

function simFourSqVertex(ncycle::Int64,nex::Int64)
    inds,s = genPlusState(9)
    inds1 = [inds[1],inds[2],inds[4],inds[5]]
    inds2 = [inds[2],inds[3],inds[5],inds[6]]
    inds3 = [inds[4],inds[5],inds[7],inds[8]]
    inds4 = [inds[5],inds[6],inds[8],inds[9]]
    inds_arr = [inds1,inds2,inds3,inds4]
    si = s
    for i = 1:ncycle
        error_flag = [false,false,false,false]
        for j = 1:4
            sp = appSqCProjN(inds_arr[j],nex,si)
            prob = rand()
            if prob < real((sp*dag(sp))[1])
                si = normalizeState(sp)
                error_flag[j] = true
            else
                sp = appSqProjN(inds_arr[j],nex,si)
                si = normalizeState(sp)
            end
        end

        if !(error_flag[1] | error_flag[2] | error_flag[3] | error_flag[4])
            return true, i
        end

        hadamard_flags = zeros(Bool,9)
        if error_flag[1]
            hadamard_flags[1] = true
            hadamard_flags[2] = true
            hadamard_flags[4] = true
            hadamard_flags[5] = true
        end
        if error_flag[2]
            hadamard_flags[2] = true
            hadamard_flags[3] = true
            hadamard_flags[5] = true
            hadamard_flags[6] = true
        end
        if error_flag[3]
            hadamard_flags[4] = true
            hadamard_flags[5] = true
            hadamard_flags[7] = true
            hadamard_flags[8] = true
        end
        if error_flag[4]
            hadamard_flags[5] = true
            hadamard_flags[6] = true
            hadamard_flags[8] = true
            hadamard_flags[9] = true
        end

        for (j,hadamard_flag) in enumerate(hadamard_flags)
            if hadamard_flag
                si = appHadamard(inds[j:j],si)
            end
        end
    end
    return false, -1
end

function simSqVertex(ncycle::Int64,nex::Int64,dx::Int64,dy::Int64)
    inds,s = genPlusState((dx+1)*(dy+1))
    inds_arr = [inds[1],inds[2],inds[dy+2],inds[dy+3]]
    for x = 0:(dx-1)
        for y = 0:(dy-1)
            if (x == 0) & (y == 0)
                continue
            end
            inds_arr = hcat(inds_arr,[inds[x*(dy+1)+y+1],inds[x*(dy+1)+y+2],inds[(x+1)*(dy+1)+y+1],inds[(x+1)*(dy+1)+y+2]])
        end
    end

    si = s
    for i = 1:ncycle
        error_flag = zeros(Bool,dx*dy)
        # @show inds_arr
        for j = 1:dx*dy
            sp = appSqCProjN(inds_arr[:,j],nex,si)
            prob = rand()
            if prob < real((sp*dag(sp))[1])
                si = normalizeState(sp)
                error_flag[j] = true
            else
                sp = appSqProjN(inds_arr[:,j],nex,si)
                si = normalizeState(sp)
            end
        end

        wrong_state = false
        for j = 1:dx*dy
            wrong_state = wrong_state | error_flag[j]
        end
        if !wrong_state
            return true,i
        end 

        hadamard_flags = zeros(Bool,(dx+1)*(dy+1))
        for x = 0:(dx-1)
            for y = 0:(dy-1)
                j = x*dy + y + 1
                if error_flag[j]
                    hadamard_flags[x*(dy+1)+y+1] = true
                    hadamard_flags[x*(dy+1)+y+2] = true
                    hadamard_flags[(x+1)*(dy+1)+y+1] = true
                    hadamard_flags[(x+1)*(dy+1)+y+2] = true
                end
            end
        end

        for (j,hadamard_flag) in enumerate(hadamard_flags)
            if hadamard_flag
                si = appHadamard(inds[j:j],si)
            end
        end
    end
    return false, -1
end

function simOneSqEdge(ncycle::Int64,nex::Int64)
    inds,s = genPlusState(4)
    inds1 = [inds[1],inds[2],inds[3],inds[4]]
    inds_arr = [inds1]
    si = s
    for i = 1:ncycle
        error_flag = [false]
        for j = 1:1
            sp = appSqCProjN(inds_arr[j],nex,si)
            prob = rand()
            if prob < real((sp*dag(sp))[1])
                si = normalizeState(sp)
                error_flag[j] = true
            else
                sp = appSqProjN(inds_arr[j],nex,si)
                si = normalizeState(sp)
            end
        end

        if !(error_flag[1])
            return true, i
        end

        hadamard_flags = zeros(Bool,4)
        if error_flag[1]
            hadamard_flags[1] = true
            hadamard_flags[2] = true
            hadamard_flags[3] = true
            hadamard_flags[4] = true
        end
        for (j,hadamard_flag) in enumerate(hadamard_flags)
            if hadamard_flag
                si = appHadamard(inds[j:j],si)
            end
        end
    end
    return false, -1
end

function simTwoSqEdge(ncycle::Int64,nex::Int64)
    inds,s = genPlusState(4)
    inds1 = [inds[1],inds[2],inds[3],inds[4]]
    inds2 = [inds[4],inds[5],inds[6],inds[7]]
    inds3 = [inds[2],inds[8],inds[9],inds[10]]
    inds2 = [inds[5],inds[10],inds[11],inds[12]]
    inds_arr = [inds1,inds2,inds3,inds4]
    si = s
    for i = 1:ncycle
        error_flag = [false]
        for j = 1:1
            sp = appSqCProjN(inds_arr[j],nex,si)
            prob = rand()
            if prob < real((sp*dag(sp))[1])
                si = normalizeState(sp)
                error_flag[j] = true
            else
                sp = appSqProjN(inds_arr[j],nex,si)
                si = normalizeState(sp)
            end
        end

        if !(error_flag[1]|error_flag[2]|error_flag[3]|error_flag[4])
            return true, i
        end

        hadamard_flags = zeros(Bool,4)
        if error_flag[1]
            hadamard_flags[1] = true
            hadamard_flags[2] = true
            hadamard_flags[3] = true
            hadamard_flags[4] = true
        end
        if error_flag[2]
            hadamard_flags[4] = true
            hadamard_flags[5] = true
            hadamard_flags[6] = true
            hadamard_flags[7] = true
        end
        if error_flag[3]
            hadamard_flags[2] = true
            hadamard_flags[8] = true
            hadamard_flags[9] = true
            hadamard_flags[10] = true
        end
        if error_flag[4]
            hadamard_flags[5] = true
            hadamard_flags[10] = true
            hadamard_flags[11] = true
            hadamard_flags[12] = true
        end

        for (j,hadamard_flag) in enumerate(hadamard_flags)
            if hadamard_flag
                si = appHadamard(inds[j:j],si)
            end
        end
    end
    return false, -1
end

function simOneKagome(ncycle::Int64,nex::Int64)
    inds,s = genPlusState(6)
    for i = 1:ncycle
        prob1 = rand()
        s1 = appHexCProjN(inds[1:6],nex,s)
        if prob1 < real((s1*dag(s1))[1]) 
            s = appHadamard(inds[1:6],s1)
        else
            return true, i, s
        end
        s = normalizeState(s)
    end
    return false, -1, s
end

function simTwoKagome(ncycle::Int64,nex::Int64)
    inds,s = genPlusState(11)
    inds1 = [inds[1],inds[2],inds[3],inds[4],inds[5],inds[6]]
    inds2 = [inds[6],inds[7],inds[8],inds[9],inds[10],inds[11]]
    inds_arr = hcat(inds1,inds2)
    si = s
    for i = 1:ncycle
        error_flag = [false,false]
        for j = 1:2
            sp = appHexCProjN(inds_arr[:,j],nex,si)
            prob = rand()
            if prob < real((sp*dag(sp))[1])
                si = normalizeState(sp)
                error_flag[j] = true
            else
                sp = appHexProjN(inds_arr[:,j],nex,si)
                si = normalizeState(sp)
            end
        end

        if !(error_flag[1] | error_flag[2])
            return true, i, s
        end

        hadamard_flags = zeros(Bool,11)
        if error_flag[1]
            hadamard_flags[1] = true
            hadamard_flags[2] = true
            hadamard_flags[3] = true
            hadamard_flags[4] = true
            hadamard_flags[5] = true
            hadamard_flags[6] = true
        end
        if error_flag[2]
            hadamard_flags[6] = true
            hadamard_flags[7] = true
            hadamard_flags[8] = true
            hadamard_flags[9] = true
            hadamard_flags[10] = true
            hadamard_flags[11] = true
        end

        for (j,hadamard_flag) in enumerate(hadamard_flags)
            if hadamard_flag
                si = appHadamard(inds[j:j],si)
            end
        end
    end
    return false, -1, s
end

function simThreeKagome(ncycle::Int64,nex::Int64)
    inds,s = genPlusState(15)
    inds1 = [inds[1],inds[2],inds[3],inds[4],inds[5],inds[6]]
    inds2 = [inds[6],inds[7],inds[8],inds[9],inds[10],inds[11]]
    inds3 = [inds[5],inds[11],inds[12],inds[13],inds[14],inds[15]]
    inds_arr = hcat(inds1,inds2,inds3)
    si = s
    for i = 1:ncycle
        error_flag = [false,false,false]
        for j = 1:3
            sp = appHexCProjN(inds_arr[:,j],nex,si)
            prob = rand()
            if prob < real((sp*dag(sp))[1])
                si = normalizeState(sp)
                error_flag[j] = true
            else
                sp = appHexProjN(inds_arr[:,j],nex,si)
                si = normalizeState(sp)
            end
        end

        if !(error_flag[1] | error_flag[2] | error_flag[3])
            return true, i, s
        end

        hadamard_flags = zeros(Bool,15)
        if error_flag[1]
            hadamard_flags[1] = true
            hadamard_flags[2] = true
            hadamard_flags[3] = true
            hadamard_flags[4] = true
            hadamard_flags[5] = true
            hadamard_flags[6] = true
        end
        if error_flag[2]
            hadamard_flags[6] = true
            hadamard_flags[7] = true
            hadamard_flags[8] = true
            hadamard_flags[9] = true
            hadamard_flags[10] = true
            hadamard_flags[11] = true
        end
        if error_flag[3]
            hadamard_flags[5] = true
            hadamard_flags[11] = true
            hadamard_flags[12] = true
            hadamard_flags[13] = true
            hadamard_flags[14] = true
            hadamard_flags[15] = true
        end

        for (j,hadamard_flag) in enumerate(hadamard_flags)
            if hadamard_flag
                si = appHadamard(inds[j:j],si)
            end
        end
    end
    return false, -1, s
end

function simFourKagome(ncycle::Int64,nex::Int64)
    inds,s = genPlusState(19)
    inds1 = [inds[1],inds[2],inds[3],inds[4],inds[5],inds[6]]
    inds2 = [inds[6],inds[7],inds[8],inds[9],inds[10],inds[11]]
    inds3 = [inds[5],inds[11],inds[12],inds[13],inds[14],inds[15]]
    inds4 = [inds[10],inds[12],inds[16],inds[17],inds[18],inds[19]]
    inds_arr = hcat(inds1,inds2,inds3,inds4)
    si = s
    for i = 1:ncycle
        error_flag = [false,false,false,false]
        for j = 1:4
            sp = appHexCProjN(inds_arr[:,j],nex,si)
            prob = rand()
            if prob < real((sp*dag(sp))[1])
                si = normalizeState(sp)
                error_flag[j] = true
            else
                sp = appHexProjN(inds_arr[:,j],nex,si)
                si = normalizeState(sp)
            end
        end

        if !(error_flag[1] | error_flag[2] | error_flag[3] | error_flag[4])
            return true, i, s
        end

        hadamard_flags = zeros(Bool,19)
        if error_flag[1]
            hadamard_flags[1] = true
            hadamard_flags[2] = true
            hadamard_flags[3] = true
            hadamard_flags[4] = true
            hadamard_flags[5] = true
            hadamard_flags[6] = true
        end
        if error_flag[2]
            hadamard_flags[6] = true
            hadamard_flags[7] = true
            hadamard_flags[8] = true
            hadamard_flags[9] = true
            hadamard_flags[10] = true
            hadamard_flags[11] = true
        end
        if error_flag[3]
            hadamard_flags[5] = true
            hadamard_flags[11] = true
            hadamard_flags[12] = true
            hadamard_flags[13] = true
            hadamard_flags[14] = true
            hadamard_flags[15] = true
        end
        if error_flag[4]
            hadamard_flags[10] = true
            hadamard_flags[12] = true
            hadamard_flags[16] = true
            hadamard_flags[17] = true
            hadamard_flags[18] = true
            hadamard_flags[19] = true
        end

        for (j,hadamard_flag) in enumerate(hadamard_flags)
            if hadamard_flag
                si = appHadamard(inds[j:j],si)
            end
        end
    end
    return false, -1, s
end

function expOneSqVertex(ncycle::Int64,nex::Int64,nsim::Int64,show_prog::Bool=false)
    exit_cycle_arr = zeros(ncycle)
    if show_prog
        @showprogress for i = 1:nsim
            suc, exit_cycle_idx = simOneSqVertex(ncycle,nex)
            if suc
                exit_cycle_arr[exit_cycle_idx] = exit_cycle_arr[exit_cycle_idx] + 1
            end
        end
    else
        for i = 1:nsim
            suc, exit_cycle_idx = simOneSqVertex(ncycle,nex)
            if suc
                exit_cycle_arr[exit_cycle_idx] = exit_cycle_arr[exit_cycle_idx] + 1
            end
        end
    end
    return exit_cycle_arr / nsim
end

function expTwoSqVertex(ncycle::Int64,nex::Int64,nsim::Int64,show_prog::Bool=false)
    exit_cycle_arr = zeros(ncycle)
    if show_prog
        @showprogress for i = 1:nsim
            suc, exit_cycle_idx = simTwoSqVertex(ncycle,nex)
            if suc
                exit_cycle_arr[exit_cycle_idx] = exit_cycle_arr[exit_cycle_idx] + 1
            end
        end
    else
        for i = 1:nsim
            suc, exit_cycle_idx = simTwoSqVertex(ncycle,nex)
            if suc
                exit_cycle_arr[exit_cycle_idx] = exit_cycle_arr[exit_cycle_idx] + 1
            end
        end
    end
    return exit_cycle_arr / nsim
end

function expFourSqVertex(ncycle::Int64,nex::Int64,nsim::Int64,show_prog::Bool=false)
    exit_cycle_arr = zeros(ncycle)
    if show_prog
        @showprogress for i = 1:nsim
            suc, exit_cycle_idx = simFourSqVertex(ncycle,nex)
            if suc
                exit_cycle_arr[exit_cycle_idx] = exit_cycle_arr[exit_cycle_idx] + 1
            end
        end
    else
        for i = 1:nsim
            suc, exit_cycle_idx = simFourSqVertex(ncycle,nex)
            if suc
                exit_cycle_arr[exit_cycle_idx] = exit_cycle_arr[exit_cycle_idx] + 1
            end
        end
    end
    return exit_cycle_arr / nsim
end

function expSqVertex(ncycle::Int64,nex::Int64,dx::Int64,dy::Int64,nsim::Int64,show_prog::Bool=false)
    exit_cycle_arr = zeros(ncycle)
    if show_prog
        @showprogress for i = 1:nsim
            suc, exit_cycle_idx = simSqVertex(ncycle,nex,dx,dy)
            if suc
                exit_cycle_arr[exit_cycle_idx] = exit_cycle_arr[exit_cycle_idx] + 1
            end
        end
    else
        for i = 1:nsim
            suc, exit_cycle_idx = simSqVertex(ncycle,nex,dx,dy)
            if suc
                exit_cycle_arr[exit_cycle_idx] = exit_cycle_arr[exit_cycle_idx] + 1
            end
        end
    end
    return exit_cycle_arr / nsim
end

function runFixedSizeExp(func,ncycle::Int64,nex::Int64,nsim::Int64,show_prog::Bool=false)
    exit_cycle_arr = zeros(ncycle)
    if show_prog
        @showprogress for i = 1:nsim
            suc, exit_cycle_idx = func(ncycle,nex)
            if suc
                exit_cycle_arr[exit_cycle_idx] = exit_cycle_arr[exit_cycle_idx] + 1
            end
        end
    else
        for i = 1:nsim
            suc, exit_cycle_idx = func(ncycle,nex)
            if suc
                exit_cycle_arr[exit_cycle_idx] = exit_cycle_arr[exit_cycle_idx] + 1
            end
        end
    end
    return exit_cycle_arr / nsim
end

function save1DData(x_l,y_l,file_name)
    df = DataFrame()
    df.x_l = x_l
    df.y_l = y_l
    CSV.write(file_name, df)
end