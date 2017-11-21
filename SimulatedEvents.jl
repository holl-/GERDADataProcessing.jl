using ROOTFramework, Cxx
using Clustering

using MJDSigGen: fieldgen, signal_calc_init, get_signal!, drift_path, read_fields, outside_detector, nearest_field_grid_index, cyl2cart, cart2cyl

type Event
    hitpositions
    hitenergies
    hitpulses
    hitepaths
    hithpaths
    hitclusters
    clusterpositions
    clustersigmas
    clusterenergies
    clusterpulses
    clusterepaths
    clusterhpaths
    totenergy
    sumpulse
end

#function to read in converted MaGe events and generate pulses from them
function create_events(setup::MJDSigGen.Struct_MJD_Siggen_Setup, root_file_name, root_tree_name; min_energy = 0, max_energy = 9999, save_hits = false, save_clusters = false, cluster_size = -1)

    #input is a converted MaGe root file that has a hittree
    #min_energy = 600 is useful to save memory and storage space
    #save_hits = true saves all hits: positions, energies, pulses, drift paths
    #save_clusters = true saves all clusters: positions, sizes, energies, pulses, drift paths
    #cluster_size <= 0 or Nhits <= 3 the sum of hit pulses will be saved, otherwise the sum of cluster pulses will be saved
    #clustering only works for more than 3 hits :S

    if save_clusters && cluster_size <= 0
        error("Error: cannot save clusters with size <= 0")
    end

    #root file input
    bindings = TTreeBindings()

    Nhits = bindings[:Nhits] = Ref(zero(Int32))
    TotEnergy = bindings[:TotEnergy] = Ref(zero(Float32))
    HitEnergy = bindings[:HitEnergy] = ROOTFramework.CxxObjWithPtrRef(icxx"std::vector<float>();")
    HitX = bindings[:HitX] = ROOTFramework.CxxObjWithPtrRef(icxx"std::vector<float>();")
    HitY = bindings[:HitY] = ROOTFramework.CxxObjWithPtrRef(icxx"std::vector<float>();")
    HitZ = bindings[:HitZ] = ROOTFramework.CxxObjWithPtrRef(icxx"std::vector<float>();")

    #initialize empty output
    allevents = Array{Event}(0)
    
    #loop over root file
    open(TChainInput, bindings, root_tree_name, root_file_name) do input

        i = 0
        for _ in input

            i += 1
            if TotEnergy.x < min_energy
                continue
            end
            if TotEnergy.x > max_energy
                continue
            end

            #initializations
            eventpulse = zero(Float32)
            goodevent = true #workaround for events with hits at the contact or groove

            setup.charge_cloud_size = 0.0
            setup.use_diffusion = 1

            if cluster_size > 0 || save_hits || save_clusters 
                positions = zeros(3, Nhits.x) #needs to stay Float64 because of the clustering algorithm
                energies = zeros(Float32, Nhits.x)
                if save_hits
                    pulses = Array{Array{Float32}}(0)
                    epaths = Array{Array{Float32,2}}(0)
                    hpaths = Array{Array{Float32,2}}(0)
                end
                if save_clusters
                    clpositions = Array{Array{Float32}}(0)
                    clsigmas = Array{Float32}(0)
                    clenergies = Array{Float32}(0)
                    clpulses = Array{Array{Float32}}(0)
                    clepaths = Array{Array{Float32,2}}(0)
                    clhpaths = Array{Array{Float32,2}}(0)
                end
            end

            #loop over hits
            for j in 0:(Nhits.x-1)
                if cluster_size > 0 || save_hits || save_clusters 
                    positions[1,j+1], positions[2,j+1], positions[3,j+1] = HitX.x[j], HitY.x[j], HitZ.x[j]
                    energies[j+1] = HitEnergy.x[j]
                end
                if cluster_size <= 0 || Nhits.x <= 3 || save_hits
                    p1 = (HitX.x[j], HitY.x[j], HitZ.x[j])
                    ind = nearest_field_grid_index(setup, p1)
                    setup.energy = HitEnergy.x[j]
                    if ind[1] != :outside
                        if ind[1] == :interpol
                            pulse1 = get_signal!(setup, p1)
                        elseif ind[1] == :extrapol
                            #println("event $i point need extrapolation $p1")
                            phi = cart2cyl(p1[1],p1[2],p1[3])[2]
                            p1 = cyl2cart(ind[2]*setup.xtal_grid, phi, ind[3]*setup.xtal_grid)
                            pulse1 = get_signal!(setup, p1)
                        else
                            println("this should not happen")
                        end
                        if cluster_size <= 0 || Nhits.x <= 3
                            eventpulse += HitEnergy.x[j] * pulse1
                        end
                        if save_hits
                            push!(pulses, HitEnergy.x[j] * pulse1)
                            push!(epaths, drift_path(setup, :e))
                            push!(hpaths, drift_path(setup, :h))
                        end
                    else
                        println("event $i skipped because point outside $p1")
                        goodevent = false
                    end
                end
            end

            #do clustering
            if cluster_size > 0 && Nhits.x > 3
                clusters = dbscan(positions, cluster_size, min_neighbors=1, min_cluster_size=1, leafsize=20)
                for j in 1:length(clusters)
                    inds = vcat(clusters[j].core_indices, clusters[j].boundary_indices)
                    Ecl = sum(energies[inds])
                    p2 = [sum(energies[inds].*positions[1,inds]); sum(energies[inds].*positions[2,inds]); sum(energies[inds].*positions[3,inds])]./Ecl
                    pulse2 = zero(Float32)
                    check = nearest_field_grid_index(setup, tuple(p2...))
                    if check[1] != :outside
                        clsize = 0
                        for k=1:length(inds)
                            clsize += energies[inds[k]] * (norm(positions[:,inds[k]]-p2)^2)
                        end
                        clsize = sqrt(clsize/Ecl)
                        setup.charge_cloud_size = clsize
                        setup.energy = Ecl
                        if check[1] == :interpol
                            pulse2 = get_signal!(setup, tuple(p2...))
                        elseif check[1] == :extrapol
                            #println("event $i cluster need extrapolation $p2")
                            phi = cart2cyl(p2[1],p2[2],p2[3])[2]
                            p2[1],p2[2],p2[3] = cyl2cart(check[2]*setup.xtal_grid, phi, check[3]*setup.xtal_grid)
                            checkagain = nearest_field_grid_index(setup, tuple(p2...))
                            if checkagain[1] == :interpol
                                pulse2 = get_signal!(setup, tuple(p2...))
                            else
                                println("event $i skipped because extrapolated cluster outside $p2")
                                goodevent = false
                            end
                        else
                            println("this should not happen")
                        end
                        if save_clusters && goodevent
                            push!(clenergies, Ecl)
                            push!(clpositions, p2)
                            push!(clsigmas, clsize)
                            push!(clpulses, Ecl * pulse2)
                            push!(clepaths, drift_path(setup, :e))
                            push!(clhpaths, drift_path(setup, :h))
                        end
                        eventpulse += Ecl * pulse2
                    else
                        println("event $i skipped because cluster outside $p2")
                        goodevent = false
                    end
                end
            end

            #save everything
            if goodevent
                event1 = Event(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, TotEnergy.x, eventpulse)
                if save_hits
                    event1.hitpositions = positions
                    event1.hitenergies = energies
                    event1.hitpulses = pulses
                    event1.hitepaths = epaths
                    event1.hithpaths = hpaths
                    if cluster_size > 0 && Nhits.x > 3
                        event1.hitclusters = clusters
                    end
                end
                if save_clusters
                    event1.clusterpositions = clpositions
                    event1.clustersigmas = clsigmas
                    event1.clusterenergies = clenergies
                    event1.clusterpulses = clpulses
                    event1.clusterepaths = clepaths
                    event1.clusterhpaths = clhpaths
                end
                push!(allevents, event1)
            end
        end #loop over tree
    end #loop over root file

    allevents

end


function plot_crystal_3d(setup::MJDSigGen.Struct_MJD_Siggen_Setup)
    R = setup.xtal_radius
    L = setup.xtal_length
    LT = setup.outer_taper_length
    WT = setup.outer_taper_width
    cyl_x = cos.(linspace(0,2π,40))
    cyl_y = sin.(linspace(0,2π,40))
    cyl_z = zeros(cyl_x)
    myfig = plot(cyl_x*R, cyl_y*R, cyl_z,
        label = "", l = (:grey, :dot),
        xaxis = ("x (mm)", (-R, R)), yaxis = ("y (mm)", (-R, R)), zaxis = ("z (mm)", (0, L)),
        size = (600,600),
    )
    for z in linspace(0,L-LT,5)
        plot!(cyl_x*R, cyl_y*R, cyl_z+z, label = "", linecolor = :grey, line = :dot)
    end
    for z in linspace(L-LT,L,5)
        r = -z*WT/LT + WT/LT*L + R-WT
        plot!(cyl_x*r, cyl_y*r, cyl_z+z, label = "", linecolor = :grey, line = :dot)
    end
    for r in linspace(5,R,5)
        plot!(cyl_x*r, cyl_y*r, cyl_z, label = "", linecolor = :grey, line = :dot)
    end
    for r in linspace(5,R-WT,5)
        plot!(cyl_x*r, cyl_y*r, cyl_z+L, label = "", linecolor = :grey, line = :dot)
    end
    myfig
end

#function to write simple root tree with energy and amplitude values
function write_aoe_tree(allevents::Array{Event}, root_file_name)

    bindings = TTreeBindings()
    Energy = bindings[:Energy] = Ref(zero(Float64))
    Amplitude = bindings[:Amplitude] = Ref(zero(Float64))

    open(TFile, root_file_name, "recreate") do tfile
        ttree = create_ttree!(tfile, "simulation", "Simulated events")
        output = TTreeOutput(ttree, bindings)
        for i in 1:length(allevents)
            Energy.x = allevents[i].totenergy
            Amplitude.x = maximum(vcat([0],diff(allevents[i].sumpulse)))
            push!(output)
        end
    end

end
