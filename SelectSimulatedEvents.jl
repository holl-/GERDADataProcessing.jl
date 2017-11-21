using JLD

include("/remote/ceph/user/a/azsigmon/Simulation/SigGen/SimulatedEvents.jl")
include("/remote/ceph/user/a/azsigmon/Simulation/Electronics/WaveformFunctions.jl")

function readfiles(detName::AbstractString, n::Integer = 99)

    setup = signal_calc_init("/remote/ceph/user/a/azsigmon/Simulation/SigGen/config_files/$(detName).config")
    treename = "eventTree_$(detName)"
    minen = 1592
    maxen = 1593
    clusteringsize = 0.3

    allevents = Array{Event}(0)

    for i in 0:n
        filename = "/remote/ceph/group/gerda/data/simulation/converted-hits/$(detName)/hittree_$(detName)_calibrations_$(i).root"
        events = create_events(setup, filename, treename, min_energy = minen, max_energy = maxen, cluster_size = clusteringsize)
        allevents = vcat(allevents, events)
    end

    allevents

end

function processchannel(detName::AbstractString, outDir::AbstractString)

    allevents = readfiles(detName, 99)
    println(detName, " events read and pulses created ...")

    waveforms = zeros(Float32, 1000, length(allevents))

    for i in 1:length(allevents)

        wf = allevents[i].sumpulse
        wf = resample_wf(wf, 10)
        if length(wf) != 200
            error("Error: Waveform length does not fit expectation")
        end
        trig = find_intersect(wf, 0.5*wf[length(wf)])
        if trig > 0
            waveforms[:, i] = vcat(zeros(Float32, 400-trig), wf, ones(Float32, 400+trig).*allevents[i].totenergy)
        else
            waveforms[:, i] = vcat(zeros(Float32, 400), wf, ones(Float32, 400).*allevents[i].totenergy)
        end

    end

    save("$outDir/selected_DEP_events_simulation_$detName.jld", "waveforms", waveforms)
    println(detName, " file saved ...")

end

