using MGDO
using ROOTFramework
using GERDAMetadata

include("WaveformFunctions.jl")

averagetaus = Dict(0=>199612.0, 1=>197587.1, 2=>192841.1, 3=>186923.3, 4=>184649.8, 5=>193312.5, 6=>196179.7, 7=>183324.6, 8=>187796.1, 9=>193977.3, 10=>185260.4, 11=>189115.3, 12=>195378.1, 13=>177770.6, 14=>180358.2, 15=>201432.6, 16=>198417.7, 17=>198128.2, 18=>213583.7, 19=>178823.5, 20=>187351.5, 21=>187275.5, 22=>192617.9, 23=>209490.7, 24=>187144.8, 25=>203311.2, 26=>186445.6, 27=>188655.0, 28=>191505.9, 29=>207221.9, 30=>209026.0, 31=>189175.2, 32=>198745.6, 33=>203861.1, 34=>192164.4, 35=>203072.2, 36=>199879.0, 37=>208351.8, 38=>184272.6, 39=>187607.0)
nchannels = length(averagetaus)

function process_filekey(gedata, key)

    #input
    filename_tier1 = data_filename(gedata, key, :ged, :tier1)

    #output
    filename_tierA = data_filename(gedata, key, :all, :tierA)
    mkpath(dirname(filename_tierA))
    bindings = TTreeBindings()
    tau1 = bindings[:tau1] = zeros(Float32, nchannels)
    tau2 = bindings[:tau2] = zeros(Float32, nchannels)
    dc1 = bindings[:dc1] = zeros(Float32, nchannels)
    dc2 = bindings[:dc2] = zeros(Float32, nchannels)
    dc3 = bindings[:dc3] = zeros(Float32, nchannels)

    #open input file
    open(MGTEventTree{JlMGTEvent}, filename_tier1) do tier1_tree

        #open output file
        open(TFile, filename_tierA, "recreate") do tfile

            ttree = create_ttree!(tfile, "tierA", "Alpha parameters from Anna")
            output = TTreeOutput(ttree, bindings)

            #loop over events
            for i in eachindex(tier1_tree)

                event = tier1_tree[i]

                #loop over channels
                for chidx in eachindex(event.digitizer_data.ch)
                    ch = event.digitizer_data.ch[chidx]

                    wf = event.waveforms.samples[chidx]
                    bl = mean(wf[1:1000])
                    wf = -(wf-bl)

                    dt = event.waveforms.delta_t[chidx]
                    trigLF = find_intersect(wf, 0.5 * maximum(wf), 5)

                    tau1[ch+1] = 0
                    tau2[ch+1] = 0
                    if trigLF < 2440
                        tau1[ch+1] = get_tau(wf, dt, trigLF+404, trigLF+1654)
                        tau2[ch+1] = get_tau(wf, dt, trigLF+4, trigLF+204)
                    end

                    wf2 = deconv_tau(wf, averagetaus[ch]/dt)
                    normalization = mean(wf2[3097:4096])
                    wf2 = wf2./normalization
                    dc1[ch+1] = 0
                    dc2[ch+1] = 0
                    dc3[ch+1] = 0
                    if trigLF > 0 && trigLF < 3830
                        dc1[ch+1] = mean(wf2[trigLF+10:trigLF+110])
                        dc2[ch+1] = mean(wf2[trigLF+60:trigLF+160])
                        dc3[ch+1] = mean(wf2[trigLF+160:trigLF+260])
                    end

                end #channel loop

                #fill output tree with event
                push!(output)

            end #event loop

        end #open output file

    end #open input file

end #function

gedata = GERDAData(read(DataConfig, "/remote/ceph/user/a/azsigmon/gerda-analysis/Tier1analysis/prod/gerda-dataflow-config.json"))
metadata_dir = meta_data_location(gedata, :gerda)

ds = read(DataSet, joinpath(metadata_dir, "data-sets/phy/run0053-run0079-phy-analysis.txt"))

for key in ds.keys
#key = first(ds.keys)
    println("Processing key ", key)
    process_filekey(gedata, key)
end

