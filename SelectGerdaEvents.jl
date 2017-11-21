using MGDO
using ROOTFramework
using GERDAMetadata
using JLD

reshape_arrays(a::Array, len) = reshape(a, len, div(length(a), len))
reshape_arrays(d::Dict, len) = Dict(ch => reshape_arrays(a, len) for (ch, a) in d)

function process_filekey(key, events, waveforms=0, waveforms_hf=0)

    filename_tier1 = data_filename(gedata, key, :ged, :tier1)
    filename_tier4 = data_filename(gedata, key, :all, :tier4)

    tier4_bindings = TTreeBindings()
    energy = tier4_bindings[:energy] = zeros(Float64, 0)
    isAoEevaluated = tier4_bindings[:isAoEevaluated] = zeros(Bool, 0)
    AoEclassifier = tier4_bindings[:AoEclassifier] = zeros(Float64, 0)
    isTP = tier4_bindings[:isTP] = Ref(zero(Int32))
    isBL = tier4_bindings[:isBL] = Ref(zero(Int32))
    multiplicity = tier4_bindings[:multiplicity] = Ref(zero(Int32))

    open(MGTEventTree{JlMGTEvent}, filename_tier1) do tier1_tree
        open(TChainInput, tier4_bindings, "tier4", filename_tier4) do tier4_tree
            for i in eachindex(tier1_tree)

                event = tier1_tree[i]
                tier4_tree[i]

                if isBL.x==1 || isTP.x==1 || multiplicity.x!=1
                    continue
                end

                for chidx in eachindex(event.digitizer_data.ch)
                    ch = event.digitizer_data.ch[chidx]
                    det = ch+1
                    if energy[det]>1591 && energy[det]<1594

                        wf = event.waveforms.samples[chidx]
                        bl = mean(wf[1:1000])
                        wf = -(wf-bl)
                        if waveforms != 0
                            append!(get!(() -> Float32[], waveforms, ch), wf)
                        end

                        awf = event.aux_waveforms.samples[chidx]
                        bl = mean(awf[1:250])
                        awf = -(awf-bl)
                        if waveforms_hf != 0
                            append!(get!(() -> Float32[], waveforms_hf, ch), awf)
                        end

                        AoE = isAoEevaluated[det] ? AoEclassifier[det] : 99

                        ev = [energy[det] AoE] 
                        append!(get!(() -> Float32[], events, ch), ev)

                    end
                end #channel loop
            end #event loop
        end #open
    end #open
    
end #function

gedata = GERDAData(read(DataConfig, "/remote/ceph/group/gerda/data/phase2/blind/v03.00/gerda-dataflow-config.json"))
metadata_dir = meta_data_location(gedata, :gerda)

ds = read(DataSet, joinpath(metadata_dir, "data-sets/cal/run0053-run0064-cal-analysis.txt"))

events = Dict{Int, Vector{Float32}}()
waveforms = Dict{Int, Vector{Float32}}()
waveforms_hf = Dict{Int, Vector{Float32}}()

for key in ds.keys
    println(key)
    process_filekey(key, events, waveforms, waveforms_hf)
end

events = reshape_arrays(events, 2)
waveforms = reshape_arrays(waveforms, 4096)
waveforms_hf = reshape_arrays(waveforms_hf, 1000)

save("selected_DEP_events_run0053-run0064-cal-analysis.jld", "events", events, "waveforms", waveforms, "waveforms_hf", waveforms_hf)

