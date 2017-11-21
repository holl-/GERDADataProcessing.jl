using MGDO
using ROOTFramework
using GERDAMetadata

include("WaveformFunctions.jl")

#average of taus measured on the last 1800 samples of the LF waveforms in the full physics data run0053-run0079
#in ns
averagetaus = Dict(0=>199612.0, 1=>197587.1, 2=>192841.1, 3=>186923.3, 4=>184649.8, 5=>193312.5, 6=>196179.7, 7=>183324.6, 8=>187796.1, 9=>193977.3, 10=>185260.4, 11=>189115.3, 12=>195378.1, 13=>177770.6, 14=>180358.2, 15=>201432.6, 16=>198417.7, 17=>198128.2, 18=>213583.7, 19=>178823.5, 20=>187351.5, 21=>187275.5, 22=>192617.9, 23=>209490.7, 24=>187144.8, 25=>203311.2, 26=>186445.6, 27=>188655.0, 28=>191505.9, 29=>207221.9, 30=>209026.0, 31=>189175.2, 32=>198745.6, 33=>203861.1, 34=>192164.4, 35=>203072.2, 36=>199879.0, 37=>208351.8, 38=>184272.6, 39=>187607.0)

reshape_arrays(a::Array, len) = reshape(a, len, div(length(a), len))
reshape_arrays(d::Dict, len) = Dict(ch => reshape_arrays(a, len) for (ch, a) in d)

function process_filekey(key, events, waveforms=0, waveforms_hf=0)

    filename_tier1 = data_filename(gedata, key, :ged, :tier1)
    filename_tier3 = data_filename(gedata, key, :all, :tier3)
    filename_tier4 = data_filename(gedata, key, :all, :tier4)

    tier3_bindings = TTreeBindings()
    rawEnergyGauss = tier3_bindings[:rawEnergyGauss] = zeros(Float64, 0)
    rawEnergyZAC = tier3_bindings[:rawEnergyZAC] = zeros(Float64, 0)
    risetime = tier3_bindings[:risetime] = zeros(Float64, 0)
    #pileup = tier3_bindings[:failedFlag_is0vbbFromCal] = zeros(Int16, 0)

    tier4_bindings = TTreeBindings()
    energy = tier4_bindings[:energy] = zeros(Float64, 0)
    isAoEevaluated = tier4_bindings[:isAoEevaluated] = zeros(Bool, 0)
    AoEclassifier = tier4_bindings[:AoEclassifier] = zeros(Float64, 0)
    isPSDVetoed = tier4_bindings[:isPSDVetoed] = zeros(Int32, 0)
    isTP = tier4_bindings[:isTP] = Ref(zero(Int32))
    isBL = tier4_bindings[:isBL] = Ref(zero(Int32))
    multiplicity = tier4_bindings[:multiplicity] = Ref(zero(Int32))
    isMuVetoed = tier4_bindings[:isMuVetoed] = Ref(zero(Int32))
    isLArVetoed = tier4_bindings[:isLArVetoed] = Ref(zero(Int32))

    open(MGTEventTree{JlMGTEvent}, filename_tier1) do tier1_tree
        open(TChainInput, tier3_bindings, "tier3", filename_tier3) do tier3_tree
            open(TChainInput, tier4_bindings, "tier4", filename_tier4) do tier4_tree
                for i in eachindex(tier1_tree)

                    event = tier1_tree[i]
                    tier3_tree[i]
                    tier4_tree[i]

                    if isBL.x==1 || isTP.x==1 || multiplicity.x!=1
                        continue
                    end
                    if isMuVetoed.x==1
                        continue
                    end

                    for chidx in eachindex(event.digitizer_data.ch)
                        ch = event.digitizer_data.ch[chidx]
                        det = ch+1
                        if energy[det]>1000 && energy[det]<9999

                            wf = event.waveforms.samples[chidx]
                            bl = mean(wf[1:1000])
                            wf = -(wf-bl)
                            if waveforms != 0
                                append!(get!(() -> Float32[], waveforms, ch), wf)
                            end

                            dt = event.waveforms.delta_t[chidx]
                            trigLF = find_intersect(wf, 0.5 * maximum(wf), 2)
                            if trigLF < 2440
                                tau1 = get_tau(wf, dt, trigLF+404, trigLF+1654)
                                tau2 = get_tau(wf, dt, trigLF+4, trigLF+204)
                            else
                                tau1 = 0
                                tau2 = 0
                            end

                            wf2 = deconv_tau(wf, averagetaus[ch]/dt)
                            normalization = mean(wf2[3097:4096])
                            wf2 = wf2./normalization
                            charge1 = 0
                            charge2 = 0
                            charge3 = 0
                            if trigLF > 0 && trigLF < 3830
                                charge1 = mean(wf2[trigLF+10:trigLF+110])
                                charge2 = mean(wf2[trigLF+60:trigLF+160])
                                charge3 = mean(wf2[trigLF+160:trigLF+260])
                            end
                            
                            awf = event.aux_waveforms.samples[chidx]
                            bl = mean(awf[1:250])
                            awf = -(awf-bl)
                            if waveforms_hf != 0
                                append!(get!(() -> Float32[], waveforms_hf, ch), awf)
                            end

                            adt = event.aux_waveforms.delta_t[chidx]
                            trigHF = find_intersect(awf, 0.5*maximum(awf), 3)
                            rt0595 = get_risetime(awf, 0.05, 0.95, 3, adt)
                            rt0550 = get_risetime(awf, 0.05, 0.50, 3, adt)
                            rt5095 = get_risetime(awf, 0.50, 0.95, 3, adt)
                            if trigHF < 500
                                tau3 = get_tau(awf, adt, trigHF+20, trigHF+500)
                            else
                                tau3 = 0
                            end

                            AoE = isAoEevaluated[det] ? AoEclassifier[det] : 99

                            ev = [rawEnergyGauss[det] rawEnergyZAC[det] risetime[det] isLArVetoed.x AoE isPSDVetoed[det] tau1 tau2 tau3 rt0595 rt0550 rt5095 normalization charge1 charge2 charge3] 
                            append!(get!(() -> Float32[], events, ch), ev)

                        end
                    end #channel loop
                end #event loop
            end #open
        end #open
    end #open
    
end #function

