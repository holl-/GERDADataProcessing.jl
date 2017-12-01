using LsqFit
using Interpolations

function find_intersect(samples::Vector, threshold::Real, noisefilter = 1)

    if noisefilter < 1 noisefilter = 1 end
    if length(samples) < 1 error("Empty array") end
    i = 1
    x = samples[i]
    counter = 0
    intersect = i
    if x != threshold
        findHigher = x < threshold
        while i < length(samples)
            i += 1
            x = samples[i]
            if (findHigher && (x >= threshold)) || (!findHigher && (x <= threshold))
                if counter == 0
                    intersect = i
                end
                if counter >= noisefilter-1
                    return intersect
                else
                    counter += 1
                end
            else
                counter = 0
            end
        end
        #println("No intersect found")
        return NaN
    end
    return NaN

end


function get_risetime(samples::Vector; lowfraction = 0.05, highfraction = 0.95, pulseheight = 0, samplingtime = 1, noisefilter = 1)

    if pulseheight == 0
        pulseheight = maximum(samples)
    end
    lowindex = find_intersect(samples, lowfraction * pulseheight, noisefilter)
    highindex = find_intersect(samples, highfraction * pulseheight, noisefilter)
    risetime = (highindex - lowindex) * samplingtime

end


function get_risetime_interpolate(samples::Vector; lowfraction = 0.05, highfraction = 0.95, pulseheight = 0, itpfactor = 0.1, samplingtime = 1, noisefilter = 1)

    if pulseheight == 0
        pulseheight = maximum(samples)
    end

    itp = interpolate(wf, BSpline(Quadratic(Flat())), OnCell())
    interpolated_samples = [itp[t] for t in itpfactor:itpfactor:length(samples)]

    lowindex = find_intersect(interpolated_samples, lowfraction * pulseheight, noisefilter)
    highindex = find_intersect(interpolated_samples, highfraction * pulseheight, noisefilter)
    risetime = (highindex - lowindex) * samplingtime * itpfactor

end


function resample_wf(samples::Vector, factor::Integer)
    if rem(length(samples), factor) != 0
        error("Error: resampling does not work with your factor")
    end
    n_out = convert(Int64, length(samples)/factor)
    wf_out = zeros(typeof(samples[1]), n_out)
    for j in 1:n_out
        wf_out[j] = sum(samples[((j-1)*factor+1):(j*factor)])/factor
    end
    wf_out
end


function deconv_tau(samples::Vector, tau::Real)

    alpha = 1 - exp(- 1/tau)
    acc = zero(eltype(samples))
    newsamples = zeros(samples)
    for i in eachindex(samples)
        x = samples[i]
        new_x = x + acc
        acc += x * alpha
        newsamples[i] = new_x
    end
    newsamples

end


function deconv_tau!(samples::Vector, tau::Real)

    samples = deconv_tau(samples, tau)

end


model_exp(x, p) = p[1]*exp.(-x.*p[2])

function get_tau(samples::Vector, samplingtime::Real, firstsample::Integer, lastsample::Integer)

    xdata = convert(Array{Float64}, linspace(firstsample, lastsample, lastsample-firstsample+1))
    ydata = convert(Array{Float64}, samples[firstsample:lastsample])
    fitresult = curve_fit(model_exp, xdata, ydata, [2000., 1/4600.])
    samplingtime / fitresult.param[2]

end


function mwa(samples::Vector, window::Integer)
    newsamples = zeros(samples)
    for i in window:length(samples)
        newsamples[i] = mean(samples[i-window+1:i])
    end
    newsamples
end
