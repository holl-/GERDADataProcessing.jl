using DSP

function applyElectronics(pulse; Ts = 10e-9, GBP = 2750e6, tau = 180e-6, Kv = 150e3, Cd = 50e-12, Cf = 0.35e-12, Rf = 500e6)

    wop = GBP / (2 * pi * Kv)
    Cmod = Cf + Cd
    wmod = 1.0 / (Rf * Cmod)
    alfa = Cmod / (Cf * GBP)

    b0 = 1.0 / alfa
    a2 = 1.0
    a1 = 1.0 / alfa + wop + wmod
    a0 = 1.0 / (tau * alfa) + wmod*wop

    # then the transfer function in the *Laplace* s-domain looks like this:
    #                       b0
    #   T(s) = ----------------------------
    #             a2 * s^2 + a1 * s + a0

    # PolynomialRatio needs z-transform paramters: s- and z-domains can be connected by
    # the bilinear transform:
    #        z - 1
    # s = K -------- , K = Ts/2  , Ts - sampling period
    #        z + 1 
    #        
    # we can then convert T(s) to T(z):
    #              bz2 * z^2 + bz1 * z + bz0
    #   T(z) = -------------------------------
    #              az2 * z^2 + az1 * z + az0
    #

    K = 2/Ts

    az2 = 1.0   # normalized
    az1 = (2*a0 - 2*K^2)/(K^2 + a1*K + a0)
    az0 = (K^2 - a1*K + a0)/(K^2 + a1*K + a0)

    bz2 = b0/(K^2 + a1*K + a0)
    bz1 = 2*b0/(K^2 + a1*K + a0)
    bz0 = b0/(K^2 + a1*K + a0)

    myfilter = PolynomialRatio([bz2, bz1, bz0], [az2, az1, az0])

    filtered = filt(myfilter, vcat([0], diff(pulse)))

end
