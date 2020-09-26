using Plots

F(x; a) = a*x*(1-x)

xlist = []
ylist = []

function main()

    a_min = 3.8
    a_max = 4.
    na = 1000
    da = (a_max - a_min) / na

    TRANS = 1000
    NPLOT = 100

    for ia = 0:na
        a = a_min + ia*da
        x = [0.5]
        for _ = 1:TRANS
            x[1] = F(x[1]; a=a)
        end
        for _ = 1:NPLOT
            x[1] = F(x[1]; a=a)
            append!(xlist, a)
            append!(ylist, x[1])
        end
    end

    plt = plot(
        xlist, ylist,
        st=:scatter, markersize=1, markercolor=:red,
        ylims=(0.0, 1.2))
    # png(plt, "logistic.png")

end


main()
