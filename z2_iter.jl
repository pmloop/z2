function main(Ns = 20, dim = 4)

    latt = ones(Int64, ntuple(i1 -> i1 <= dim ? Ns : dim, dim + 1))
    links = CartesianIndices(latt)

    function mvup(x::CartesianIndex, d)
        # key routine to optimize: KISS
        f1(i1) = x[i1]
        f2(i1) = x[i1] == Ns ? 1 : x[i1] + 1
        f(i1) = i1 != d ? f1(i1) : f2(i1)
        return CartesianIndex(ntuple(f, length(x)))

    end

    function mvdown(x::CartesianIndex, d)
        # key routine to optimize: KISS
        f1(i1) = x[i1]
        f2(i1) = x[i1] == 1 ? Ns : x[i1] - 1
        f(i1) = i1 != d ? f1(i1) : f2(i1)
        return CartesianIndex(ntuple(f, length(x)))

    end

    function coldstart()
        latt .= 1
    end

    function randomstart()
        latt .= rand([-1, 1])
    end

    function staplecal(link::CartesianIndex)

        x = CartesianIndex(link.I[1:end-1])
        d = link.I[end]

        # following M. Creutz
        # staples around link(1->4)
        #    dperp        6--5
        #    ^            |  |
        #    |            1--4
        #    |            |  |
        #    -----> d     2--3

        staplesum = 0.0
        for dperp = 1:length(x)
            if dperp != d

                # plaquette 1234
                x = mvdown(x, dperp)
                link1 = latt[x, dperp]
                link2 = latt[x, d]
                x = mvup(x, d)
                link3 = latt[x, dperp]
                x = mvup(x, dperp)
                # plaquette 4561
                link4 = latt[x, dperp]
                x = mvup(x, dperp)
                x = mvdown(x, d)
                link5 = latt[x, d]
                x = mvdown(x, dperp)
                link6 = latt[x, dperp]
                staplesum += link1 * link2 * link3 + link4 * link5 * link6

            end
        end
        return staplesum
    end

    function update(beta)

        action = 0.0
        for link in links

            staplesum = staplecal(link)
            # calculate the Boltzmann weight
            bplus = exp(beta * staplesum * 1)
            bplus /= (bplus + 1 / bplus)

            # the heatbath algorithm
            if rand() < bplus
                latt[link] = 1
                action += staplesum
            else
                latt[link] = -1
                action -= staplesum
            end
        end

        return 1.0 - action / Ns^dim / dim / 6.0
    end

    beta_arr1 = range(1, stop = 0, length = 20)
    beta_arr2 = range(0, stop = 1, length = 20)

    coldstart()
    # randomstart()
    println("#cold -> hot")
    for beta in beta_arr1
        println(beta, " ", update(beta))
    end

    # already hot
    println("#hot -> cold")
    for beta in beta_arr2
        println(beta, " ", update(beta))
    end

end

@time main()
