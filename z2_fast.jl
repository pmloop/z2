function main(N = 50)

    # Z(2) gauge theory in 4D (improved)
    latt = ones(Int64, (N, N, N, N, 4))
    links = Iterators.product(1:N, 1:N, 1:N, 1:N, 1:4)

    # utility
    function mvup(x, d)
        # tuple -> tuple
        f1(i1) = x[i1]
        f2(i1) = x[i1] == N ? 1 : x[i1] + 1
        f(i1) = i1 != d ? f1(i1) : f2(i1)
        return ntuple(f, 4)
    end
    function mvdown(x, d)
        # tuple -> tuple
        f1(i1) = x[i1]
        f2(i1) = x[i1] == 1 ? N : x[i1] - 1
        f(i1) = i1 != d ? f1(i1) : f2(i1)
        return ntuple(f, 4)
    end

    function coldstart!()
        latt .= 1
    end

    function randomstart!()
        latt .= rand([-1, 1])
    end

    function staplecal(link)

        # x is a tuple
        x = link[1:4]
        d = link[5]

        # following M. Creutz
        # staples around link(1->4)
        #    dperp        6--5
        #    ^            |  |
        #    |            1--4
        #    |            |  |
        #    -----> d     2--3

        staplesum = 0.0
        for dperp = 1:4
            if dperp != d

                # plaquette 1234
                x = mvdown(x, dperp)
                # x is tuple, and hence x... is fast (forum question)
                link1 = latt[x..., dperp]
                link2 = latt[x..., d]
                x = mvup(x, d)
                link3 = latt[x..., dperp]
                x = mvup(x, dperp)
                # plaquette 4561
                link4 = latt[x..., dperp]
                x = mvup(x, dperp)
                x = mvdown(x, d)
                link5 = latt[x..., d]
                x = mvdown(x, dperp)
                link6 = latt[x..., dperp]

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
            bplus = exp(beta * staplesum)
            bplus /= (bplus + 1.0 / bplus)

            # the heatbath algorithm
            if rand() < bplus
                latt[link...] = 1
                action += staplesum
            else
                latt[link...] = -1
                action -= staplesum
            end
        end

        return 1.0 - action / N^4 / 4.0 / 6.0
    end

    # beta_c = 0.44
    beta_arr1 = range(1, stop = 0, length = 50)
    beta_arr2 = range(0, stop = 1, length = 50)

    coldstart!()
    # randomstart()
    println("#cold -> hot")
    for beta in beta_arr1
        println([beta, " ", update(beta)]...)
    end

    # already hot
    println("#hot -> cold")
    for beta in beta_arr2
        println([beta, " ", update(beta)]...)
    end

end

@time main()
