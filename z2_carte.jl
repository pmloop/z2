function main(latt)

    # Z(2) gauge theory in 4D (improved)
    # using CartesianIndex
    N = size(latt, 1)
    links = CartesianIndices(latt)

    # utility
    function mvup(x::CartesianIndex, d)
        # tuple -> tuple
        function _f(i1::Int)
            if i1 != d
                return x[i1]
            elseif i1 == d && x[i1] != N
                return x[i1] + 1
            else
                return 1
            end
        end
        return CartesianIndex(ntuple(_f, 4))
    end

    function mvdown(x::CartesianIndex, d)
        # tuple -> tuple
        function _f(i1::Int)
            if i1 != d
                return x[i1]
            elseif i1 == d && x[i1] != 1
                return x[i1] - 1
            else
                return N
            end
        end
        return CartesianIndex(ntuple(_f, 4))
    end

    function coldstart!()
        latt .= 1
    end

    function randomstart!()
        for link in links
            rand(0:1) == 0 ? latt[link] = -1 : latt[link] = 1
        end
    end

    function staplecal(link::CartesianIndex)

        x = CartesianIndex(link.I[1:4])
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
            bplus = exp(beta * staplesum * 1.0)
            bplus /= (bplus + 1.0 / bplus)

            # the heatbath algorithm
            if rand() < bplus
                latt[link] = 1
                action += staplesum
            else
                latt[link] = -1
                action -= staplesum
            end
        end

        return 1.0 - action / N^4 / 4.0 / 6.0
    end

    # beta_c = 0.44
    beta_arr1 = range(1, stop = 0, length = 100)
    beta_arr2 = range(0, stop = 1, length = 100)

    coldstart!()
    # randomstart()
    println("#cold -> hot")
    for beta in beta_arr1
        println([beta, update(beta)])
    end

    # already hot
    println("#hot -> cold")
    for beta in beta_arr2
        println([beta, update(beta)])
    end

end


const N = 40
latt1 = ones(Int64, (N, N, N, N, 4))
@time main(latt1)
