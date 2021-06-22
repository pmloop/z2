function main()

N = 20
latt = Array{Int64}(undef, N, N, N, N, 4)
# latt .= 1
# latt = ones(Int64, (N, N, N, N, 4))

# utility
function moveup(xvec, d)
    #  xvec is mutable
    if xvec[d] == N
        xvec[d] = 1
    else
        xvec[d] += 1
    end
    # xvec[d] = xvec[d] % N
end


function movedown(xvec, d)
    #  xvec is mutable
    if xvec[d] == 1
        xvec[d] = N
    else
        xvec[d] -= 1
    end
end


function coldstart()
    latt .= 1
end


function randomstart()

    sites = collect(Base.product(1:N, 1:N, 1:N, 1:N, 1:4))

    for site in sites
        x1, x2, x3, x4, d = site

        spin = rand(0:1)
        if spin == 0
            spin = -1
        end

        latt[x1, x2, x3, x4, d] = spin
    end

end


function update(beta)

    sites = collect(Base.product(1:N, 1:N, 1:N, 1:N, 1:4))

    action = 0.

    for site in sites

        x1, x2, x3, x4, d = site
        # print("$x")
        x = [x1, x2, x3, x4]
        # d = site[5]

        # following M. Creutz
        # staples around link(1->4)
        #    dperp        6--5
        #    ^            |  |
        #    |            1--4
        #    |            |  |
        #    -----> d     2--3

        staplesum = 0.
        for dperp in 1:4
            if dperp != d

                # plaquette 1234
                movedown(x, dperp)

                # print("here")

                staple = latt[x[1], x[2], x[3], x[4], dperp]
                staple *= latt[x[1], x[2], x[3], x[4], d]
                moveup(x, d)
                staple *= latt[x[1], x[2], x[3], x[4], dperp]
                moveup(x, dperp)
                staplesum += staple

                # plaquette 4561
                staple = latt[x[1], x[2], x[3], x[4], dperp]
                moveup(x, dperp)
                movedown(x, d)
                staple *= latt[x[1], x[2], x[3], x[4], d]
                movedown(x, dperp)
                staple *= latt[x[1], x[2], x[3], x[4], dperp]
                staplesum += staple
            end

        end

        # calculate the Boltzmann weight
        bplus = exp(beta*staplesum)
        bminus = exp(-beta*staplesum)
        bplus = bplus/(bplus+bminus)

        # the heatbath algorithm
        r = rand()

        if r <= bplus
            latt[x1, x2, x3, x4, d] = 1
            action += staplesum
        else
            latt[x1, x2, x3, x4, d] = -1
            action -= staplesum
        end
    end

    return 1. - action / N^4 / 4. / 6.
end


    # beta_c = 0.44
    beta_arr1 = range(1., stop=0., length=100)
    beta_arr2 = range(0., stop=1., length=100)

    coldstart()
    # randomstart()
    println("cold -> hot")
    for beta in beta_arr1
        action = update(beta)
        println("$beta $action")
    end

    # already hot
    println("hot -> cold")
    for beta in beta_arr2
        action = update(beta)
        println("$beta $action")
    end

end


main()
