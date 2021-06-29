function main()

N = 40
latt = ones(Int64, (N, N, N, N, 4))

# utility
function mvup!(xvec, d)
    #  xvec would be modified
    xvec[d] += 1
    if xvec[d] == N + 1
        xvec[d] = 1
    end
    return nothing
end


function mvdown!(xvec, d)
    #  xvec would be modified
    xvec[d] -= 1
    if xvec[d] == 0
        xvec[d] = N
    end
    return nothing
end


function coldstart()
    latt .= 1
end


function randomstart()

    # col. based
    for d in 1:4,
        i4 in 1:N,
        i3 in 1:N,
        i2 in 1:N,
        i1 in 1:N

        spin = rand(0:1)
        if spin == 0
            spin = -1
        end
        latt[i1, i2, i3, i4, d] = spin

    end
end


function update(beta)

    action = 0.

    # col. based
    for d in 1:4,
        i4 in 1:N,
        i3 in 1:N,
        i2 in 1:N,
        i1 in 1:N

        x = [i1, i2, i3, i4]

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
                mvdown!(x, dperp)
                staple = latt[x[1], x[2], x[3], x[4], dperp]
                staple *= latt[x[1], x[2], x[3], x[4], d]
                mvup!(x, d)
                staple *= latt[x[1], x[2], x[3], x[4], dperp]
                mvup!(x, dperp)
                staplesum += staple

                # plaquette 4561
                staple = latt[x[1], x[2], x[3], x[4], dperp]
                mvup!(x, dperp)
                mvdown!(x, d)
                staple *= latt[x[1], x[2], x[3], x[4], d]
                mvdown!(x, dperp)
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
            latt[i1, i2, i3, i4, d] = 1
            action += staplesum
        else
            latt[i1, i2, i3, i4, d] = -1
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
    println("#cold -> hot")
    for beta in beta_arr1
        action = update(beta)
        println("$beta $action")
    end

    # already hot
    println("#hot -> cold")
    for beta in beta_arr2
        action = update(beta)
        println("$beta $action")
    end

end


@time main()
