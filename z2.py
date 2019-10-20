#!/usr/bin/env python3
# 4D Z2 gauge theory
# basic python adaptation of Creutz's C code
# by Pok Man Lo
# email: pmlo@gsi.de
# 2019 Oct

import numpy as np
from itertools import product

N = 20
latt = np.ones(N**4*4).reshape([N, N, N, N, 4])


# utility
def moveup(xvec, d):
    #  xvec is mutable
    xvec[d] += 1
    xvec[d] = xvec[d] % N   # stays between [0, 1, ..., N-1]
    return None


def movedown(xvec, d):
    #  xvec is mutable
    xvec[d] -= 1
    xvec[d] = xvec[d] % N
    return None


def coldstart():
    _sitegen = [range(N), range(N), range(N), range(N), range(4)]
    for site in product(*_sitegen):
        latt[site] = 1
    return None


def randomstart():
    _sitegen = [range(N), range(N), range(N), range(N), range(4)]
    for site in product(*_sitegen):
        spin = np.random.randint(0, 2)
        if spin == 0:
            spin = -1
        latt[site] = spin
    return None


def update(beta):
    _sitegen = [range(N), range(N), range(N), range(N), range(4)]
    action = 0.
    for site in product(*_sitegen):
        *x, d = site

        # following M. Creutz
        # staples around link(1->4)
        #    dperp        6--5
        #    ^            |  |
        #    |            1--4
        #    |            |  |
        #    -----> d     2--3

        staplesum = 0.
        for dperp in range(4):
            if dperp != d:

                # plaquette 1234
                movedown(x, dperp)
                staple = latt[x[0], x[1], x[2], x[3], dperp]
                staple *= latt[x[0], x[1], x[2], x[3], d]
                moveup(x, d)
                staple *= latt[x[0], x[1], x[2], x[3], dperp]
                moveup(x, dperp)
                staplesum += staple

                # plaquette 4561
                staple = latt[x[0], x[1], x[2], x[3], dperp]
                moveup(x, dperp)
                movedown(x, d)
                staple *= latt[x[0], x[1], x[2], x[3], d]
                movedown(x, dperp)
                staple *= latt[x[0], x[1], x[2], x[3], dperp]
                staplesum += staple

        # calculate the Boltzmann weight
        bplus = np.exp(beta*staplesum)
        bminus = np.exp(-beta*staplesum)
        bplus = bplus/(bplus+bminus)
        # the heatbath algorithm
        r = np.random.random()
        if r < bplus:
            latt[x[0], x[1], x[2], x[3], d] = 1
            action += staplesum
        else:
            latt[x[0], x[1], x[2], x[3], d] = -1
            action -= staplesum

    return 1. - action/N**4/4./6.


def main():

    # beta_c = 0.44

    beta_arr1 = [1. - i1*0.01 for i1 in range(100)]
    beta_arr2 = [0. + i1*0.01 for i1 in range(100)]

    coldstart()
    print('# cold -> hot')
    for beta in beta_arr1:
        action = update(beta)
        print(beta, action)

    # already hot
    print('# hot -> down')
    for beta in beta_arr2:
        action = update(beta)
        print(beta, action)

    return None


main()
