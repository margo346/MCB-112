#! /usr/bin/env python3

# Usage:
#   w09-visualize.py <infile> <fig.png>
#
# Example:
#   w09-visualize.py w09-data.tbl foo.png
#
# Reads w09-data.tbl; produces a figure in <fig.png> that visualizes the
# five true clusters.

import matplotlib.pyplot as plt
import numpy             as np
import sys

infile = sys.argv[1]
outfig = sys.argv[2]


def read_data(infile):
    '''
    read_data(infile)
    Read Wiggins' input file, w09-data.tbl, or a file in that format.
    Return:
       ctype[0..N-1] : cell types 0..Q-1 for each cell i
       data[i,g]     : array of count data; rows = cells i; cols = genes g
       N             : number of cells (rows of the data file)
       G             : number of genes (cols of the data file, after the first)
       Q             : number of cell types
    '''
    ctype = []
    data  = []
    with open(infile) as f:
        for line in f:
            if line[0] == '#': continue   # skip comment lines
            fields = line.split()
            ctype.append(int(fields[1]))
            data.append( [int(fields[2]), int(fields[3])])  # assumes exactly 2 genes!!
    ctype = np.array(ctype)
    data  = np.array(data)
    N, G  = np.shape(data)
    Q     = np.max(ctype) + 1
    return ctype, data, N, G, Q


def initialize_at_true():
    '''
    initialize_at_true():
    Returns the true mu centroids, and the true proportions;
    don't say I never gave you anything.
       mu[q,g]  : array of means for mixture q, gene g
       qp[q]    : mixture coefficient for mixture q
    '''
    qp = np.array([ 0.1, 0.2, 0.4, 0.2, 0.1 ])
    mu = np.array([[    3., 100. ],
                   [  100., 100. ],
                   [   30.,  30. ],
                   [    3.,   3. ],
                   [  100.,   3. ]])
    return mu, qp

def visualize_data(data, mu, C, outpng):
    '''
    visualize_data():

    This might give you a starting point that saves some matplotlib
    machinations; you can certainly spiff this up from here.

    Input:
       data[i,g] : count data for each cell i, for each gene g
       mu[q,g]   : array of mean counts for mixture q, gene g
       C[i]      : assignment of cell i to a cluster 0..Q-1
       outpng    : save figure to PNG file (must end in .png; example 'foo.png')

    '''
    N, G  = np.shape(data)
    Q, G2 = np.shape(mu)
    assert G == G2
    assert len(C) == N

    genes = [ 'defA', 'kilA' ]

    # We can assign colors to up to Q=10 components. If you want more, add more.
    colormap = ['xkcd:orange', 'xkcd:olive',     'xkcd:azure',    'xkcd:rose', 'xkcd:mustard', 
                'xkcd:peach',  'xkcd:turquoise', 'xkcd:lavender', 'xkcd:rust', 'xkcd:red']

    fig, ax = plt.subplots()
    for i in range(N):
        edgecolor = colormap[ C[i]]
        fillcolor = 'w'
        shape     = 'o'
        ax.loglog( data[i,0], data[i,1], marker=shape, mec=edgecolor, mfc=fillcolor, mew=1.5)

    for q in range(Q):
        ax.loglog(mu[q,0], mu[q,1], '*k', ms=10)

    ax.set_xlabel('{} (counts)'.format(genes[0]))
    ax.set_ylabel('{} (counts)'.format(genes[1]))

    fig.savefig(outpng)


def main():
    ctype, data, N, G, Q = read_data(infile)
    centroids, qp        = initialize_at_true()
    visualize_data(data, centroids, ctype, outfig)

if __name__ == "__main__":
    main()


