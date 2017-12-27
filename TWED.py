'''
Created on Oct 19, 2017
@author: Qiangeng Xu
'''

import math


def init_matrix(data):
    for i in xrange(len(data)):
        data[i][0] = float('inf')
    for i in xrange(len(data[0])):
        data[0][i] = float('inf')
    data[0][0] = 0
    return data


def LpDist(time_pt_1, time_pt_2):
    if (type(time_pt_1) == int and type(time_pt_2) == int):
        return abs(time_pt_1 - time_pt_2)
    else:
        return sum(abs(time_pt_1 - time_pt_2))


def TWED(t1_data, t2_data, lam, nu):
    """"Requires: t1: multivariate time series in numpy matrix format.
    t2: multivariate time series in numpy matrix format. lam: penalty lambda parameter, nu: stiffness coefficient"""
    # """Returns the TWED distance between the two time series. """
    n = t1_data.shape[0]
    m = t2_data.shape[0]
    result = [[0] * m for row in xrange(n)]
    result = init_matrix(result)
    for i in xrange(1, n):
        for j in xrange(1, m):
            insertion = result[i - 1][j] + LpDist(t1_data[i - 1], t1_data[i]) + \
                         nu  + lam
            deletion = result[i][j - 1] + LpDist(t2_data[j - 1], t2_data[j]) + \
                        nu + lam
            # print i, j, n , m, t1_time[i], t2_time[j]
            match = result[i - 1][j - 1] + LpDist(t1_data[i], t2_data[j]) + \
                     2 * nu * (abs(i - j)) + LpDist(t1_data[i - 1], t2_data[j - 1])
            # print nu, lam, insertion, deletion, match
            result[i][j] = min(insertion, deletion, match)
    return result[n - 1][m - 1]
