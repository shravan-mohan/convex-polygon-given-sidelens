import numpy as np
import cvxpy as cvx
from sympy.geometry import Segment
import matplotlib.pyplot as plt

r = np.random.rand(10)

def checkValidity(r):
    """
    This function checks if a polygon exists with the given side lengths.
    :param r: Side lengths as a vector.
    :return: True or False.
    """
    s = np.sum(r)/2
    for k in range(len(r)):
        if(s<r[k]):
            return False
    return True

def checkDists(X, r):
    """
    This function checks if the distance between consecutive points are in agreement
    with the given side lengths.
    :param X: Coordinate matrices of the polygon.
    :param r: Side lengths as a vector.
    :return: True or False.
    """
    dists = []
    for k in range(len(r)-1):
        dists = dists + [np.linalg.norm(X[k,:]-X[k+1,:])]
    dists = dists + [np.linalg.norm(X[-1,:]-X[0,:])]

    if(np.linalg.norm(dists-r)<=1e-6):
        return True
    else:
        print('Something went wrong!')


def checkConvexity(A, r):
    """
    This function checks if the given polygon with coordinates is convex.
    :param A: Coordinate matrices of the polygon.
    :param r: Side lengths as a vector.
    :return: True or False
    """
    An = np.vstack((A, A[0,:], A[1,:]))
    t = []
    for k in range(len(r)):
        t = t + [np.cross(An[k+1,:]-An[k,:], An[k+2,:]-An[k+1,:])]

    if((np.array(t)<=0).all() or (np.array(t)>=0).all()):
        return True
    else:
        print('Something went wrong!')

def checkSimple(A,r):
    """
    This function checks if the given polygon with coordinates is simple.
    :param A: Coordinate matrices of the polygon.
    :param r: Side lengths as a vector.
    :return: True or False
    """
    n = len(r)
    An = np.vstack((A, A[0,:], A[1,:]))
    flag = 0
    for k in range(n):
        if(k<n-1):
            for l in range(max(0,k-2)):
                a = Segment(tuple(An[l, :]), tuple(An[l + 1, :]))
                b = Segment(tuple(An[k, :]), tuple(An[k + 1, :]))
                pts = a.intersect(b)
                if (len(pts) == 1):
                    flag = 1
                    break
            if (flag == 1):
                break
        else:
            for l in range(1, max(0,k-2)):
                a = Segment(tuple(An[l, :]), tuple(An[l + 1, :]))
                b = Segment(tuple(An[k, :]), tuple(An[k + 1, :]))
                pts = a.intersect(b)
                if (len(pts) == 1):
                    flag = 1
                    break
            if (flag == 1):
                break

    if(flag==1):
        print('The polygon is not simple!')
        return False
    else:
        return True

def getConvexPolygon(r=np.array([0.43202382, 0.53261951, 0.24482898, 0.37902022, 0.31978369,
                                0.76482198, 0.53154303, 0.22408833, 0.03323694, 0.63061874]),
                     solver='ECOS'):
    """
    This function generates the coordinates of the polygon with the given side lengths.
    :param r: Side lengths as a vector.
    :param solver: One of the solvers in CVXPY. Default set to 'ECOS'.
    :return: The coordinates matrix of the generated polygon and its side lengths.
    """
    if(checkValidity(r)==False):
        print('Such a polygon does not exist!')
        return -1

    X = cvx.Variable((len(r), 2))

    constraints = []
    for k in range(len(r) - 1):
        constraints = constraints + [cvx.norm2(X[k, :] - X[k + 1, :]) <= r[k]]
    constraints = constraints + [cvx.norm2(X[-1, :] - X[0, :]) <= r[-1]]
    constraints = constraints + [np.ones(10) @ X == 0]
    
    Y = np.zeros(2)
    theta = 2 * np.pi / len(r)
    for k in range(1, len(r)):
        Y = np.vstack((Y, np.array([np.cos(k * theta), np.sin(k * theta)])))

    prob = cvx.Problem(cvx.Maximize(cvx.trace(X @ (Y.T))), constraints)
    prob.solve(solver=solver)

    if(prob.status=='optimal_inaccurate'):
        print('The solution is inaccurate. Please try with another solver!')

    plotPoly(X.value)
    return X.value, r

def plotPoly(A):
    """
    This function plots a polygon with the given coordinates.
    :param A: Coordinate matrices of the polygon.
    :return: Plot of the polygon.
    """
    coord = [[A[k,0],A[k,1]] for k in range(A.shape[0])]
    coord.append(coord[0])
    xs, ys = zip(*coord)

    plt.figure()
    plt.plot(xs, ys)
    plt.scatter(xs, ys)
    plt.show()

if(__name__=='__main__'):
    A, r = getConvexPolygon()
    print(checkDists(A, r))
    print(checkSimple(A, r))
    print(checkConvexity(A, r))








