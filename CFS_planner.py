import numpy as np
import time
from CFS_problem import *
import cvxopt
from cvxopt import matrix, solvers
from numpy import linalg as LA
import matplotlib.pyplot as plt
from cvxpy import *
import matplotlib.pyplot as plt



#Quadratic programming solver.
def cvxopt_solve_qp(P, q, G=None, h=None, A=None, b=None):
    P = .5 * (P + P.T)  # make sure P is symmetric
    args = [matrix(P), matrix(q)]
    if G is not None:
        args.extend([matrix(G), matrix(h)])
        if A is not None:
            args.extend([matrix(A), matrix(b)])
    sol = cvxopt.solvers.qp(*args)
    if 'optimal' not in sol['status']:
        return None
    return np.array(sol['x']).reshape((P.shape[1],))



def get_tayler_expansion_points(high_pri_pt, low_pri_pt, mini_distance):
    new_ref = np.zeros((4, 1))
    new_ref[0:2] = high_pri_pt
    dir = (low_pri_pt - high_pri_pt) / np.linalg.norm(low_pri_pt - high_pri_pt)
    new_ref[2:4] = new_ref[0:2] + dir * mini_distance
    return new_ref

def two_car_distance(path1, path2):
    nstep = path1.shape[0]
    distance = []
    for i in range(nstep):
        pts1 = path1[i, :]
        pts2 = path2[i, :]
        dist = np.linalg.norm(pts1-pts2)
        distance.append(dist)
    plt.plot(distance)
    plt.show()
    print(np.min(distance))

def get_line(x1, y1, x2, y2):
    '''
    Args:
        x1, y1: Point 1 on the line.
        x2, y2: Point 2 on the line.
    Return:
        coe: coe[0:2] is the k vector, coe[2] is the b, so that k[0]x+k[1]y+b=0.
    '''
    coe = np.zeros((3,1))
    coe[0] = y2 - y1
    coe[1] = -(x2 - x1)
    coe[2] = y2*(x2-x1) - x2*(y2-y1)
    return coe

    
def Plan_trajectory(MAX_ITER, multi_path, mini_distance):
	#mini_distance = 20
    Qref, Qabs, nstep, dim, oripath, I_2 = Setup_problem(multi_path)
    refpath = oripath
    print("refpath shape is:{}".format(refpath.shape))
    # print("refpath and oripath:{}".format(refpath.squeeze()))
    print("Qref shape:{}, Qabs shape:{}".format(Qref.shape, Qabs.shape))
    Qe = Qref + Qabs
    n = nstep * dim
    print("n is: {}".format(n))

    
    for i in range(MAX_ITER):
        print(i)
        x = Variable(n)
        objective = Minimize(0.5*quad_form(x,Qe) + (np.matmul(-Qref,oripath)).T @ x)
        constraints = []
        for j in range(nstep):
            x_ref_1 = refpath[dim*j:dim*(j+1)]
            if j < nstep-1:
                x_ref_2 = refpath[dim*(j+1):dim*(j+2)]
            if j <= 1 or j >= nstep - 2:
            	constraints.append(x[dim*j:dim*(j+1)] - oripath.squeeze()[dim*j:dim*(j+1)] == 0)
            
            # Define distance constraint
            for l in range(int(dim/2)):                              # Loop through multiple cars
                for m in range(l+1, int(dim/2)):                     # Loop through other cars to define the distance constraint
                    ref_point_1 = x_ref_1[2*l : 2*(l+1)].reshape(2,1)
                    ref_point_2 = x_ref_1[2*m : 2*(m+1)].reshape(2,1)
                    ref_points = get_tayler_expansion_points(ref_point_1, ref_point_2, mini_distance)
                    A_step = 2 * ref_points.T @ I_2
                    b = -mini_distance**2 - 1/2 * (A_step @ ref_points)
                    A = -A_step

                    #print(A[0, 0:2].shape, x[dim*j:dim*(j+1)][2*l:2*(l+1)].shape, b.shape)
                    cons = A[0, 0:2].reshape(1,2)@x[dim*j:dim*(j+1)][2*l:2*(l+1)] + A[0, 2:4].reshape(1,2)@x[dim*j:dim*(j+1)][2*m:2*(m+1)] <= b
                    constraints.append(cons)
                    #break
                #break
            
            # Define priority constraint
            # ...
            #print("x_ref_1 shape: {}, x_ref_2 shape: {}".format(x_ref_1.shape, x_ref_2.shape))
            if j < nstep-1:
                coe1 = get_line(x_ref_1[0], x_ref_1[1], x_ref_2[0], x_ref_2[1])
                coe2 = get_line(x_ref_1[2], x_ref_1[3], x_ref_2[2], x_ref_2[3])
                k1 = coe1[0:2]
                b1 = coe1[2]
                k2 = coe2[0:2]
                b2 = coe2[2]
                #print("k shape: {}".format(k1.shape))
                if k1.T @ x_ref_1[2:4] + b1 > 0 and k1.T @ x_ref_2[2:4] + b1 < 0 and (k2.T @ x_ref_1[0:2] + b2)*(k2.T @ x_ref_2[0:2] + b2) < 0:
                    cons1 = -k1.T @ x[dim*j+dim/2 : dim*(j+1)] <= b1
                    cons2 = k1.T @ x[dim*(j+1)+dim/2 : dim*(j+2)] <= -b1
                    constraints.append(cons1)
                    constraints.append(cons2)
                    if k2.T @ x_ref_1[0:2] + b2 > 0:
                        cons3 = -k2.T @ x[dim*j : dim*j + dim/2] <= b2
                        cons4 = -k2.T @ x[dim*(j+1) : dim*(j+1) + dim/2] <= b2
                        constraints.append(cons3)
                        constraints.append(cons4)
                    else:
                        cons3 = k2.T @ x[dim*j : dim*j + dim/2] <= -b2
                        cons4 = k2.T @ x[dim*(j+1) : dim*(j+1) + dim/2] <= -b2
                        constraints.append(cons3)
                        constraints.append(cons4)
                
                elif k1.T @ x_ref_1[2:4] + b1 < 0 and k1.T @ x_ref_2[2:4] + b1 > 0 and (k2.T @ x_ref_1[0:2] + b2)*(k2.T @ x_ref_2[0:2] + b2) < 0:
                    cons1 = k1.T @ x[dim*j+dim/2 : dim*(j+1)] <= -b1
                    cons2 = -k1.T @ x[dim*(j+1)+dim/2 : dim*(j+2)] <= b1
                    constraints.append(cons1)
                    constraints.append(cons2)
                    if k2.T @ x_ref_1[0:2] + b2 > 0:
                        cons3 = -k2.T @ x[dim*j : dim*j + dim/2] <= b2
                        cons4 = -k2.T @ x[dim*(j+1) : dim*(j+1) + dim/2] <= b2
                        constraints.append(cons3)
                        constraints.append(cons4)
                    else:
                        cons3 = k2.T @ x[dim*j : dim*j + dim/2] <= -b2
                        cons4 = k2.T @ x[dim*(j+1) : dim*(j+1) + dim/2] <= -b2
                        constraints.append(cons3)
                        constraints.append(cons4)

        # constraints.append(-2 <= 0)
        p = Problem(objective, constraints)
        primal_result = p.solve(solver = CVXOPT)

        pathnew = x.value
        print("pathnew is of shape: {}".format(pathnew.shape))
        print("pathnew is {}".format(pathnew))

        diff = LA.norm(refpath - pathnew)
        print("diff is: ", diff)
        if diff < 0.001*nstep*dim:
            print("Converged at step {}".format(i))
            break

        refpath = pathnew

        # Show planned path
        for i in range(2):
            x = pathnew[2*i : : 4]
            y = pathnew[2*i+1 : : 4]
            plt.plot(y, x)
        plt.legend(('1', '2', '3', '4', '5', '6', '7', '8', '9'))
        plt.show()




    print("Loop finished!")


    # #Calculate the distance between two agents at every step.
    # distance = np.zeros((nstep, 1))
    # for i in range(nstep):
        # distance[i] = LA.norm(pathnew[dim*i : dim*i+2] - pathnew[dim*i+2 : dim*i+4])

    return pathnew



def main():
    MAX_ITER = 10

    multi_path = []
    # Define path 0
    path_seg_0 = np.array([[0, 0], [20, 20]])
    path_0 = get_path(path_seg_0, 1)
    # print(path_0.shape)
    # print(path_0)
    multi_path.append(path_0)
    # Define path 1
    path_seg_1 = np.array([[0, 20], [20, 0]])
    path_1 = get_path(path_seg_1, 1)
    # print(path_1.shape)
    # print(path_1)
    multi_path.append(path_1)

    # # Define path 2
    # path_seg_2 = np.array([[0, 30], [0, 130]])
    # path_2 = get_path(path_seg_2, 3)
    # # print(path_2.shape)
    # # print(path_2)
    # multi_path.append(path_2)
    # # Define path 3
    # path_seg_3 = np.array([[3.5, 15], [3.5, 115]])
    # path_3 = get_path(path_seg_3, 3)
    # # print(path_3.shape)
    # # print(path_3)
    # multi_path.append(path_3)
    # # Define path 4
    # path_seg_4 = np.array([[3.5, 0], [3.5, 100]])
    # path_4 = get_path(path_seg_4, 3)
    # # print(path_4.shape)
    # # print(path_4)
    # multi_path.append(path_4)
    # # Define path 5
    # path_seg_5 = np.array([[3.5, -20], [3.5, 80]])
    # path_5 = get_path(path_seg_5, 3)
    # # print(path_5.shape)
    # # print(path_5)
    # multi_path.append(path_5)
    # # Define path 6
    # path_seg_6 = np.array([[0, -15], [0, 0], [3.5, 0], [3.5, 85]])
    # path_6 = get_path(path_seg_6, 3)
    # # print(path_6.shape)
    # # print(path_6)
    # multi_path.append(path_6)
    # # Define path 7
    # path_seg_7 = np.array([[-3.5, -25], [-3.5, -10], [0, -10], [0, 75]])
    # path_7 = get_path(path_seg_7, 3)
    # # print(path_7.shape)
    # # print(path_7)
    # multi_path.append(path_7)
    # # Define path 8
    # path_seg_8 = np.array([[-3.5, -5], [-3.5, 95]])
    # path_8 = get_path(path_seg_8, 3)
    # # print(path_8.shape)
    # # print(path_8)
    # multi_path.append(path_8)
    # # path_seg_2 = np.array()


    multi_path = np.array(multi_path)
    # print("Converted shape is: {}".format(multi_path.shape))

    # for i in range(len(multi_path)):
    #     multi_path[i] = multi_path[i][0:35, :]
    #     print(multi_path[i].shape)




    pathnew = Plan_trajectory(MAX_ITER, multi_path, 3)
    print(type(pathnew), pathnew.shape)
    nstep = int(pathnew.shape[0] / 4)
    pathnew1 = np.empty((nstep, 2))
    pathnew2 = np.empty((nstep, 2))
    pathnew1[:, 0] = pathnew[0 : : 4]
    pathnew1[:, 1] = pathnew[1 : : 4]
    pathnew2[:, 0] = pathnew[2 : : 4]
    pathnew2[:, 1] = pathnew[3 : : 4]

    for i in range(2):
        x = pathnew[2*i : : 4]
        y = pathnew[2*i+1 : : 4]
        plt.plot(y, x)

    plt.legend(('1', '2', '3', '4', '5', '6', '7', '8', '9'))
    plt.show()
    
    two_car_distance(pathnew1, pathnew2)
    

    


if __name__ == "__main__":
    main()

#The iteration
# MAX_ITER = 10
# nstep = 50
# mini_distance = 10
# dim = 4
# start_1 = [0, 0];
# ending_1 = [200, 0]; 
# start_2 = [100, 0.5];
# ending_2 = [100, 0.5];

# start_time = time.time()
# pathnew, distance = Plan_trajectory(MAX_ITER, start_1, ending_1, start_2, ending_2, mini_distance, nstep, dim)
# print("Calculating time: {}".format(time.time() - start_time))






# plt.axis([-10, 210, -100, 100])
# plt.plot(pathnew[2*0 :  : dim], pathnew[2*0+1 :  : dim], 'bo')
# plt.plot(pathnew[2*1 :  : dim], pathnew[2*1+1 :  : dim], 'rs')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title('Cozmo Trajectory Generation!')
# plt.show()

















