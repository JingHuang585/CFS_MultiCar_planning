import numpy as np
import time
from CFS_problem import *
import cvxopt
from cvxopt import matrix, solvers
from numpy import linalg as LA
import matplotlib.pyplot as plt
from cvxpy import *
import matplotlib.pyplot as plt
#print(Q3)

# print("A is:", A)
# print("A is of shape:", A.shape)

# x_ref = refpath[dim*0:dim*(0+1)]
# A_step = 2 * (x_ref.T @ I_2)

# print(x_ref)
# print(A_step)

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



def get_tyler_expansion_points(high_pri_pt, low_pri_pt, mini_distance):
    new_ref = np.zeros((2, 1))
    dx = high_pri_pt[0] - low_pri_pt[0]
    dy = high_pri_pt[1] - low_pri_pt[1]
    new_ref[0] = high_pri_pt[0] - mini_distance * dx / np.sqrt(dx**2 + dy**2)
    new_ref[1] = high_pri_pt[1] - mini_distance * dy / np.sqrt(dx**2 + dy**2)
    
    return new_ref

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

    

    # for i in range(9):
    #     x = oripath[2*i : : 18]
    #     y = oripath[2*i+1 : : 18]

    #     plt.plot(x, y)
    # plt.show()
    for i in range(MAX_ITER):
        print(i)

        x = Variable(n)
        # print((Qref.T@oripath).shape)
        objective = Minimize(0.5*quad_form(x,Qe) + (np.matmul(-Qref,oripath)).T*x)
        constraints = []
        for j in range(nstep):
            x_ref = refpath[dim*j:dim*(j+1)]
            if j <= 1 or j >= nstep - 2:
            	constraints.append(x[dim*j:dim*(j+1)] - oripath.squeeze()[dim*j:dim*(j+1)] == 0)

            # for m in range(int(dim/2)):
            #     for n in range(m+1, int(dim/2)):
            #         ref_point_1 = x_ref[2*m : 2*(m+1)].reshape(2,1)
            #         point_2 = x_ref[2*n : 2*(n+1)]
            #         ref_point_2 = get_tyler_expansion_points(ref_point_1, point_2, mini_distance)
            #         ref_points = np.vstack((ref_point_1, ref_point_2))
            #         A = 2 * ref_points.T @ I_2
            #         b = -mini_distance**2 - 1/2 * A @ ref_points

            #         print(A[0, 0:2].shape, x[dim*j:dim*(j+1)][2*m:2*(m+1)].shape, b.shape)
            #         cons = A[0, 0:2].reshape(1,2)@x[dim*j:dim*(j+1)][2*m:2*(m+1)] + A[0, 2:4].reshape(1,2)@x[dim*j:dim*(j+1)][2*n:2*(n+1)] <= b

            #         constraints.append(cons)
            #         break
            #     break





        # constraints.append(-2 <= 0)
        p = Problem(objective, constraints)
        primal_result = p.solve(solver = CVXOPT)

        pathnew = x.value



        print("pathnew is of shape: {}".format(pathnew.shape))
        print("pathnew is {}".format(pathnew))



        #Check if the distance constraint is satisfied
        # test = A @ pathnew;
        # print("A is of shape: ", A.shape)
        # print("pathnew is of shape: ", pathnew.shape)
        # print("test is of shape: ", test.shape)
        # a = (test<=b)
        # for e in a:
        # 	print (e)


        diff = LA.norm(refpath - pathnew)
        print("diff is: ", diff)
        if diff < 0.001*nstep*dim:
            print("Converged at step {}".format(i))
            break


        refpath = pathnew


    print("Loop finished!")


    # #Calculate the distance between two agents at every step.
    # distance = np.zeros((nstep, 1))
    # for i in range(nstep):
        # distance[i] = LA.norm(pathnew[dim*i : dim*i+2] - pathnew[dim*i+2 : dim*i+4])

    return pathnew



MAX_ITER = 10

multi_path = []
# Define path 0
path_seg_0 = np.array([[0, 0], [0, 100]])
path_0 = get_path(path_seg_0, 3)
# print(path_0.shape)
# print(path_0)
multi_path.append(path_0)
# Define path 1
path_seg_1 = np.array([[-3.5, 15], [-3.5, 30], [0, 30], [0, 115]])
path_1 = get_path(path_seg_1, 3)
# print(path_1.shape)
# print(path_1)
multi_path.append(path_1)
# Define path 2
path_seg_2 = np.array([[0, 30], [0, 130]])
path_2 = get_path(path_seg_2, 3)
# print(path_2.shape)
# print(path_2)
multi_path.append(path_2)
# Define path 3
path_seg_3 = np.array([[3.5, 15], [3.5, 115]])
path_3 = get_path(path_seg_3, 3)
# print(path_3.shape)
# print(path_3)
multi_path.append(path_3)
# Define path 4
path_seg_4 = np.array([[3.5, 0], [3.5, 100]])
path_4 = get_path(path_seg_4, 3)
# print(path_4.shape)
# print(path_4)
multi_path.append(path_4)
# Define path 5
path_seg_5 = np.array([[3.5, -20], [3.5, 80]])
path_5 = get_path(path_seg_5, 3)
# print(path_5.shape)
# print(path_5)
multi_path.append(path_5)
# Define path 6
path_seg_6 = np.array([[0, -15], [0, 0], [3.5, 0], [3.5, 85]])
path_6 = get_path(path_seg_6, 3)
# print(path_6.shape)
# print(path_6)
multi_path.append(path_6)
# Define path 7
path_seg_7 = np.array([[-3.5, -25], [-3.5, -10], [0, -10], [0, 75]])
path_7 = get_path(path_seg_7, 3)
# print(path_7.shape)
# print(path_7)
multi_path.append(path_7)
# Define path 8
path_seg_8 = np.array([[-3.5, -5], [-3.5, 95]])
path_8 = get_path(path_seg_8, 3)
# print(path_8.shape)
# print(path_8)
multi_path.append(path_8)
# path_seg_2 = np.array()


multi_path = np.array(multi_path)
# print("Converted shape is: {}".format(multi_path.shape))

for i in range(len(multi_path)):
    multi_path[i] = multi_path[i][0:35, :]
    print(multi_path[i].shape)




pathnew = Plan_trajectory(MAX_ITER, multi_path, 4)


print(type(pathnew), pathnew.shape)

for i in range(9):
	x = pathnew[2*i : : 18]
	y = pathnew[2*i+1 : : 18]
	plt.plot(y, x)

plt.legend(('1', '2', '3', '4', '5', '6', '7', '8', '9'))
plt.show()



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

















