import numpy as np
from numpy import linalg as LA


def define_path():
    '''
    Define initial reference paths for multiple cars.
    Return:
        multi_path: Shape: num_cars x nsteps x 2. Paths for number of cars.
    '''
    multi_path = []
    # Define path 0
    path_seg_0 = np.array([[10, 0], [30, 0]])
    path_0 = get_path(path_seg_0, 1)
    multi_path.append(path_0)
    # Define path 1
    path_seg_1 = np.array([[0, 0.01], [40, 0.01]])
    path_1 = get_path(path_seg_1, 2)
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
    return multi_path


def get_interpolate(loc_1, loc_2, resolution = 3):
    '''
    Inputs:
        loc_1, loc_2: Fisrt and last point location. Type: List. Dimension: 1x2
    '''
    # print(loc_1, loc_2)
    dist = LA.norm(loc_1 - loc_2)
    points_num = int(dist // resolution) + 1 # From loc_1 point to the last point between loc_1 and loc_2. Not including loc_2.
    points_array = np.zeros((points_num, 2))
    for i in range(points_num):
        points_array[i][0] = (points_num-i)/points_num * loc_1[0] + i/points_num * loc_2[0]
        points_array[i][1] = (points_num-i)/points_num * loc_1[1] + i/points_num * loc_2[1]

    return points_array

def get_path(path_segments_array, resolution):
    '''
    Inputs:
        path_segments_array: numpy array that contains path points that segement the path. 
        resolution: Interpolation resolution
    Outputs:
        path: numpy array that demonstrate the path. Dimension: n x 2
    '''
    num_seg = path_segments_array.shape[0]
    path = np.zeros((1, 2))
    for i in range(num_seg - 1):
        subpath = get_interpolate(path_segments_array[i], path_segments_array[i+1], resolution)
        if i == 0:
            path = subpath
        else:
            path = np.concatenate((path, subpath), axis = 0)
		
    last_point = np.array([path_segments_array[num_seg-1]])
    path = np.concatenate((path, last_point), axis = 0)
    return path


def Setup_problem(multi_path):
    # nstep = 100;
    # dim = 4;

    # start_1 = [0, 0];
    # ending_1 = [200, 0]; 


    # start_2 = [100, 0.5];
    # ending_2 = [100, 0.5];
    # margin = 0.25; 
    I_2 = np.array([[1, 0, -1, 0],
        [0, 1, 0, -1],
        [-1, 0, 1, 0],
        [0, -1, 0, 1]])

    dim = multi_path.shape[0] * 2
    nstep = multi_path[0].shape[0]


    #print(I_2)

    refpath = []
    for i in range(nstep):
        for j in range(multi_path.shape[0]):
            subpath = multi_path[j]
            refpath.append(subpath[i, 0])
            refpath.append(subpath[i, 1])

    refpath = np.array([refpath]).T
    oripath = refpath

    # for i in refpath:
    # 	print(i, '\n')

    # print("Path 1 is: ", path_1)
    # print("refpath is of shape: ", refpath.shape)
    # print("refpath is: ", refpath)

    # nstep = 4
    # dim = 2


    #Define the cost matrix
    dim = dim
    Q1 = np.eye(nstep*dim)
    Vdiff = -np.eye(nstep*dim) + np.diag(np.ones((1, (nstep-1)*dim)).squeeze(), dim)
    Q2 = Vdiff[0:(nstep-1)*dim, :].T @ Vdiff[0:(nstep-1)*dim, :]
    # print("Vdiff is:", Vdiff)
    # print("Vdiff is of shape:", Vdiff.shape)
    Adiff = -Vdiff - np.diag(np.ones((1, (nstep-1)*dim)).squeeze(), dim) + np.diag(np.ones((1, (nstep-2)*dim)).squeeze(), dim*2)
    # print("Adiff is:", Adiff)
    # print("Adiff is of shape:", Adiff.shape)

    Q3 = Adiff[0:(nstep-2)*dim, :].T @ Adiff[0:(nstep-2)*dim, :]
    # print("Q3 is:", Q3[0,:])
    # print("Q3 is of shape:", Q3.shape)


    #Define the weights
    cref = np.array([0.1, 0, 0])
    cabs = np.array([0, 0, 100])

    Qref = Q1*cref[0] + Q2*cref[1] + Q3*cref[2]
    Qabs = Q1*cabs[0] + Q2*cabs[1] + Q3*cabs[2]

    # print(Qref)
    # print(Qref.shape)

    #Fix the start and end points of all nine agents. Also, fix the path of agent 2 to 
    #be always at the initial position.
    # Aeq = np.zeros((2*dim + (nstep-2)*(int(dim/2)), nstep*dim))
    # beq = np.zeros((2*dim + (nstep-2)*(int(dim/2)), 1))

    # # Aeq = np.zeros((2*dim), nstep*dim))
    # # beq = np.zeros((2*dim), 1))
    
    # Aeq[0:dim, 0:dim] = np.eye(dim)
    # Aeq[dim:2*dim, (nstep-1)*dim:nstep*dim] = np.eye(dim)
    # beq[0:2*dim] = np.concatenate((oripath[0:4], oripath[len(oripath)-4:]), axis = 0)
    # #print(Aeq.shape)

    # for i in range(1, nstep-1):
    #     Aeq[2*dim + int(dim/2)*(i-1), dim*i + 2] = 1
    #     Aeq[2*dim + int(dim/2)*(i-1) + 1, dim*i + 3] = 1
    #     beq[2*dim + int(dim/2)*(i-1) : 2*dim + int(dim/2)*(i-1) + 2] = oripath[dim*i+2 : dim*i+4]



    return Qref, Qabs, nstep, dim, oripath, I_2




