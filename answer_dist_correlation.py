import series_matching as sm
import os
import scipy.io
import numpy as np
from matplotlib import pyplot as plt
import scipy.io

# def draw_line()

def calculate_answer_matrix(info):
    dimension = len(info)
    answer_dis=np.zeros((dimension,dimension),dtype=np.float)
    answer_dis_penalty=np.zeros((dimension,dimension),dtype=np.float)
    for i in range(len(info)):
        for j in range(len(info)):
            for k in range(4):
                if info[i]['score'][k]==1 and info[j]['score'][k] == 1:
                    answer_dis[i][j] += 1
                    answer_dis_penalty[i][j] += 1
                elif info[i]['score'][k]+info[j]['score'][k] == 0:
                    answer_dis_penalty[i][j] += 0.5
    answer_dis = answer_dis / np.max(answer_dis)
    answer_dis_penalty = answer_dis_penalty / np.max(answer_dis_penalty)
    print answer_dis, answer_dis_penalty
    return answer_dis, answer_dis_penalty

def draw_fig(x,y,labels):

    plt.subplots_adjust(bottom=0.1)
    plt.scatter(
        x, y, marker='o',
        cmap=plt.get_cmap('Spectral'))

    for label, x, y in zip(labels, x, y):
        plt.annotate(
            label,
            xy=(x, y), xytext=(-20, 20),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    plt.show()


if __name__ == "__main__":
    # fig = PLT.figure()
    # ax1 = fig.add_subplot(111)
    # ax1.fill_between(x, 0, y_stack[0, :], facecolor="#CC6666", alpha=.7)
    # ax1.fill_between(x, y_stack[0, :], y_stack[1, :], facecolor="#1DACD6", alpha=.7)
    # ax1.fill_between(x, y_stack[1, :], y_stack[2, :], facecolor="#6E5160")

    # 405-435
    simon=[
        {'sub':'an1','score':[0,0,0,0,1]},
        {'sub':'at1','score':[0,0,1,1,0]},
        {'sub':'azz1','score':[1,0,0,0,1]},
        {'sub':'cd1','score':[1,0,0,0,1]},
        {'sub':'klt1','score':[0,1,1,1,0]},
        {'sub':'lr1','score':[0,0,0,0,0]},
        {'sub':'lv1','score':[1,0,1,0,0]},
        {'sub':'msb1','score':[1,0,0,0,1]},
        {'sub':'pc1','score':[0,1,1,0,1]},
        {'sub':'pp1','score':[0,0,0,0,0]},
        {'sub':'sw1','score':[0,1,1,1,0]}
    ]
    answer_dis, answer_dis_penalty = calculate_answer_matrix(simon)
    TWED_matrix = scipy.io.loadmat('result/30s/s_location/lam10000_nu10000/dis_14.mat')['matrix'][0]
    print TWED_matrix
    x = []
    y=[]
    labels = []
    for i in range(len(TWED_matrix)):
        for j in range(len(TWED_matrix[0])):
            if i > j:
                x.append(answer_dis[i][j])
                y.append(TWED_matrix[i][j])
                labels.append(simon[i]['sub'] + '_' + simon[j]['sub'])
    # x = []
    # y=[]
    # labels=[]
    # for i in range(len(TWED_matrix)):
    #     for j in range(len(TWED_matrix[0])):
    #         if i>=j:
    #             x.append(answer_dis_penalty[i][j])
    #             y.append(TWED_matrix[i][j])
    #             labels.append(simon[i]['sub']+'_'+simon[j]['sub'])
    draw_fig(x,y,labels)
    # plt.scatter(x_p, y_p)
    # plt.show()