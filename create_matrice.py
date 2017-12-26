
import series_matching as sm
import os
import scipy.io

def check_create_dir(dir):
    try:
        os.stat(dir)
    except:
        os.mkdir(dir)
    return dir


if __name__ == "__main__":
    stride = 300
    window_duraing= 30
    frame_per_second = 32
    window_length = frame_per_second * window_duraing
    time_length = 30

    for data_source in ['c_location','s_location']:
        for lam in [0, 1000, 10000, 20000, 50000, 100000]:
            for nu in [0, 1000, 10000, 20000, 50000,100000]:
                result_dir="result/{}s/{}/lam{}_nu{}/".format(window_duraing,data_source,lam,nu)
                check_create_dir(result_dir)
                for file_num in range(1,21):
                    x_file = 'data/{}/xo{}.txt'.format(data_source,file_num)
                    y_file = 'data/{}/yo{}.txt'.format(data_source,file_num)
                    matrix,_,_ = sm.sliding_window(x_file=x_file, y_file=y_file, lam=lam, nu=nu,
                                   stride=stride, window_length=window_length,
                                      frame_per_second=frame_per_second, time_length=time_length)
                    scipy.io.savemat(result_dir+'dis_{}.mat'.format(file_num), mdict={'matrix': matrix})
                    print("saved:"+result_dir+'dis_{}.mat'.format(file_num))





