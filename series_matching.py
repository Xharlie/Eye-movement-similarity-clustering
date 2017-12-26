import numpy as np
import preprocessing
import math
import TWED


def sliding_window(x_file=None,y_file=None,lam=None,nu=None,
                   stride=None,window_length=None,frame_per_second=None,time_length=30):

    x = preprocessing.load(x_file)
    y = preprocessing.load(y_file)
    sample_rate = x.shape[0] / (frame_per_second * time_length)
    if ( x.shape[0] !=  y.shape[0] or x.shape[1] != y.shape[1]):
        print("shape mismatch",x.shape,y.shape)
    x = preprocessing.down_sample(x, sample_rate)
    y = preprocessing.down_sample(y, sample_rate)
    periods = int(math.floor((x.shape[0] - window_length) / stride) + 1)

    series_distance = np.zeros((periods,x.shape[1],x.shape[1]))
    for t in range(periods):
        print "periods:{}".format(t)
        x_window = x[t * stride:t * stride + window_length,:]
        y_window = y[t * stride:t * stride + window_length,:]
        for i in range(x.shape[1]):
            for j in range(x.shape[1]):
                p1_data = np.stack((x_window[:,i], y_window[:,i]))
                p2_data = np.stack((x_window[:,j], y_window[:,j]))
                series_distance[t,i,j] = TWED.TWED(p1_data,p2_data,lam,nu)

    series_similarity = 1 - series_distance / np.max(series_distance)

    # print series_distance[0,:,:]
    # correlation_matrix(series_similarity[10,:,:])
    return series_similarity, x, y

if __name__ == "__main__":
    x_file = 'data/c_location/xo1.txt'
    y_file = 'data/c_location/yo1.txt'
    # the x,y value is 489871.7
    lam = 5000
    nu = 5000
    stride = 32
    window_length = 32 * 5
    frame_per_second = 32
    time_length = 30
    sliding_window(x_file=x_file, y_file=y_file, lam=lam, nu=nu,
        stride=stride, window_length=window_length, frame_per_second=frame_per_second, time_length=time_length)

