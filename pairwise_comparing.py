from matplotlib import pyplot as plt
from matplotlib import cm as cm
import matplotlib.colors
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
import parula
import numpy as np
import series_matching
import matplotlib.patches as patches
import os.path
import scipy.io
import preprocessing

def highResPoints(x,y,factor=10, RESFACT= 10):
    '''
    Take points listed in two vectors and return them at a higher
    resultion. Create at least factor*len(x) new points that include the
    original points and those spaced in between.

    Returns new x and y arrays as a tuple (x,y).
    '''

    # r is the distance spanned between pairs of pointsd
    if RESFACT == 1:
        return x,y
    r = [0]
    for i in range(1,len(x)):
        dx = x[i]-x[i-1]
        dy = y[i]-y[i-1]
        r.append(np.sqrt(dx*dx+dy*dy))
    r = np.array(r)

    # rtot is a cumulative sum of r, it's used to save time
    rtot = []
    for i in range(len(r)):
        rtot.append(r[0:i].sum())
    rtot.append(r.sum())
    dr = rtot[-1]/(len(x)*RESFACT-1)
    xmod=[x[0]]
    ymod=[y[0]]
    rPos = 0 # current point on walk along data
    rcount = 1
    while rPos < r.sum():
        x1,x2 = x[rcount-1],x[rcount]
        y1,y2 = y[rcount-1],y[rcount]
        dpos = rPos-rtot[rcount]
        theta = np.arctan2((x2-x1),(y2-y1))
        rx = np.sin(theta)*dpos+x1
        ry = np.cos(theta)*dpos+y1
        xmod.append(rx)
        ymod.append(ry)
        rPos+=dr
        while rcount+1 < len(rtot) and rPos > rtot[rcount+1]:
            rPos = rtot[rcount+1]
            rcount+=1
            if rcount>rtot[-1]:
                break

    return xmod,ymod

def correlation_matrix(df, x, y, stride,window_length, o1=0, o2=1, threeD=False, t = 0, time_range=None):
    RESFACT = 1
    df=df[t]
    x_window = x[t * stride:t * stride + window_length, :] / 1000.0
    y_window = y[t * stride:t * stride + window_length, :] / 1000.0
    fig = plt.figure(figsize=(15,5))
    if threeD:
        tra = fig.add_subplot(121, projection='3d')
        tra.set_xlabel('X(pixel)')
        tra.set_ylabel('Y(pixel)')
        tra.set_zlabel('T(s)')
    else:
        tra = fig.add_subplot(121)
        tra.set_xlabel('X(pixel)')
        tra.set_ylabel('Y(pixel)')
    xHiRes, yHiRes = highResPoints(x_window[:,o1], y_window[:,o1], RESFACT = RESFACT)
    npointsHiRes = len(xHiRes)
    for i in range(npointsHiRes - 1):
        if threeD:
            tra.plot(xHiRes[i:i + 2], yHiRes[i:i + 2], time_range[0] + (time_range[1]-time_range[0]) * i / (npointsHiRes - 1),
                     alpha=float(i) / (npointsHiRes - 1),
                     color='r', linewidth=1 + i * 4.0 / (npointsHiRes - 1))
        else:
            tra.plot(xHiRes[i:i + 2], yHiRes[i:i + 2],
                     alpha=float(i) / (npointsHiRes - 1),
                     color='r', linewidth= 1 + i * 4.0 / (npointsHiRes - 1))

    xHiRes, yHiRes = highResPoints(x_window[:,o2], y_window[:,o2], RESFACT = RESFACT)
    npointsHiRes = len(xHiRes)
    for i in range(npointsHiRes - 1):
        if threeD:
            tra.plot(xHiRes[i:i + 2], yHiRes[i:i + 2], time_range[0] + (time_range[1]-time_range[0]) * i / (npointsHiRes - 1),
                     alpha=float(i) / (npointsHiRes - 1),
                     color='b', linewidth=1 + i * 4.0 / (npointsHiRes - 1))
        else:
            tra.plot(xHiRes[i:i + 2], yHiRes[i:i + 2],
                     alpha=float(i) / (npointsHiRes - 1),
                     color='b', linewidth= 1 + i * 4.0 / (npointsHiRes - 1))

    ax1 = fig.add_subplot(122)

    parula_map = LinearSegmentedColormap.from_list('parula', parula.ger_parula())
    # For use of "viscm view"
    test_cm = parula_map

    # viscm(parula_map)

    # cmap = cm.get_cmap('jet', 30)
    norm = matplotlib.colors.Normalize(vmin=0.3, vmax=1)

    cax = ax1.imshow(df, interpolation="nearest", cmap=parula_map, norm=norm)
    ax1.grid(True)
    plt.title('TWED normalized')
    labels=range(df.shape[0])
    ax1.set_xticks(np.arange(0, df.shape[0], df.shape[0] * 1.0 / len(labels)))
    ax1.set_yticks(np.arange(0, df.shape[1], df.shape[1] * 1.0 / len(labels)))
    ax1.set_xticklabels(labels,fontsize=6)
    ax1.set_yticklabels(labels,fontsize=6)
    # Add colorbar, make sure to specify tick locations to match desired ticklabels

    fig.colorbar(cax, ticks=[.25, .45, .65, .85, .95, 1], spacing='proportional')
    ax1.add_patch(
        patches.Circle(
            (o1 * 1.0 , o2*1.0),
            df.shape[0] * 1.0 / len(labels) / 4,
            color=(1, 1 ,1)
        )
    )
    ax1.add_patch(
        patches.Circle(
            (o2 * 1.0, o1 * 1.0),
            df.shape[0] * 1.0 / len(labels) / 4,
            color=(1, 1, 1)
        )
    )
    plt.show()



if __name__ == "__main__":
    lam = 5000
    nu = 100000
    stride = 32 * 5
    window_length = 32 * 5
    frame_per_second = 32
    time_length = 30
    segment=2
    t = 3
    o1 = 1
    o2 = 8
    threeD=True

    x_file = 'data/c_location/xo{}.txt'.format(segment)
    y_file = 'data/c_location/yo{}.txt'.format(segment)
    # the x,y value is 489871.7
    dest = 'result/5s/c_location/lam{}_nu{}/dis_{}.mat'.format(lam,nu,segment)
    if os.path.exists(dest):
        series_similarity = scipy.io.loadmat(dest)['matrix']
        # print series_similarity
        x = preprocessing.load(x_file)
        y = preprocessing.load(y_file)
        sample_rate = x.shape[0] / (frame_per_second * time_length)
        if (x.shape[0] != y.shape[0] or x.shape[1] != y.shape[1]):
            print("shape mismatch", x.shape, y.shape)
        x = preprocessing.down_sample(x, sample_rate)
        y = preprocessing.down_sample(y, sample_rate)
    else:
        series_similarity,x,y = series_matching.sliding_window(x_file=x_file, y_file=y_file, lam=lam, nu=nu,
            stride=stride, window_length=window_length, frame_per_second=frame_per_second, time_length=time_length)
    time_range = [time_length * (segment-1) + (t - 1) * window_length / frame_per_second,
                 time_length * (segment - 1) + t * window_length / frame_per_second]
    correlation_matrix(series_similarity,x,y,stride,window_length,o1=o1, o2=o2, threeD=threeD, t = t,time_range=time_range)
