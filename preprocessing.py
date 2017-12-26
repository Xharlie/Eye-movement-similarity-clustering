import numpy as np

def load(file):
    return np.loadtxt(file)

def down_sample(data, sample_rate):
    interval = int(sample_rate)
    return data[::interval]

if __name__ == "__main__":
    x = load("data/c_location/yo10.txt")
    print x.shape
    x = down_sample(x,0.05)
    print x[0,0]