import time
from functools import wraps
import skimage.filters as filters

# ------------ Python File to experiment and improve performance -----------


def measure(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = f(*args, **kwargs)
        end_time = time.time()
        print("This took: " + str(end_time-start_time) + "s")
        return result
    return wrapper


# vol = utils.load_volume("Lunge10um_mini/")
#
# vol = vol[::2, ::2, ::2]

@measure
def gaussian(vol, sigma=1):
    filters.gaussian(vol, sigma=sigma)


#gaussian(vol, 1)
#gaussian(vol, 100)

from multiprocessing import Process

if __name__ == '__main__':
    start_time = time.time()
    processes = []
    for i in range(8):
        processes.append(
            Process(target=gaussian, args=(vol, 20))
        )
        processes[i].start()
        processes[i].join()

    for process in processes:
        process.join()

    end_time = time.time()
    print(str(end_time-start_time) + "s")
