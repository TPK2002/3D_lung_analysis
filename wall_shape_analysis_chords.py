import sys
import numpy as np
import matplotlib.pyplot as plt
import napari

def do_analysis(path, shape):
    d_min = np.fromfile(path + '_min_inverse.raw', dtype=np.uint16)
    d_max = np.fromfile(path + '_max_inverse.raw', dtype=np.uint16)
    divided = np.divide(d_max, d_min, where=d_min != 0)
    divided = np.reshape(divided, shape)

    # Edges are not reliable, so we remove them
    divided = divided[5:-5, 5:-5, 5:-5]

    divided[divided <= 50] = 0
    thresholded = np.copy(divided)
    thresholded[thresholded > 0] = 1

    print('Mean:', np.mean(divided))
    print('Median:', np.median(divided))
    print('Std:', np.std(divided))

    print("Number of edge pixels: ", np.sum(thresholded))

    plt.hist(divided[divided > 0].flatten())
    plt.show()

    viewer = napari.Viewer()
    viewer.add_image(divided)
    viewer.add_image(thresholded)
    napari.run()


if __name__ == '__main__':
    if (len(sys.argv) != 5):
        print("Usage: python wall_shape_analysis_rays.py <path> <nx> <ny> <nz>")
        sys.exit(1)
    path = sys.argv[1]
    shape = (int(sys.argv[4]), int(sys.argv[3]), int(sys.argv[2]))
    do_analysis(path, shape)
