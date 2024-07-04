import numpy as np
import napari
import sys

# Function to read uint8 binary raw file
def read_uint8_binary(file_path, shape):
    with open(file_path, 'rb') as f:
        data = np.fromfile(f, dtype=np.uint64)
        return data.reshape(shape)

if __name__ == "__main__":
    file_path = sys.argv[1]
    shape = tuple(map(int, sys.argv[2].strip('()').split(',')))

    # Read data from file
    data = read_uint8_binary(file_path, shape)

    # Display data in napari
    viewer = napari.Viewer()
    viewer.add_image(data, name="Raw Image")
    napari.run()