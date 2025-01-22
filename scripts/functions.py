import numpy as np


# the theory of cooling oceanic litosphere in geo3 script1 P. 124

# i need a reference depth at specific age for the sqrt(t) relation


# the not ideal way to downsample data. Lot of processing.
def compress_grid(grid, block_size):
    compressed_grid = np.zeros((grid.shape[0] // block_size, grid.shape[1] // block_size))
    compressed_lat = np.linspace(-90, 90, grid.shape[0] // block_size)
    compressed_lon = np.linspace(-180, 180, grid.shape[1] // block_size)
    print(f'Shape of compressed grid: {np.shape(compressed_grid)}')
    for i in range(0, grid.shape[0] - 1, block_size):  # calculate blocks
        for j in range(0, grid.shape[1] - 1, block_size):
            # fill compressed_grid with values
            compressed_grid[i // block_size, j // block_size] = grid[i:i + block_size, j:j + block_size].mean()

    return compressed_lon, compressed_lat, compressed_grid


# Does not work properly right now
def downsample_grid(grid, block_size):
    y, x = grid.shape
    new_y = (y // block_size) * block_size
    new_x = (x // block_size) * block_size
    trimmed_grid = grid[:new_y, :new_x] # trimmed grid to match the divisible size
    downsampled_lat = np.linspace(-90, 90, new_y // block_size)
    downsampled_lon = np.linspace(-180, 180, new_x // block_size)

    # Reshape into 4 dimensions-> each new pixel contains a 2D block -> 4D
    reshaped = trimmed_grid.reshape(new_x // block_size, block_size, new_y // block_size, block_size)
    downsampled_grid = reshaped.mean(axis=(1, 3))  # Average across block rows and block columns
    print(np.shape(downsampled_grid))

    return downsampled_lat, downsampled_lon, downsampled_grid


def age_to_depth(grid):
    m = 0.3
    depth_grid = m * np.sqrt(1000000*grid)+2000
    return depth_grid
