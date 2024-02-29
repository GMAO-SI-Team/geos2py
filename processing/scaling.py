import numpy as np

def bytscl(image: np.ndarray, low=None, high=None) -> np.ndarray:
    """
    Scale the pixel values of an image to the range [0, 1].

    Parameters:
    - image (np.ndarray): Input image.
    - low (float, optional): Lower bound of input data. If None, it's the minimum of the input data. Default is None.
    - high (float, optional): Upper bound of input data. If None, it's the maximum of the input data. Default is None.

    Returns:
    - np.ndarray: Image with pixel values scaled to the range [0, 1].
    """
    if low is None:
        low = np.min(image)
    if high is None:
        high = np.max(image)

    scaled_image = ((image - low) / (high - low)) * 255
    scaled_image = np.clip(scaled_image, 0, 255)
    return scaled_image / 255