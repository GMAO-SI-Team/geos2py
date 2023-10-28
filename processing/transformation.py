import numpy as np

def flip_hemispheres(image: np.ndarray, frac=None) -> np.ndarray:
    image_dump = image
    shape = image.shape
    split = int(shape[0] * frac) if frac else int(shape[0] / 2.0)

    if len(shape) == 1:
        image_dump[:split] = image[split:shape[0] - 1]
        image_dump[split:shape[0] - 1] = image[:split]
        image = image_dump

    if len(shape) == 2:
        image_dump[:,:split] = image[:,split:shape[0] - 1]
        image_dump[:,split:shape[0] - 1] = image[:,:split]
        image = image_dump

    if len(shape) == 3:
        for k in range(shape[0]):
            image_dump[k,:,:split] = image[k,:,split:shape[0] - 1]
            image_dump[k,:,split:shape[0] - 1] = image[k,:,:split]
        image = image_dump

    return image