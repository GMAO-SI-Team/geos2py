from sunpy.image import resample as rs

def resample(data, shape, center=True):
    return rs.resample(data, shape, center=center)