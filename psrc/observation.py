


class Observation:
    """
    An observation of a field value at a certain time.
    """

    def __init__(self,tm,lat,lon,obs,var,ngp):
        """
        Constructs an observation packet from the info dictionary.
        """
        self.tm = tm
        self.lon = lon
        self.lat = lat
        self.obs_val = obs
        self.obs_var = var
        self.ngp = ngp


    def get_value(self):
        """
        Return the observation value.
        """
        return self.obs_val


    def get_measurement_variance(self):
        """
        Return the variance of the measurement of the given field.
        """
        return self.obs_var


    def get_position(self):
        """
        Longitude and lattitude of the originating station (shortcut).
        """
        return (self.lat, self.lon)


    def get_nearest_grid_point(self):
        """
        Return the indices that identify the nearest grid point.
        """
        return self.ngp

    def get_time(self):
        """
        Return the GMT time of the observations.
        """
        return self.tm
