import os

class Path(object):
    def __init__(self, data_dir: str, f_date: str):
        self.f_year = f_date
        self.f_month = f_date
        self.f_day = f_date.split('_')[0][-2:]
        self.f_hour = f_date.split('_')[1][:2]
        self.s_year = 0
        self.s_month = 0
        self.s_day = 0
        self.s_hour = 0
        self.s_minute = 0
        self.data_dir_head = f'{data_dir}/Y{self.f_year}/M{self.f_month}/D{self.f_day}/H{self.f_hour}'
        self.data_dir_tail = f'{f_date}+{self.s_year}{self.s_month}{self.s_day}_{self.s_hour}{self.s_minute}.V01.nc4'

class IRPath(Path):
    def __init__(self, data_dir, f_date):
        super().__init__(self)
        self.full_data_dir = f'{self.data_dir_head}/GEOS.fp.fcst.inst_30mn_met_c0720sfc.{self.data_dir_tail}'

class RadarPath(Path):
    def __init__(self):
        super().__init__(self)
        self.full_data_dir = f'{self.data_dir_head}/GEOS.fp.fcst.inst_30mn_met_c0720sfc.{self.data_dir_tail}'

class WxtypePath(Path):
    def __init__(self):
        super().__init__(self)
        self.full_data_dir = (
            f'{self.data_dir_head}/GEOS.fp.fcst.tavg1_2d_slv_Nx.{self.data_dir_tail}',
            f'{self.data_dir_head}/GEOS.fp.fcst.tavg1_2d_flx_Nx.{self.data_dir_tail}',
            f'{self.data_dir_head}/GEOS.fp.fcst.inst3_3d_asm_Np.{self.data_dir_tail}'
        )

class AerosolPath(Path):
    def __init__(self):
        super().__init__(self)
        self.full_data_dir = f'{self.data_dir_head}/GEOS.fp.fcst.inst1_2d_hwl_Nx.{self.data_dir_tail}'

class PrecrainPath(Path):
    def __init__(self):
        super().__init__(self)
        self.full_data_dir = f'{self.data_dir_head}/GEOS.fp.fcst.tavg1_2d_flx_Nx.{self.data_dir_tail}'

class PrecsnowPath(Path):
    def __init__(self):
        super().__init__(self)
        self.full_data_dir = f'{self.data_dir_head}/GEOS.fp.fcst.tavg1_2d_flx_Nx.{self.data_dir_tail}'

class SLPPath(Path):
    def __init__(self):
        super().__init__(self)
        self.full_data_dir = f'{self.data_dir_head}/GEOS.fp.fcst.tavg1_2d_slv_Nx.{self.data_dir_tail}'

class T2MPath(Path):
    def __init__(self):
        super().__init__(self)
        self.full_data_dir = f'{self.data_dir_head}/GEOS.fp.fcst.tavg1_2d_slv_Nx.{self.data_dir_tail}'

class CAPEPath(Path):
    def __init__(self):
        super().__init__(self)
        self.full_data_dir = f'{self.data_dir_head}/GEOS.fp.fcst.inst3_2d_met_Nx.{self.data_dir_tail}'

class TPWPath(Path):
    def __init__(self):
        super().__init__(self)
        self.full_data_dir = '/discover/nobackup/projects/gmao/gmao_ops/pub/'+pTag+'/forecast/Y'+fYear+'/M'+fMon+'/D'+fDay+'/H'+fHour+'/GEOS.fp.fcst.tavg1_2d_slv_Nx.'+FDATE+'+'+sYear0+sMon0+sDay0+'_'+sHour0+'30.V01.nc4'

class VortPath(Path):
    def __init__(self):
        super().__init__(self)
        self.full_data_dir = f'{self.data_dir_head}/GEOS.fp.fcst.tavg1_2d_slv_Nx.{self.data_dir_tail}'

class WindPath(Path):
    def __init__(self):
        super().__init__(self)
        f'{self.data_dir_head}/GEOS.fp.fcst.tavg1_2d_flx_Nx.{self.data_dir_tail}'