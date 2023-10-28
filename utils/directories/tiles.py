def find_tile_file(nx, ny):
    n_lons = nx * 4
    n_lats = nx * 2 + 1
    if nx in [768, 768 * 2, 768 * 4]:
        nx_str = str(nx).zfill(4)
        n_lons_str = str(n_lons).zfill(5) if nx == 3072 else str(n_lons).zfill(4)
        n_lats_str = str(n_lats).zfill(5) if nx == 3072 else str(n_lats).zfill(4)
        old_tile = '.old' if nx == 3072 else ''
        tile_file = f'/discover/nobackup/projects/gmao/osse2/stage/BCS_FILES/C{str(nx).strip()}{old_tile}/DC{n_lons_str}xPC{n_lats_str}_CF{nx_str}x6c.bin'
    else:
        if n_lons < 10000:
            nx_str = str(nx).zfill(4)
            n_lons_str = str(n_lons).zfill(4)
            n_lats_str = str(n_lats).zfill(4)
            tile_file = f'/discover/nobackup/ltakacs/bcs/Ganymed-4_0/Ganymed-4_0_Ostia/Shared/DC{n_lons_str}xPC{n_lats_str}_CF{nx_str}x6C.bin'
        else:
            nx_str = str(nx).zfill(4) if nx == 5760 else str(nx).zfill(5)
            n_lons_str = str(n_lons).zfill(4)
            n_lats_str = str(n_lats).zfill(4)
            extension = 'bin' if nx == 5760 else 'til'
            tile_file = f'/discover/nobackup/projects/gmao/osse2/stage/BCS_FILES/C{str(nx).strip()}/DC{n_lons_str}xPC{n_lats_str}_CF{nx_str}x6C.{extension}'
    return tile_file