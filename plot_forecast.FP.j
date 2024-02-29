#!/bin/csh -f
#SBATCH --time=4:00:00
#SBATCH --job-name=PLOT_FP
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH -A s1062
#SBATCH --partition=preops --qos=benchmark
##SBATCH --partition=scutest --constraint=cas --account=k3002
#SBATCH --ntasks=7

source /home/wputman/.cshrc-new
set PROJDIR = /discover/nobackup/projects/gmao/g6dev/
set PROJDIR2 = /discover/nobackup/projects/gmao/osse2/
set  MAIN_DIR = /discover/nobackup/projects/gmao/osse2/IDL_FP
# set MAIN_DIR = /discover/nobackup/projects/gmao/g6dev/pub/qcambrel
set PYTHON_DIR = /discover/nobackup/qcambrel/geospy
hostname
cd $MAIN_DIR
pwd
source /usr/share/modules/init/tcsh
# if ( -e /etc/os-release ) then
#   source idlstart_sles12
# else
  # source idlstart
# source ml use -a /home/mathomp4/modulefiles-SLES12
# source ml python/MINIpyD
source /discover/nobackup/qcambrel/newpythonstart
endif
set host = `hostname`
@ n = 0
module list

# ./plot_forecast.FP.j f5271_fp 12KM plotall_precsnow all 20211207_00z

# arguments: 
set LANG    = $1
set FPTAG   = $2
set RESTAG  = $3
set PLOT    = $4
set REGION  = $5
set FDATE   = $6 
                 set PDATE = ''
if ($#argv >= 6) set PDATE = $6

set ANIMATE_ONLY = 0
if ($PDATE == 'ANIMATE') then
  set PDATE = ''
  set ANIMATE_ONLY = 1
endif
if ($host == 'NOdiscover12') set ANIMATE_ONLY = 1

set PROCESSES = 12
set   MINSIZE = 750000
set OVERWRITE = 0

                    set VIEW = 0
if ($?PBS_NODEFILE) set VIEW = 0
if ($#argv == 7)    set VIEW = 1
if ($VIEW == 1) set OVERWRITE = 1

echo $PROCESSES
echo $OVERWRITE

if (${FPTAG} == 'GFS') then
 set   ETAG = "NCEP"
 set   FTAG = "${FPTAG}"
 set STREAM = "portal"
 set PUBTAG = "${FPTAG}"
 set DATADIR = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/"
 set fcstdate = `echo $FDATE | cut -c1-8`
 set fcsthour = `echo $FDATE | cut -c10-11`
 set interval = 3600
 echo $fcstdate $fcsthour
 set fYear = `echo $FDATE | cut -c1-4`
 set fMon  = `echo $FDATE | cut -c5-6`
 set fDay  = `echo $FDATE | cut -c7-8`
 set fHour = `echo $FDATE | cut -c10-11`
 @ nframes = 240
else
 set   ETAG = "G5GMAO"
 set   FTAG = "${FPTAG}"
 set STREAM = "inst1_2d_asm_Mx"
 set PUBTAG = "${FPTAG}"
endif

set TAG = ${ETAG}-${FTAG}-${PLOT}-${REGION}

mkdir ${MAIN_DIR}/logs/LOGS_${FDATE}
set fFile = "${MAIN_DIR}/logs/LOGS_${FDATE}/.${TAG}.$host.file"
set lFile = "${MAIN_DIR}/logs/LOGS_${FDATE}/.${TAG}.$host.locs"
echo $fFile
echo $lFile
/bin/rm ${fFile}
/bin/rm ${lFile}

set PROJDIR2 = "/discover/nobackup/projects/gmao/osse2"
set GEOSUTIL = "$PROJDIR2/GEOS/DYAMOND/GEOSagcm/src/GMAO_Shared/GEOS_Util"

if ($STREAM != 'portal') then

if ($FDATE == 'ANALYSIS') then

set interval = 3600
set fYear = `echo $PDATE | cut -c1-4`
set pYear = $fYear
set fMon  = `echo $PDATE | cut -c5-6`
set pMon = $fMon
if ($fMon == 'MM') set pMon = '*'
if ($fMon == 'MM') set fMon = '08'
set fDay  = `echo $PDATE | cut -c7-8`
set pDay = $fDay
if ($fDay == 'DD') set pDay = '*'
if ($fDay == 'DD') set fDay = '01'
set fHour = `echo $PDATE | cut -c10-11`
set pHour = $fHour
if ($fHour == 'HH') set pHour = '*'
if ($fHour == 'HH') set fHour = '00'
set fcstdate = "${fYear}${fMon}${fDay}"
set fcsthour = "${fHour}"
echo $fcstdate $fcsthour
echo "/archive/dao_ops/$FPTAG/diag/Y${pYear}/M${pMon}/${FPTAG}.tavg1_2d_slv_Nx.${pYear}${pMon}${pDay}_${pHour}*z.nc4"
set nframes = `ls /archive/dao_ops/$FPTAG/diag/Y${pYear}/M${pMon}/${FPTAG}.tavg1_2d_slv_Nx.${pYear}${pMon}${pDay}_${pHour}*z.nc4 | grep -c nc4`
echo "/home/dao_ops/$FPTAG/run/.../pub/$FPTAG/das/Y${pYear}/M${pMon}/D${pDay}/GEOS.fp.asm.tavg1_2d_slv_Nx.${pYear}${pMon}${pDay}_${pHour}*.nc4"
set nframes = `ls /home/dao_ops/$FPTAG/run/.../pub/$FPTAG/das/Y${pYear}/M${pMon}/D${pDay}/GEOS.fp.asm.tavg1_2d_slv_Nx.${pYear}${pMon}${pDay}_${pHour}*.nc4 | grep -c nc4`
set DATADIR =    "/home/dao_ops/$FPTAG/run/.../pub/$FPTAG/das"
echo $nframes
echo $DATADIR

else

set fcstdate = `echo $FDATE | cut -c1-8`
set fcsthour = `echo $FDATE | cut -c10-11`
set interval = 3600
if (${PLOT} == 'plotall_ir8') set interval = 1800
if (${PLOT} == 'plotall_radar') set interval = 1800
echo $fcstdate $fcsthour
set fYear = `echo $FDATE | cut -c1-4`
set fMon  = `echo $FDATE | cut -c5-6`
set fDay  = `echo $FDATE | cut -c7-8`
set fHour = `echo $FDATE | cut -c10-11`
set nframes = `ls /discover/nobackup/projects/gmao/gmao_ops/pub/$FPTAG/forecast/Y${fYear}/M${fMon}/D${fDay}/H${fHour}/GEOS.fp.fcst.tavg1_2d_slv_Nx.${fYear}${fMon}${fDay}_${fHour}+*.V01.nc4 | grep -c nc4`
set DATADIR =    "/discover/nobackup/projects/gmao/gmao_ops/pub/$FPTAG/forecast"
if ($nframes == 0) then
set nframes = `ls /home/dao_ops/$FPTAG/run/.../pub/$FPTAG/forecast/Y${fYear}/M${fMon}/D${fDay}/H${fHour}/GEOS.fp.fcst.tavg1_2d_slv_Nx.${fYear}${fMon}${fDay}_${fHour}+*.V01.nc4 | grep -c nc4`
set DATADIR =    "/home/dao_ops/$FPTAG/run/.../pub/$FPTAG/forecast"
endif
if ($nframes == 0) then
set nframes = `ls /archive/dao_ops/$FPTAG/forecast/Y${fYear}/M${fMon}/D${fDay}/H${fHour}/GEOS.fp.fcst.tavg1_2d_slv_Nx.${fYear}${fMon}${fDay}_${fHour}+*.V01.nc4 | grep -c nc4`
set DATADIR =    "/archive/dao_ops/$FPTAG/forecast"
endif

endif
endif

echo $nframes
echo $DATADIR

if (${PLOT} == 'stats_120hour') @ nframes = 1

if ( $REGION == 'globe' ) then
  set REGs = ("globe" "goeseast_proj" "goeswest_proj" "meteosat8_proj" "meteosat10_proj" "himawari_proj" "nh_proj" "sh_proj" "northamerica_proj" "westpacific_proj")
else if ( $REGION == 'all' ) then
  set REGs = ("usa_mapset" "centralusa_mapset" "midatlantic_mapset" "northatlantic_mapset" "maryland_mapset" "epacific_mapset" "europe_mapset" "asia_mapset" "australia_mapset" "africa_mapset" "southamerica_mapset" "northamerica_mapset" "indianocean_mapset" "westatlantic_mapset")
else if ( $REGION == 'usa' ) then
  set REGs = ("usa_mapset" "centralusa_mapset" "midatlantic_mapset" "maryland_mapset")
else
  set REGs = ($REGION)
endif

@ frameNo = 0
#if ($PLOT == 'prectotal') @ frameNo = 1
while ($frameNo <= $nframes)

@ fhour = $frameNo * $interval
set framedate = `$GEOSUTIL/post/tick $fcstdate ${fcsthour}0000 $fhour`
echo $framedate

set YEAR    = `echo $framedate | cut -c01-04`
set MONTH   = `echo $framedate | cut -c05-06`
set DAY     = `echo $framedate | cut -c07-08`
set HOUR    = `echo $framedate | cut -c10-11`
set MINUTE  = `echo $framedate | cut -c12-13`
set dateseg = "${YEAR}${MONTH}${DAY}_${HOUR}${MINUTE}"

 # Get DATE information
  set sYR    = `echo $dateseg | cut -c1-4`
  set sMON   = `echo $dateseg | cut -c5-6`
  set sDAY   = `echo $dateseg | cut -c7-8`
  set sHOUR  = `echo $dateseg | cut -c10-11`
  set sMIN   = `echo $dateseg | cut -c12-13`

 # Get IMAGE file information
  if ($PLOT =~ 'plot_storm*') then
   set regionCode = "'${REGION}'"
   set prjBase = "${PLOT}-${REGION}_${ETAG}-${FTAG}"
  else
   if ($LANG == 'python') then
    set regionCode = `grep -B 1 "${REGION}" ${PYTHON_DIR}/regions.json`
    set regionCode = `echo $regionCode[1] | sed 's/:/ /g' | sed 's/"/ /g' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
   else
    set regionCode = `/usr/bin/grep -B 1 "'${REGION}'" ~/IDL_BASE/setup_region.pro`
    set regionCode = `echo $regionCode[1] | sed "s/:/ /g"`
   endif
   set regionCode = $regionCode[1]
   set prjBase = "${PLOT}-${REGION}_${ETAG}-${FTAG}"
  endif
  set pngFile = "${prjBase}-${fYear}${fMon}${fDay}_${fHour}z_${dateseg}.png"


 # Remove any local copy of the image
 #set removeLocal = 1
 #if ($removeLocal) /bin/rm -rf $pngFile

 # Make pub dir available
  set uPLOT = `echo $PLOT | tr "[a-z]" "[A-Z]" `
  set PUB_DIR = $PROJDIR/pub/${PUBTAG}/$uPLOT
  set pubLoc  = "$PUB_DIR/Y$sYR/M$sMON/D$sDAY"
  if (! -e $pubLoc) then
     mkdir -p $pubLoc
     chmod 755 ${PUB_DIR}
     chmod 755 ${PUB_DIR}/Y${sYR}
     chmod 755 ${PUB_DIR}/Y${sYR}/M${sMON}
     chmod 755 ${PUB_DIR}/Y${sYR}/M${sMON}/D${sDAY}
  endif

if (! $ANIMATE_ONLY) then

# Now done in cron jobs
## clear old images
# set nDAYS = `ls -ld $PUB_DIR/Y*/M*/D* | grep -c D`
# if ($nDAYS > 46) then
#   echo $nDAYS
#   set DAY_DIRS = `ls -d $PUB_DIR/Y*/M*/D*`
#   @ nDAY = 1
#   foreach DIR ($DAY_DIRS)
#     if ($nDAY < $nDAYS - 45) then
#        echo $DIR
#        /bin/rm -rf $DIR &
#     endif
#     @ nDAY++
#   end
# endif

 # Check if IMAGE exists and is the proper size
  set doImage = 0
  set pubFile = "$pubLoc/$pngFile"
  set pngFile0 = $pngFile
  set pubFile0 = $pubFile 
  foreach REG ($REGs) 
  if ($REGION == 'globe') then
  set pngFile = `echo $pngFile0 | sed "s/-globe_/-${REG}_/g"`
  set pubFile = `echo $pubFile0 | sed "s/-globe_/-${REG}_/g"`
  endif
  if ($REGION == 'all') then
  set pngFile = `echo $pngFile0 | sed "s/-all_/-${REG}_/g"`
  set pubFile = `echo $pubFile0 | sed "s/-all_/-${REG}_/g"`
  endif
  if ($REGION == 'usa') then
  set pngFile = `echo $pngFile0 | sed "s/-usa_/-${REG}_/g"`
  set pubFile = `echo $pubFile0 | sed "s/-usa_/-${REG}_/g"`
  endif
  if (! -e $pubFile) then
      echo $pubFile
      if (! -e $pngFile) set doImage = 1
  else
    /bin/ls -l  $pubFile
    set size = `stat -c %s $pubFile`
    if ($size < $MINSIZE) then
       set doImage = 1
    endif
  endif
  set pngFile = $pngFile0
  set pubFile = $pubFile0
  end

  if ($OVERWRITE) set doImage = 1
  echo $dateseg $PDATE
  if ( ($#argv >= 6) && ($dateseg != $PDATE) ) set doImage = 0
  echo $#argv $dateseg $PDATE
  if ( ($PLOT == 'plot_anom') && (($sMIN != '30') || (${regionCode} != 1)) ) set doImage = 0
  if ( ($PLOT == 'precrain') && ($sMIN != '00') ) set doImage = 0

  echo $doImage
  if ($doImage) then
    echo $pubFile

    if ($VIEW) then
        echo "${PLOT}, '$sYR', '$sMON', '$sDAY', '$sHOUR', '$sMIN', '${ETAG}', '${FTAG}-${fYear}${fMon}${fDay}_${fHour}z', '${DATADIR}', '${STREAM}', FDATE='${FDATE}', REGION=${regionCode}"
        if $LANG == 'idl'
         idl -e "${PLOT}, '$sYR', '$sMON', '$sDAY', '$sHOUR', '$sMIN', '${ETAG}', '${FTAG}-${fYear}${fMon}${fDay}_${fHour}z', '${DATADIR}', '${STREAM}', FDATE='${FDATE}', REGION=${regionCode}"
      # python -e "${PLOT}, '$sYR', '$sMON', '$sDAY', '$sHOUR', '$sMIN', '${ETAG}', '${FTAG}-${fYear}${fMon}${fDay}_${fHour}z', '${DATADIR}', '${STREAM}', FDATE='${FDATE}', REGION=${regionCode}"
        else
         python $PLOT $sYR $sMON $sDAY $sHOUR $sMIN $ETAG '${FTAG}-${fYear}${fMon}${fDay}_${fHour}z' $DATADIR $STREAM --f_date=$FDATE --region=$regionCode
        endif
    else
      if $LANG == 'idl'
       idl -e "${PLOT}, '$sYR', '$sMON', '$sDAY', '$sHOUR', '$sMIN', '${ETAG}', '${FTAG}-${fYear}${fMon}${fDay}_${fHour}z', '${DATADIR}', '${STREAM}', FDATE='${FDATE}', REGION=${regionCode}" >& /dev/null & #.${PLOT}_${dateseg}.log  &
      # python -e "${PLOT}, '$sYR', '$sMON', '$sDAY', '$sHOUR', '$sMIN', '${ETAG}', '${FTAG}-${fYear}${fMon}${fDay}_${fHour}z', '${DATADIR}', '${STREAM}', FDATE='${FDATE}', REGION=${regionCode}" >& /dev/null & #.${PLOT}_${dateseg}.log  &
      else
       python $PLOT $sYR $sMON $sDAY $sHOUR $sMIN $ETAG '${FTAG}-${fYear}${fMon}${fDay}_${fHour}z' $DATADIR $STREAM --f_date=$FDATE --region=$regionCode >& /dev/null
      endif
    endif
    @ n++

    foreach REG ($REGs)
      set file2 = `echo $pngFile | sed "s/globe/$REG/g"`
      echo ${file2} >> ${fFile}
      echo ${pubLoc} >> ${lFile}
    end

  endif

  if ($n == $PROCESSES) then
  wait
  @ n = 0
  @ i = 1
  set ifiles = `cat ${fFile}`
  set ilocs  = `cat ${lFile}`
  while ($i <= $#ifiles)
     set iloc  = $ilocs[$i]
     set ifile = $ifiles[$i]
     if (-e ${ifile}) then
       echo ${ifile} $iloc 
       /bin/mv ${ifile} $iloc
       chmod 755 $iloc/${ifile}
     endif
     @ i = $i + 1
  end
  if (-e ${fFile}) /bin/rm ${fFile}
  if (-e ${lFile}) /bin/rm ${lFile}
  endif
  @ f++

  if ($doImage) then
     if ($VIEW) exit 0
  endif

@ frameNo++
end
wait
if (-e ${fFile}) then
  @ i = 1
  set ifiles = `cat ${fFile}`
  set ilocs  = `cat ${lFile}`
  while ($i <= $#ifiles)
     set iloc  = $ilocs[$i]
     set ifile = $ifiles[$i]
     echo ${ifile} $iloc
     if (-e ${ifile}) then
       /bin/mv ${ifile} $iloc
       chmod 755 $iloc/${ifile}
     endif
     @ i = $i + 1
  end
  if (-e ${fFile}) /bin/rm ${fFile}
  if (-e ${lFile}) /bin/rm ${lFile}
  @ n = 0
endif

endif # ANIMATE_ONLY

if (1) then 
set prjBase0 = ${prjBase}
@ m = 0
foreach reg ($REGs)
set prjBase = `echo $prjBase0 | sed "s/-${REGION}_/-${reg}_/g"`
echo ANIMATE ${prjBase}
if ( -e /etc/os-release ) then
 if ($ANIMATE_ONLY) then
 /bin/rm -rf ANIMATE_LOGS/*log
 mkdir ANIMATE_LOGS
 ./animate_frames_fp_sles12.j "$RESTAG" "$REGION" "$FDATE" "$FPTAG" "$prjBase" "$PUB_DIR" #</dev/null  #>& /dev/null & #>& ANIMATE_LOGS/${prjBase}.log &
 else
 sbatch animate_frames_fp_sles12.j "$RESTAG" "$REGION" "$FDATE" "$FPTAG" "$prjBase" "$PUB_DIR" 
 endif
else
 ./animate_frames_fp.j "$RESTAG" "$REGION" "$FDATE" "$FPTAG" "$prjBase" "$PUB_DIR" >& /dev/null &
endif
@ m++
#if ($m == $PROCESSES) then
if ($m == 1) then
  wait
  @ m = 0
endif
end # Loop over regions
wait
endif

if (0) then
if ( ($regionCode == 1) || ($regionCode == 0) ) then
if ( (${PLOT} == 'plot_olr') || (${PLOT} == 'plot_ir') || (${PLOT} == 'plot_cldvis') ) then
 convert  -resize 20%  -delay 20   -loop 0   $PUB_DIR/Y*/M*/D*/${prjBase}-${fYear}${fMon}${fDay}_${fHour}z*20170821_16*.png \
                                             $PUB_DIR/Y*/M*/D*/${prjBase}-${fYear}${fMon}${fDay}_${fHour}z*20170821_17*.png \
                                             $PUB_DIR/Y*/M*/D*/${prjBase}-${fYear}${fMon}${fDay}_${fHour}z*20170821_18*.png \
                                             $PUB_DIR/Y*/M*/D*/${prjBase}-${fYear}${fMon}${fDay}_${fHour}z*20170821_19*.png \
                                             $PUB_DIR/Y*/M*/D*/${prjBase}-${fYear}${fMon}${fDay}_${fHour}z*20170821_20*.png \
  eclipse_${prjBase}-${fYear}${fMon}${fDay}_${fHour}z.gif
  mkdir -p                                                        $PROJDIR/pub/MOVIES_NWP_${FPTAG}/${FDATE}
  /bin/mv eclipse_${prjBase}-${fYear}${fMon}${fDay}_${fHour}z.gif $PROJDIR/pub/MOVIES_NWP_${FPTAG}/${FDATE}
endif
endif
endif

exit 0

