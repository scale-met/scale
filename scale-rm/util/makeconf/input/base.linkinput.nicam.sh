#!/bin/bash

TARGET_FILE="base.linkinput"

rm -f ${TARGET_FILE}
echo "#!/bin/bash" > ${TARGET_FILE}

# convert integer from string (TIME_DURATION)
LINK_DURATION=`echo ${TIME_DURATION}| cut -d \. -f 1`

# not to need fix if the unit is SEC
if   [ ${TIME_UNIT} = "MIN"  ]; then
  LINK_DURATION=`expr ${LINK_DURATION} \* 60`
elif [ ${TIME_UNIT} = "HOUR" ]; then
  LINK_DURATION=`expr ${LINK_DURATION} \* 3600`
elif [ ${TIME_UNIT} = "DAY"  ]; then
  LINK_DURATION=`expr ${LINK_DURATION} \* 86400`
fi

start_unix_sec=`date -u -d "${STARTDATE[0]}-${STARTDATE[1]}-${STARTDATE[2]} ${STARTDATE[3]}:${STARTDATE[4]}:${STARTDATE[5]}" +%s`
end_unix_sec=`expr ${start_unix_sec} + ${LINK_DURATION}`

fn=0
unix_sec=${start_unix_sec}
while [ ${unix_sec} -le `expr ${end_unix_sec} + ${LINK_BOUNDARY_DT}` ];
do
  year=`date -u -d "@${unix_sec}" +%Y`
  month=`date -u -d "@${unix_sec}" +%m`
  day=`date -u -d "@${unix_sec}" +%d`
  hour=`date -u -d "@${unix_sec}" +%H`
  min=`date -u -d "@${unix_sec}" +%M`
  sec=`date -u -d "@${unix_sec}" +%S`

  fmtd_fn=`printf "%05d" $fn`

  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/la_tg.nc      ${LINK_OUTPUTDIR}/la_tg_${fmtd_fn}.${LINK_FTYPE}.nc     " >> ${TARGET_FILE}
  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/la_wg.nc      ${LINK_OUTPUTDIR}/la_wg_${fmtd_fn}.${LINK_FTYPE}.nc     " >> ${TARGET_FILE}
  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/lsmask.nc     ${LINK_OUTPUTDIR}/lsmask_${fmtd_fn}.${LINK_FTYPE}.nc    " >> ${TARGET_FILE}
  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/ms_pres.nc    ${LINK_OUTPUTDIR}/ms_pres_${fmtd_fn}.${LINK_FTYPE}.nc   " >> ${TARGET_FILE}
  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/ms_qv.nc      ${LINK_OUTPUTDIR}/ms_qv_${fmtd_fn}.${LINK_FTYPE}.nc     " >> ${TARGET_FILE}
  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/ms_tem.nc     ${LINK_OUTPUTDIR}/ms_tem_${fmtd_fn}.${LINK_FTYPE}.nc    " >> ${TARGET_FILE}
  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/ms_u.nc       ${LINK_OUTPUTDIR}/ms_u_${fmtd_fn}.${LINK_FTYPE}.nc      " >> ${TARGET_FILE}
  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/ms_v.nc       ${LINK_OUTPUTDIR}/ms_v_${fmtd_fn}.${LINK_FTYPE}.nc      " >> ${TARGET_FILE}
  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/oa_sst.nc     ${LINK_OUTPUTDIR}/oa_sst_${fmtd_fn}.${LINK_FTYPE}.nc    " >> ${TARGET_FILE}
  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/ss_q2m.nc     ${LINK_OUTPUTDIR}/ss_q2m_${fmtd_fn}.${LINK_FTYPE}.nc    " >> ${TARGET_FILE}
  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/ss_slp.nc     ${LINK_OUTPUTDIR}/ss_slp_${fmtd_fn}.${LINK_FTYPE}.nc    " >> ${TARGET_FILE}
  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/ss_t2m.nc     ${LINK_OUTPUTDIR}/ss_t2m_${fmtd_fn}.${LINK_FTYPE}.nc    " >> ${TARGET_FILE}
  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/ss_tem_sfc.nc ${LINK_OUTPUTDIR}/ss_tem_sfc_${fmtd_fn}.${LINK_FTYPE}.nc" >> ${TARGET_FILE}
  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/ss_u10m.nc    ${LINK_OUTPUTDIR}/ss_u10m_${fmtd_fn}.${LINK_FTYPE}.nc   " >> ${TARGET_FILE}
  echo "ln -svf ${LINK_BASEDIR}/${year}${month}${day}/ss_v10m.nc    ${LINK_OUTPUTDIR}/ss_v10m_${fmtd_fn}.${LINK_FTYPE}.nc   " >> ${TARGET_FILE}

  fn=`expr ${fn} \+ 1`
  unix_sec=`expr ${unix_sec} \+ ${LINK_BOUNDARY_DT}`
done
