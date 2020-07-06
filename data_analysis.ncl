;*************************************************
;load "statistics2.ncl"
;load "select_region.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
;************************************************
undef("coswgt_area_average")
function coswgt_area_average(field:numeric)
;
; Arguments
;    field: numeric field.
;
; Return value
;    The area average using cosine lat weights.
;
; Description
;    Computes the area average using cosine lat weights and lon weights=1.
;
local funcname, scriptname, lat, wgt_lat, lon, lon_size, wgt_lon, ave
begin

  funcname = "coswgt_area_average"
  scriptname = "diag_scripts/shared/statistics.ncl"

  lat = field&lat
  wgt_lat = tofloat(NormCosWgtGlobe(lat))

  lon = field&lon
  lon_size = dimsizes(lon)
  wgt_lon = new((/lon_size(0)/), float)
  wgt_lon = 1.0

  ave = wgt_areaave_Wrap(field, wgt_lat, wgt_lon, 0)
  return(ave)

end
; #############################################################################
undef("coswgt_area_std")
function coswgt_area_std(field:numeric)
;
; Arguments
;    field: numeric field.
;
; Return value
;    The area average using cosine lat weights.
;
; Description
;    Computes the area average using cosine lat weights and lon weights=1.
;
begin

  funcname = "coswgt_area_std"

   lat = field&lat
   lon = field&lon
   re   = 6.37122e06
   rad  = 4.0 * atan(1.0) / 180.0
   con  = re * rad
   clat = cos(lat * rad)           ; cosine of latitude
   dlon = (lon(2) - lon(1))        ; assume dlon is constant
   dlat = (lat(2) - lat(1))        ; assume dlat is constant

   dx   = con * dlon * clat        ; dx at each latitude
   dy   = con * dlat               ; dy is constant
   dydx = dy * dx                  ; dydx(nlat)
   nlat=dimsizes(lat)
   mlon=dimsizes(lon)
   wgt  = new((/nlat, mlon/), typeof(dydx))
   wgt=conform (wgt, dydx, 0)



   ave = wgt_areaave2(field, wgt, 0)
   d2 = (field - ave) ^ 2 ;
     ; Calculate weighted standard deviation
   ss=sum(d2)
   arg = wgt_areasum2(d2, wgt, 0)

   v1 = sum(where(ismissing(field), wgt@_FillValue, wgt))
   v2 = sum(where(ismissing(field), wgt@_FillValue, wgt) ^ 2)
   dd=(v1 / (v1 ^ 2 - v2))

   out = sqrt( arg*dd)

  return(out)


end

begin
;************************************************
; read in netCDF file
;************************************************
dir="/Users/irenecionni/" 
a = addfile(dir+"/Delverables_C3S_511/ERA5/grid_o3_197901-202005.nc" ,"r")
;************************************************
; read in zonal winds
;************************************************
  ;vmro3=a->o3(:491,:,:,:)      
  
  vmro31=a->o3  (:491,:,:,:)                         
  vmro3=vmro31
  vmro3=(/vmro31*1e6*28.9644 / 47.9982 /)
  vmro3!1="plev"
  plev=a->level
  plev@long_name="Pressure"
  plev@units="hPa"
  vmro3&plev=plev
  lat=a->lat

 
  lon=a->lon
  dim=dimsizes(vmro3)
  yyyymm=  cd_calendar(vmro3&time, 1) 
  time=vmro3&time
  date= cd_calendar(time, 0)
  yy=toint(date(:,0)) 
  ts_avg=coswgt_area_average(vmro3)
  plev=vmro3&plev
  ts1_avg=new((/dimsizes(ts_avg&time),6/),typeof(ts_avg))
  ts1_avg(:,0)=ts_avg(:,{2.})
  ts1_avg(:,1)=ts_avg(:,{10.})
  ts1_avg(:,2)=ts_avg(:,{70.})   
  ts1_avg(:,3)=ts_avg(:,{250.})
  ts1_avg(:,4)=ts_avg(:,{500.})
  ts1_avg(:,5)=ts_avg(:,{1000.})
  map_lev=new((/dimsizes(lat),dimsizes(lon),6/),typeof(ts_avg))
  map_lev(:,:,0)=dim_avg_n_Wrap(vmro3(:,{2.},:,:),0)
  map_lev(:,:,1)=dim_avg_n_Wrap(vmro3(:,{10.},:,:),0)
  map_lev(:,:,2)=dim_avg_n_Wrap(vmro3(:,{70.},:,:),0)
  map_lev(:,:,3)=dim_avg_n_Wrap(vmro3(:,{250.},:,:),0)
  map_lev(:,:,4)=dim_avg_n_Wrap(vmro3(:,{500.},:,:),0)
  map_lev(:,:,5)=dim_avg_n_Wrap(vmro3(:,{1000.},:,:),0)
   
  map_lev!0="lat"
  map_lev!1="lon"
  map_lev&lat=vmro3&lat
  map_lev&lon=vmro3&lon
  map_lev!2="lev"
  avg_clim1=new(6,typeof(map_lev))
  std_clim1=new(6,"double")

  do i=0,5
    std_clim1(i)=coswgt_area_std(map_lev(:,:,i)) 
    avg_clim1(i)=coswgt_area_average(map_lev(:,:,i))
  end do

  ts1_lat=new((/dimsizes(vmro3&time),dimsizes(lat),6/),typeof(ts_avg))
  ao=vmro3(:,{2.},:,:)
  ts1_lat(:,:,0)=dim_avg_Wrap(vmro3(:,{2.},:,:))
  ts1_lat(:,:,1)=dim_avg_Wrap(vmro3(:,{10.},:,:))
  ts1_lat(:,:,2)=dim_avg_Wrap(vmro3(:,{70.},:,:))
  ts1_lat(:,:,3)=dim_avg_Wrap(vmro3(:,{250.},:,:))
  ts1_lat(:,:,4)=dim_avg_Wrap(vmro3(:,{500.},:,:))
  ts1_lat(:,:,5)=dim_avg_Wrap(vmro3(:,{1000.},:,:))
  ts1_lat!0="time"
  ts1_lat!1="lat"
  ts1_lat!2="lev"
  ts1_lat&time=vmro3&time
  ts1_lat&lat=vmro3&lat
  ts1_lat_avg=dim_avg_n_Wrap(ts1_lat,0)
  lev_label=(/"2 hPa","10 hPa","70 hPa","250 hPa","500 hPa","1000 hPa"/)

  ;regions
  region_label=(/"Global","Trop","NH Midlatitudes","SH Midlatitudes","NH Polar","SH Polar"/)
  region_max=(/90.,10.,60.,-35.,90.,-60./)
  region_min=(/-90.,-10.,35.,-60.,60.,-90./)
  lon_vmro3=dim_avg_Wrap(vmro3)

  ts_regions=new((/dimsizes(ts_avg&time),dimsizes(plev),dimsizes(region_label)/),typeof(ts_avg))
  ts_regions1=new((/dimsizes(ts_avg&time),dimsizes(plev),dimsizes(region_label)/),typeof(ts_avg))
  do ir=0,dimsizes(region_label)-1
     b=coswgt_area_average(vmro3(:,:,{region_min(ir):region_max(ir)},:))
     ts_regions1(:,:,ir)=b
     b_mean=dim_avg_n_Wrap(b,0)
     do ilev=0,dimsizes(plev)-1
        ts_regions(:,ilev,ir)=(/b(:,ilev)-b_mean(ilev)/)
     end do   
  end do
  ts_regions!0="time"
  ts_regions!1="plev"
  ts_regions!2="region"

  ts_regions&time=ts_avg&time
  ts_regions&plev=plev

  ts_region_avg=new((/3,dimsizes(plev),dimsizes(region_label)/),typeof(ts_avg))
  ts_region_avg(0,:,:)=dim_avg_n_Wrap(ts_regions1(:,:,:),0)
  ts_region_avg(1,:,:)=dim_avg_n_Wrap(ts_regions1(0:179,:,:),0)
  ts_region_avg(2,:,:)=dim_avg_n_Wrap(ts_regions1(180:,:,:),0)

  seas=(/"ANN","ANN","ANN","ANN","ANN","ANN"/)


  ts_region_avg!0="line"
  ts_region_avg!1="plev"
  ts_region_avg!2="region"
  ts_region_avg&plev=plev

;***********************time serie************************************************
  plot=new(12,"graphic")
  wks1 = gsn_open_wks("ps","ts_region_plev_time"+"CDS-ERA5") 
  gsn_define_colormap(wks1,"MPL_RdYlGn")
  gsn_reverse_colormap(wks1)
  k=0
  do ir=0,5
     res2 = True                                    ; plot mods desired
     res2@gsnDraw                  = False            ; don't draw yet
     res2@gsnFrame                 = False            ; don't advance frame yet

     res2@cnFillOn         = True               ; color Fill
     res2@cnFillMode       = "RasterFill"       ; Raster Mode
     res2@cnLinesOn        =  False             ; Turn off contour lines
     res2@cnLineLabelsOn   =  False
     res2@cnLinesOn        =  False
     res2@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res2@cnMinLevelValF       = 0.05               ; set min contour level
     res2@cnMaxLevelValF       = 0.3               ; set max contour level
     res2@cnLevelSpacingF      = 0.05               ; set contour spacing
     res2@cnLevelSelectionMode = "ExplicitLevels"
     res2@lbLabelBarOn = False
     res2@tmYRMode             = "Automatic"
     res2@tiMainString         = region_label(ir)    ; plot title
     fmt = "%3.2f"
     res2@gsnRightString     = "[ppmv]"
     res2@vpHeightF= 0.4                    ; change aspect ratio of plot
     res2@vpWidthF = 0.6
     res2@cnInfoLabelOn=False
     resTick                = True
     resTick@ttmFormat      = "%Y"
     resTick@ttmAxis        = "XB"
     resTick@ttmMajorStride = 36
     resTick@ttmMinorStride = 12
     time_axis_labels( ts_regions&time,  res2, resTick )

     res2@tmYLMode        = "Explicit"
     res2@tmYLValues      = (/1000.,700.,500.,250.,150.,100.,70.,50.,30.,10./)
     res2@tmYLLabels      = (/"1000.","700.","500.","250","150","100","70","50","30","10"/) 
     ts_regions&plev=plev
     plot(k) =gsn_csm_pres_hgt(wks1,ts_regions(plev|:,time|:,region|ir), res2)  ; create plot
     k=k+1
     res21=True
     res21@gsnDraw                  = False            ; don't draw yet
     res21@gsnFrame                 = False            ; don't advance frame yet
     res21@vpHeightF= 0.4                    ; change aspect ratio of plot
     res21@vpWidthF = 0.3
     res21@tiXAxisString = "O3 [ppmv]"                 ; x axis title
     res21@trXMinF      = 0.
     res21@trXMaxF      = 10.
     res21@tiMainString         = seas(ir)
     res21@tmYLMode        = "Explicit"
     res21@tmYLValues      = (/1000.,700.,500.,250,150,100,70,50,30,10/)
     res21@tmYLLabels      = (/"1000.","700.","500.","250","150","100","70","50","30","10"/)
     res21@tmXBMode        = "Explicit"
     res21@tmXBValues      = (/0.,1.,2.,3.,4.,5.,6.,7.,8.,9./)
     res21@tmXBLabels      = (/"0.","1.","2.","3.","4.","5.","6.","7.","8.","9."/)
     res21@tmXBLabelFontHeightF=0.015
     res21@xyYStyle          = "Log"
     res21@trYReverse =True
     res21@xyLineColors     = (/"NavyBlue"/)
     res21@xyLineThickness = 3
     plot(k) = gsn_csm_xy(wks1,ts_region_avg(line|0,plev|:,region|ir),ts_region_avg&plev, res21)  ; create plot
     k=k+1
  end do

  resP=True
  resP@gsnPanelLabelBar    = True
  resP@gsnPanelXWhiteSpacePercent =0.5 
  gsn_panel(wks1,plot,(/3,4/),resP)
;***********************time serie 2************************************************
  delete(res2)
  delete(res21)
  plota=new(12,"graphic")
  wksa = gsn_open_wks("ps","ts2_region_plev_time"+"ERA5")
  gsn_define_colormap(wksa,"MPL_RdYlGn")
  gsn_reverse_colormap(wksa)
  k=0
  do ir=0,5
     res2 = True                                    ; plot mods desired
     res2@tiMainString         = region_label(ir)    ; plot title
     fmt = "%3.2f"
     res2@gsnDraw                  = False            ; don't draw yet
     res2@gsnFrame                 = False            ; don't advance frame yet
     res2@vpHeightF= 0.4                    ; change aspect ratio of plot
     res2@vpWidthF = 0.6
     res2@tiYAxisString = "O3 [ppmv]"                 ; x axis title
     res2@tmYLValues      = (/0.,1.,2.,3.,4.,5.,6.,7.,8.,9./)
     res2@tmYLLabels      = (/"0.","1.","2.","3.","4.","5.","6.","7.","8.","9."/)

     resTick                = True
     resTick@ttmFormat      = "%Y"
     resTick@ttmAxis        = "XB"
     resTick@ttmMajorStride = 36
     resTick@ttmMinorStride = 12
     time_axis_labels( ts_regions&time,  res2, resTick )

     ts_regions&plev=plev
     res2@xyDashPattern  = 0
     res2@xyLineColors     = (/"Black","Green","Red","Blue","Orange"/)
     res2@xyMonoLineThickness =False
     res2@xyLineThickness = 2.    
     time1=ts_regions&time
     res2@trXMinF      = time1(0)
     res2@trXMaxF      = time1(dimsizes(time1)-1)
     if (ir.eq.0)then
        res2@pmLegendDisplayMode    = "Always"            ; turn on legend
     end if

     if (ir.eq.1)then
        delete(res2@pmLegendDisplayMode)
     end if

     res2@pmLegendSide           = "Top"               ; Change location of
     res2@pmLegendParallelPosF   = 0.7                  ; move units right
     res2@pmLegendOrthogonalPosF = -0.6                ; move units down
     res2@pmLegendWidthF         = 0.2                ; Change width and
     res2@pmLegendHeightF        = 0.18                ; height of legend.
     res2@lgPerimOn              = False               ; turn off box around
     res2@lgLabelFontHeightF     = .02                 ; label font height
     res2@xyExplicitLegendLabels = (/"1000 hPa","500 hPa","250 hPa","70 hPa","10 hPa","2 hPa"/)         ;
     sel_lev=(/1000.,500.,250.,70.,10.,2./)
     ta_regionsa=ts_regions({plev|sel_lev},time|:,region|ir)
     ta_regionsa(0,:)=rmAnnCycle1D(ta_regionsa(0,:))
     ta_regionsa(1,:)=rmAnnCycle1D(ta_regionsa(1,:))
     ta_regionsa(2,:)=rmAnnCycle1D(ta_regionsa(2,:))
     ta_regionsa(3,:)=rmAnnCycle1D(ta_regionsa(3,:))
     ta_regionsa(4,:)=rmAnnCycle1D(ta_regionsa(4,:))
     ta_regionsa(5,:)=rmAnnCycle1D(ta_regionsa(5,:))


     plota(k) =gsn_csm_xy(wksa,ts_regions&time,ta_regionsa, res2)  ; create plot
     k=k+1
     res21=True
     res21@gsnDraw                  = False            ; don't draw yet
     res21@gsnFrame                 = False            ; don't advance frame yet
     res21@vpHeightF= 0.4                    ; change aspect ratio of plot
     res21@vpWidthF = 0.3
     res21@tiXAxisString = "O3 [ppmv]"                 ; x axis title
     res21@tiYAxisString = "Pressure (hPa)"
     res21@trYMinF      = 1
     res21@trYMaxF      = 250
     res21@trXMinF      = 0.
     res21@trXMaxF      = 10.
     res21@tiMainString         = seas(ir)+" "+ region_label(ir)
     res21@tmYLMode        = "Explicit"
     res21@tmYLValues      = (/250,200,150,100,70,50,30,10/)
     res21@tmYLLabels      = (/"250","200","150","100","70","50","30","10"/)
     res21@tmXBMode        = "Explicit"
     res21@tmXBValues      = (/0.,1.,2.,3.,4.,5.,6.,7.,8.,9./)
     res21@tmXBLabels      = (/"0.","1.","2.","3.","4.","5.","6.","7.","8.","9."/)
     res21@tmXBLabelFontHeightF=0.015
     res21@xyYStyle          = "Log"
     res21@trYReverse =True
     res21@xyLineColors     = (/"NavyBlue"/)
     res21@xyLineThickness = 3
     plota(k) = gsn_csm_xy(wksa,ts_region_avg(line|0,plev|:,region|ir),ts_region_avg&plev, res21)  ; create plot
     k=k+1
  end do

  resP=True
  resP@gsnPanelLabelBar    = True
  resP@gsnPanelXWhiteSpacePercent =0.1
  gsn_panel(wksa,plota((/0,2,4,6,8,10/)),(/3,2/),resP)
  gsn_panel(wksa,plota((/1,3,5,7,9,11/)),(/3,2/),resP)
;*********èèèè3333333time series++++++++++++


;************************************************
  plot0=new(6,"graphic")
  wks0 = gsn_open_wks("png","map_"+"CDS-ERA5")
  gsn_define_colormap(wks0,"MPL_RdYlGn")
  gsn_reverse_colormap(wks0)
  res = True                                    ; plot mods desired
  res@gsnDraw                  = False            ; don't draw yet
  res@gsnFrame                 = False            ; don't advance frame yet
  res@lbLabelBarOn = True
  res@cnLineLabelsOn   =  False
  res@lbOrientation = "vertical"
  res@cnLinesOn        =  False
  res@cnFillOn         = True               ; color Fill
  res@cnFillMode       = "RasterFill"       ; Raster Mode
  res@cnLinesOn        =  False             ; Turn off contour lines
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = 0.5              ; set min contour level
  res@cnMaxLevelValF       = 5.               ; set max contour level
  res@cnLevelSpacingF      = 0.5               ; set contour spacing
  res@gsnLeftString      = "ppmv"
  res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1(0))+" std:"+sprintf(fmt,std_clim1(0))           ; plot center string
  res@tiMainString         = lev_label(0)    ; plot title
  plot0(0) = gsn_csm_contour_map(wks0,map_lev(:,:,0), res)  ; create plot
  res@tiMainString         = lev_label(1)    ; plot title
  res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1(1))+" std:"+sprintf(fmt,std_clim1(1))   
  res@cnMinLevelValF       = 1.             ; set min contour level
  res@cnMaxLevelValF       = 10.               ; set max contour level
  res@cnLevelSpacingF      = 1.        
  plot0(1) = gsn_csm_contour_map(wks0,map_lev(:,:,1), res)  ; create plot
  res@tiMainString         = lev_label(2)    ; plot title
  res@cnMinLevelValF       = 0.125               ; set min contour level
  res@cnMaxLevelValF       = 2.5               ; set max contour level
  res@cnLevelSpacingF      = 0.125
  res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1(2))+" std:"+sprintf(fmt,std_clim1(2)) 
  plot0(2) = gsn_csm_contour_map(wks0,map_lev(:,:,2), res)  ; create plot
  res@tiMainString         = lev_label(3)    ; plot title
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = 0.025               ; set min contour level
  res@cnMaxLevelValF       = 0.5               ; set max contour level
  res@cnLevelSpacingF      = 0.025               ; set contour spacing
  res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1(3))+" std:"+sprintf(fmt,std_clim1(3)) 
  plot0(3) = gsn_csm_contour_map(wks0,map_lev(:,:,3), res)  ; create plot
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = 0.005              ; set min contour level
  res@cnMaxLevelValF       = 0.1               ; set max contour level
  res@cnLevelSpacingF      = 0.005               ; set contour spacing
  res@tiMainString         = lev_label(4)    ; plot title
  res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1(4))+" std:"+sprintf(fmt,std_clim1(4)) 
  plot0(4) = gsn_csm_contour_map(wks0,map_lev(:,:,4), res)  ; create plot
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = 0.0025               ; set min contour level
  res@cnMaxLevelValF       = 0.05               ; set max contour level
  res@cnLevelSpacingF      = 0.0025               ; set contour spacing
  res@tiMainString         = lev_label(5)    ; plot title
  res@gsnRightString     = "mean: "+sprintf(fmt,avg_clim1(5))+" std:"+sprintf(fmt,std_clim1(5))
  plot0(5) = gsn_csm_contour_map(wks0,map_lev(:,:,5), res)  ; create plot
  resP=True
  resP=True
  resP@gsnPanelLabelBar    = False
  gsn_panel(wks0,plot0,(/3,2/),resP)

  delete(res)
;**********************************************
  plot3=new(12,"graphic")
  wks3 = gsn_open_wks("png","zon_"+"CDS-ERA5")       

  gsn_define_colormap(wks3,"MPL_RdYlGn")
  gsn_reverse_colormap(wks3)
  res = True                                    ; plot mods desired
  res@gsnDraw                  = False            ; don't draw yet
  res@gsnFrame                 = False            ; don't advance frame yet
  res@lbLabelBarOn = True
  res@cnLineLabelsOn   =  False
  res@cnLinesOn        =  False 
  res@cnFillOn         = True               ; color Fill
  res@cnFillMode       = "RasterFill"       ; Raster Mode
  res@cnLinesOn        =  False             ; Turn off contour lines
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = 0.5               ; set min contour level
  res@cnMaxLevelValF       = 5.               ; set max contour level
  res@cnLevelSpacingF      = 0.5               ; set contour spacing

  fmt = "%3.2f"
  res@vpHeightF= 0.4                    ; change aspect ratio of plot
  res@vpWidthF = 0.7
  res@tmYLMode        = "Explicit"
  res@tmYLValues      = (/-60,-30,0,30,60/)
  res@tmYLLabels      = (/"60S","30S","Eq","30N","60N"/)
  res@tmYLLabelFontHeightF=0.015
  res@tmXBLabelFontHeightF=0.015
  res@cnInfoLabelOn=False
  resTick                = True
  resTick@ttmFormat      = "%Y"
  resTick@ttmAxis        = "XB"
  resTick@ttmMajorStride = 36
  resTick@ttmMinorStride = 12
  res@lbOrientation = "vertical" 
  res@gsnRightString     = "[ppmv]"
; Set resources necessary to customize Y axis labels
  time_axis_labels( ts1_lat&time,  res, resTick )
  res@tiMainString         = lev_label(0) 
  res@cnMinLevelValF       = 1.             ; set min contour level
  res@cnMaxLevelValF       = 10.               ; set max contour level
  res@cnLevelSpacingF      = 1.

  plot3(0) = gsn_csm_lat_time(wks3,ts1_lat(lat|:,time|:,lev|0), res)  ; create plot
  res@lbOrientation = "vertical"
  res@tiMainString         = lev_label(1) 
  plot3(2) = gsn_csm_lat_time(wks3,ts1_lat(lat|:,time|:,lev|1), res)  ; create plot
  res@tiMainString         = lev_label(2)
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = 0.125               ; set min contour level
  res@cnMaxLevelValF       = 2.5               ; set max contour level
  res@cnLevelSpacingF      = 0.125
  plot3(4) = gsn_csm_lat_time(wks3,ts1_lat(lat|:,time|:,lev|2), res)  ; create plot
  res@tiMainString         = lev_label(3)
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = 0.025               ; set min contour level
  res@cnMaxLevelValF       = 0.5               ; set max contour level
  res@cnLevelSpacingF      = 0.025               ; set contour spacing

  plot3(6) = gsn_csm_lat_time(wks3,ts1_lat(lat|:,time|:,lev|3), res)  ; create plot
  res@tiMainString         = lev_label(4)
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = 0.005               ; set min contour level
  res@cnMaxLevelValF       = 0.1               ; set max contour level
  res@cnLevelSpacingF      = 0.005               ; set contour spacing

  plot3(8) = gsn_csm_lat_time(wks3,ts1_lat(lat|:,time|:,lev|4), res)  ; create plot
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = 0.0025               ; set min contour level
  res@cnMaxLevelValF       = 0.05               ; set max contour level
  res@cnLevelSpacingF      = 0.0025               ; set contour spacing
  plot3(10) = gsn_csm_lat_time(wks3,ts1_lat(lat|:,time|:,lev|5), res)  ; create plot
  delete(res21)
  res21=True
  res21@gsnDraw                  = False            ; don't draw yet
  res21@gsnFrame                 = False            ; don't advance frame yet
  res21@vpHeightF= 0.4                    ; change aspect ratio of plot
  res21@vpWidthF = 0.3
  res21@trYMinF      = -90.
  res21@trYMaxF      = 90.
  res21@trXMinF      = 0.
  res21@trXMaxF      = 10.
  res21@tmYLMode        = "Explicit"
  res21@tmYLValues      = (/-60,-30,0,30,60/)
  res21@tmYLLabels      = (/"60S","30S","Eq","30N","60N"/)
  res21@tmYLLabelFontHeightF=0.015
  res21@tmXBMode        = "Explicit"
  res21@tmXBValues      = (/0.,2.,4.,6.,8.,10./)
  res21@tmXBLabels      = (/"0.","2.","4.","6.","8.","10."/)
  res21@tmXBLabelFontHeightF=0.015

  res21@xyLineColor      = "NavyBlue"
  res21@xyLineThicknessF  =3.5
  res21@tiXAxisString = "O3 [ppm]"
  res21@tiMainString = lev_label(0)
  plot3(1) = gsn_csm_xy(wks3,ts1_lat_avg(lat|:,lev|0),ts1_lat_avg&lat, res21)  ; create plot
  res21@tiMainString = lev_label(1)
  plot3(3) = gsn_csm_xy(wks3,ts1_lat_avg(lat|:,lev|1),ts1_lat_avg&lat, res21)
  res21@tiMainString = lev_label(2) 
  res21@tmXBValues      = (/0.,1.,2.,3.,4.,5./)
  res21@tmXBLabels      = (/"0.","1.","2.","3.","4.","5."/)
  res21@trXMinF      = 0.
  res21@trXMaxF      = 5.
  plot3(5) = gsn_csm_xy(wks3,ts1_lat_avg(lat|:,lev|2),ts1_lat_avg&lat, res21)
  res21@tiMainString = lev_label(3)
  res21@trXMinF      = 0.
  res21@trXMaxF      = 0.5
  res21@tmXBValues      = (/0.,0.1,0.2,0.3,0.4,0.5/)
  res21@tmXBLabels      = (/"0.","0.1","0.2","0.3","0.4","0.5"/)
  plot3(7) = gsn_csm_xy(wks3,ts1_lat_avg(lat|:,lev|3),ts1_lat_avg&lat, res21)
  res21@tiMainString = lev_label(4)
  res21@trXMinF      = 0.
  res21@trXMaxF      = 0.1
  res21@tmXBValues      = (/0.,0.02,0.04,0.06,0.08,0.1/)
  res21@tmXBLabels      = (/"0.","0.02","0.04","0.06","0.08","0.1"/)
  plot3(9) = gsn_csm_xy(wks3,ts1_lat_avg(lat|:,lev|4),ts1_lat_avg&lat, res21)
  res21@tiMainString = lev_label(5) 
  res21@trXMaxF      = 0.1
  res21@tmXBValues      = (/0.,0.02,0.04,0.06,0.08,0.1/)
  res21@tmXBLabels      = (/"0.","0.02","0.04","0.06","0.08","0.1"/)
  plot3(11) = gsn_csm_xy(wks3,ts1_lat_avg(lat|:,lev|5),ts1_lat_avg&lat, res21)
  resP=True
  resP@gsnPanelLabelBar    = False
  gsn_panel(wks3,plot3,(/3,4/),resP)

;************************************************
  delete(res)
  wks4 = gsn_open_wks("png","lat_lev"+"CDS-ERA5")                ; send graphics to PNG file
; Set resources necessary to customize Y axis labels
  media=dim_avg_Wrap(dim_avg_n_Wrap(vmro3,(0)))

  gsn_define_colormap(wks4,"MPL_RdYlGn")
  gsn_reverse_colormap(wks4)
  res = True                                    ; plot mods desired
  res@cnLineLabelsOn   =  False
  res@cnLinesOn        =  False
  res@cnFillOn         = True               ; color Fill
  res@cnFillMode       = "RasterFill"       ; Raster Mode
  res@cnLinesOn        =  False             ; Turn off contour lines
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = (/1e-3,2*1e-3,3*1e-3,5*1e-3,7*1e-3,1e-2,2*1e-2,3*1e-2,5*1e-2,7*1e-2,0.1,0.3,0.5,0.7,1.,1.5,2.,3.,4.,5.,6.,7.,8.,9.,10./)
  fmt = "%3.2f"
  res@vpHeightF= 0.5                    ; change aspect ratio of plot
  res@vpWidthF = 0.5
  res@gsnRightString     = "[ppmv]"
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/-60,-30,0,30,60/)
  res@tmXBLabels      = (/"60S","30S","Eq","30N","60N"/)
  res@tmXBLabelFontHeightF=0.015
  res@tmXBLabelFontHeightF=0.015
  res@cnInfoLabelOn=False
 
  p5=gsn_csm_pres_hgt(wks4,media(plev|:,lat|:), res)  ; create plot

;****************************************************
;************************************************
  delete(res)
  plot6=new(12,"graphic")
  month_label=(/"Jannuary","February","March","April","May","June","July","August","September","October","November","December"/)
  wks5 = gsn_open_wks("png","seasonal_lat_lev"+"CDS-OZONE-L3-LP")                ; send graphics to PNG file
; Set resources necessary to customize Y axis labels
zon_vmro3=dim_avg_Wrap(vmro3)
  media_seas=clmMonTLL(zon_vmro3)

  gsn_define_colormap(wks5,"MPL_RdYlGn")
  gsn_reverse_colormap(wks5)
  do im=0,11
     res = True                                    ; plot mods desired
     res@gsnDraw                  = False            ; don't draw yet
     res@gsnFrame                 = False            ; don't advance frame yet
     res@lbLabelBarOn = False
     res@cnLineLabelsOn   =  False
     res@cnLinesOn        =  False
     res@cnFillOn         = True               ; color Fill
     res@cnFillMode       = "RasterFill"       ; Raster Mode
     res@cnLinesOn        =  False             ; Turn off contour lines
     res@cnLevelSelectionMode = "ExplicitLevels"
     res@cnLevels = (/1e-3,2*1e-3,3*1e-3,5*1e-3,7*1e-3,1e-2,2*1e-2,3*1e-2,5*1e-2,7*1e-2,0.1,0.3,0.5,0.7,1.,1.5,2.,3.,4.,5.,6.,7.,8.,9.,10./)
     res@gsnRightString     = "[ppmv]"
     fmt = "%3.2f"
     res@vpHeightF= 0.5                    ; change aspect ratio of plot
     res@vpWidthF = 0.5
     res@tmXBMode        = "Explicit"
     res@tmXBValues      = (/-60,-30,0,30,60/)
     res@tmXBLabels      = (/"60S","30S","Eq","30N","60N"/)
     res@tmXBLabelFontHeightF=0.015
     res@tmXBLabelFontHeightF=0.015
     res@cnInfoLabelOn=False
     res@tiMainString         = month_label(im)
     plot6(im)=gsn_csm_pres_hgt(wks5,media_seas(plev|:,lat|:,month|im), res)  ; create plot
  end do
  resP=True
  resP@gsnPanelLabelBar    = True
  resP@gsnPanelXWhiteSpacePercent =0.5
  gsn_panel(wks5,plot6,(/4,3/),resP)
;*******************************trend***********************************
  dim=dimsizes(vmro3)
 
  opt = False

  vmro32=rmMonAnnCycTLL(zon_vmro3)   ;remove annual cycle
  trend1  = trend_manken(vmro32, opt, 0)  ; ===> px(2,nlat,mlon)
  trend1(1,:,:)=(/trend1(1,:,:)*12*10/)
  trend1!0= "prob_trend" 
  trend1!1="plev"
  trend1!2="lat"
  trend1&plev=vmro32&plev
  trend1&lat=vmro32&lat
  wkst = gsn_open_wks("png","trend_1_"+"CDS-ERA5")         ; send graphics to PNG file
  gsn_define_colormap(wkst,"NEO_div_vegetation_c")
  plott=new(2,"graphic")
  delete(res)
  res = True                                    ; plot mods desired
  res@gsnDraw                  = False            ; don't draw yet
  res@gsnFrame                 = False            ; don't advance frame yet

  res@cnFillOn         = True               ; color Fill
  res@cnFillMode       = "RasterFill"       ; Raster Mode
  res@cnLinesOn        =  False             ; Turn off contour lines
  res@cnLineLabelsOn       =False
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = -0.5               ; set min contour level
  res@cnMaxLevelValF       = 0.5              ; set max contour level
  res@cnLevelSpacingF      = 0.02               ; set contour spacing
  res@tmYRMode             = "Automatic"
  res@tiMainString         = "Trend 1985-2017"    ; plot title
  fmt = "%3.2f"
  res@gsnLeftString      = "[ppmv/decade]"
  res@lbOrientation = "vertical"
  plott(0) = gsn_csm_pres_hgt(wkst,trend1(prob_trend|1,plev|:,lat|:), res)  ; create plot
  res@cnMinLevelValF       = 0               ; set min contour level
  res@cnMaxLevelValF       = 1               ; set max contour level
  res@cnLevelSpacingF      = 0.05               ; set contour spacing
  res@tiMainString         = "Trend Significance"
  res@gsnLeftString      = "[1]"
  plott(1) = gsn_csm_pres_hgt(wkst,trend1(prob_trend|0,plev|:,lat|:), res)
  delete(resP)
  resP=True
  gsn_panel(wkst,plott,(/2,1/),resP)
;*******************************trend 2***********************************
  red_vmro3=zon_vmro3(252:,:,:)
  delete(dim)  
  dim=dimsizes(red_vmro3)

  opt = False

  red_vmro32=rmMonAnnCycTLL(red_vmro3)   ;remove annual cycle
  red_trend1  = trend_manken(red_vmro32, opt, 0)  ; ===> px(2,nlat,mlon)
  red_trend1(1,:,:)=(/red_trend1(1,:,:)*12*10/)
  copy_VarCoords(red_vmro32(0,:,:),red_trend1(0,:,:))
  red_trend1!0= "prob_trend"               ; ===> size 2
  wkst1 = gsn_open_wks("png","trend_2000_2019_"+"CDS-ERA5")         ; send graphics to PNG file
  gsn_define_colormap(wkst1,"NEO_div_vegetation_c")
  plott_red=new(2,"graphic")
  delete(res)
  res = True                                    ; plot mods desired
  res@gsnDraw                  = False            ; don't draw yet
  res@gsnFrame                 = False            ; don't advance frame yet

  res@cnFillOn         = True               ; color Fill
  res@cnFillMode       = "RasterFill"       ; Raster Mode
  res@cnLinesOn        =  False             ; Turn off contour lines
  res@cnLineLabelsOn       =False
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = -0.5               ; set min contour level
  res@cnMaxLevelValF       = 0.5              ; set max contour level
  res@cnLevelSpacingF      = 0.02               ; set contour spacing
  res@tmYRMode             = "Automatic"
  res@tiMainString         = "Trend 2000-2019 "    ; plot title
  fmt = "%3.2f"
  res@gsnLeftString      = "[ppmv/decade]"
  res@lbOrientation = "vertical"
  plott_red(0) = gsn_csm_pres_hgt(wkst1,red_trend1(prob_trend|1,plev|:,lat|:), res)  ; create plot
  res@cnMinLevelValF       = 0               ; set min contour level
  res@cnMaxLevelValF       = 1               ; set max contour level
  res@cnLevelSpacingF      = 0.05               ; set contour spacing
  res@tiMainString         = "Trend Significance"
  res@gsnLeftString      = "[1]"
  plott_red(1) = gsn_csm_pres_hgt(wkst1,red_trend1(prob_trend|0,plev|:,lat|:), res)
  delete(resP)
  resP=True
  gsn_panel(wkst1,plott_red,(/2,1/),resP)


end
