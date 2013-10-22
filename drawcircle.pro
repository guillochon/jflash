;+-------------------------------------------------------------
; NAME:
;  	drawcircle
;
; PURPOSE:
;  	draw a circle centered on [x0, y0]. 
;  	YOU HAVE TO HAVE YOUR AXES ESTABLISHED BEFORE CALLING THIS.
;  	Return out arrays for later drawing on a surface, too. These are in two
;	parts, so you won't have a line down the middle when drawing.
;
; CATEGORY:
;       2-D Plotting,
;
; CALLING SEQUENCE:
;	IDL> (some plot command to estable coordinates)
;	IDL> drawcircle, x0, y0, radius
;
;  INPUTS:
;	x0,y0 - origin of circle
;	radius - radius of circle
;
;  KEYWORDS:
;    Inputs:
;	color - 
;	dashed - if set, will draw curve with dashed line
;	angle1 - angle in radians to start drawing arc
;	angle2 - angle in radians to stop drawing arc
;	deg1 - angle in degrees to start drawing arc
;	deg2 - angle in degrees to stop drawing arc
;	bottom - draw arc in bottom of circle
;    Outputs:
;	outX1 - X1 output for other ploting
;	outY1 - Y1  "
;	outX2 - X2  "
;	outY2 - Y2  "
;	npts - # of points to use for drawing circle (all the way around)
;		(default may be too coarse, so curve not smooth)
;  EXAMPLE:
;	IDL> x = findgen(100)  & y = x
;	IDL> Right_Aspect, x, y, IN_POS=[.1,.1,.9,.9], POSITION=position
;	IDL> plot, findgen(100), /nodata, POSITION=position
;	IDL> drawcircle, 40, 60, 20
;	IDL> drawcircle, 40, 60, 15, /dashed
;	IDL> drawcircle, 40, 60, 15, deg1=-60,deg2=-30, /bot
;  LIMITATIONS:
;	angle1 & angle2 can't cross over 0 or !PI
;	ANGLES NEED FIXING (determine where starting from)
;
;  04-Feb-2010 - added angle1 & angle2 keywords (in radians), bottom and npts
;  WRITTEN by Bill Davis
;--------------------------------------------------------------------
pro drawcircle, x0, y0, radius, color=color, dashed=dashed, $
                outX1=outX1, outY1=outY1,  outX2=outX2, outY2=outY2, $
		npts=npts, bottom=bottom, $
		angle1=angle1, angle2=angle2, $
		deg1=deg1, deg2=deg2, $
		noDraw=noDraw, noPlot=noPlot, _extra=_extra

if n_elements( x0 ) eq 0 then x0 = 0
if n_elements( y0 ) eq 0 then y0 = 0
if n_elements( radius ) eq 0 then radius = 5
if keyword_set( dashed ) then line=4 else line=0

if n_elements( deg1 ) ne 0 then angle1 = deg1*!PI/180
if n_elements( deg2 ) ne 0 then angle2 = deg2*!PI/180

if n_elements( angle1 ) ne 0 then angle1 = angle1 MOD (!PI/2)
if n_elements( angle2 ) ne 0 then angle2 = angle2 MOD (!PI/2)

if n_elements( npts ) eq 0 then npts = NINT(!PI*radius) > 100

x = (findgen(npts)+1)/npts*radius
x = [-1*reverse(x), 0, x]
nx = n_elements(x)

y = sqrt( radius^2 - x^2)

xtop=x
ytop=y
if keyword_set( bottom ) eq 0 then begin
   if n_elements( angle1 ) gt 0 OR n_elements( angle2 ) gt 0 then begin
      thetas = atan(x,y)
      if n_elements( angle1 ) eq 0 then angle1 = -1*!PI
      if n_elements( angle2 ) eq 0 then angle2 = !PI
      inds = where( thetas ge angle1 AND thetas LE angle2, nfound )
      if nfound gt 0 then begin
	 ;;;thetas = thetas[inds]
	 xtop = x[inds]
	 ytop = y[inds]
      endif
   endif
endif

outX1 = xtop+x0
outY1 = ytop+y0

if keyword_set( noDraw ) or keyword_set( noPlot ) then return

   ; draw each half separately, so don't get connecting line down the middle
if keyword_set( bottom ) eq 0 then begin
   oplot, xtop+x0, ytop+y0, color=color, line=line, _extra=_extra
   if n_elements( angle1 ) gt 0 OR n_elements( angle2 ) gt 0 then return
endif

   ; don't understand why this is necessary, or if it should be applied in other cases
if n_elements( angle1 ) gt 0 AND n_elements( angle2 ) gt 0 then begin
   if angle1 lt 0 and angle2 lt 0 then begin
      angle1=-1*(!PI/2+angle1)
      angle2=-1*(!PI/2+angle2)
   endif
endif

xbot=x
ybot=y
if keyword_set( bottom ) then begin
   if n_elements( angle1 ) gt 0 OR n_elements( angle2 ) gt 0 then begin
      thetas = atan(x,y)
      if n_elements( angle1 ) eq 0 then angle1 = -1*!PI
      if n_elements( angle2 ) eq 0 then angle2 = !PI
      if angle1 lt angle2 then begin
         inds = where( thetas ge angle1 AND thetas LE angle2, nfound )
      endif else begin
        inds = where( thetas ge angle2 AND thetas LE angle1, nfound )
      endelse
      if nfound gt 0 then begin
	 ;;;thetas = thetas[inds]
	 xbot = x[inds]
	 ybot = y[inds]
      endif
   endif
endif

outX2 = xbot+x0
outY2 = (-1*ybot)+y0
oplot, xbot+x0, (-1*ybot)+y0, color=color, line=line, _extra=_extra

;;;oplot,x,y, color=color
;;;oplot, [x0,x0], [y0-y[npts],y0+y[npts]], color=color

end
