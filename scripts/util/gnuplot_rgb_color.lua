
gnuplot = gnuplot or {}

--[[! --------------------------------------------------------------------------
Converts an HSL triplet to RGB
(see http://homepages.cwi.nl/~steven/css/hsl.html or http://easyrgb.com)

@param H              hue 		 (0-360)
@param S              saturation (0.0-1.0)
@param L              lightness  (0.0-1.0)
@return               RGB color, R,G,B in (0-255, 0-255, 0-255)
-------------------------------------------------------------------------------]]
function gnuplot.HSLtoRGB(H, S, L)
   H = H/360
   local g
   if L<=0.5 then g = L*(S+1)
   else 	      g = L+S-L*S
   end
   local r = L*2-g

   local function HUEtoRGB(r, g, H)

      if H<0 then H = H+1 end
      if H>1 then H = H-1 end      
      if 	 H*6<1 then	return r+(g-r)*H*6
      elseif H*2<1 then return g 
      elseif H*3<2 then return r+(g-r)*(2/3-H)*6
      else		        return r
      end
   end

   return HUEtoRGB(r, g, H+1/3)*255, HUEtoRGB(r, g, H)*255, HUEtoRGB(r, g, H-1/3)*255
end

--[[! --------------------------------------------------------------------------
Converts an HSL triplet to RGB
(see http://homepages.cwi.nl/~steven/css/hsl.html or http://easyrgb.com)

@param H              hue 		 (0-360)
@param S              saturation (0.0-1.0)
@param L              lightness  (0.0-1.0)
@return               RGB color as hex rep. (e.g., "#ff6a8b")
-------------------------------------------------------------------------------]]
function gnuplot.HSLtoRGBhex(H, S, L)
	local R,G,B = gnuplot.HSLtoRGB(H,S,L)
	local s = "#"
	s = s..string.format("%02x", R)
	s = s..string.format("%02x", G)
	s = s..string.format("%02x", B)
	return s
end

--[[! --------------------------------------------------------------------------
Converts an RGB triplet to HSL
(see http://easyrgb.com)

@param r              red   (0-255) 	or r = "#ff6a8b" (hex representation)
@param g              green (0-255) 	or nil
@param b              blue  (0-255) 	or nil
@return               HSL color, H,S,L in (0-360, 0-1, 0-1)
-------------------------------------------------------------------------------]]
function gnuplot.RGBtoHSL(r, g, b)
   if type(r) == "string" and g == nil and b == nil then
   		if r:sub(1,1)=="#" and r:len() == 7 then
     		r, g, b = tonumber(r:sub(2,3), 16), tonumber(r:sub(4,5), 16), tonumber(r:sub(6,7), 16)
     	else
     		print(r.." is not a rgb hey representation"); exit()
     	end
   elseif type(r) ~= "number" or type(g) ~= "number" or type(b) ~= "number" then
     	print("Pass rgb hey representation or r,g,b values"); exit()   
   end
   
   if r < 0 or r > 255 or g < 0 or g > 255  or b < 0 or b > 255 then
   	   print("R,G,B values are expected within 0-255"); exit()
   end
   
   r, g, b = r/255, g/255, b/255
   local min = math.min(r, g, b)
   local max = math.max(r, g, b)
   local delta = max - min

   local H, S, L = 0, 0, ((min+max)/2)

   if L > 0 and L < 0.5 then S = delta/(max+min) end
   if L >= 0.5 and L < 1 then S = delta/(2-max-min) end

   if delta > 0 then
      if max == r and max ~= g then H = H + (g-b)/delta end
      if max == g and max ~= b then H = H + 2 + (b-r)/delta end
      if max == b and max ~= r then H = H + 4 + (r-g)/delta end
      H = H / 6;
   end

   if H < 0 then H = H + 1 end
   if H > 1 then H = H - 1 end

   return H * 360, S, L
end

--[[! --------------------------------------------------------------------------
Returns by heu (in HSL color space) equally spaced colors 

@param n              number of desired colors
@param S              saturation (0.0-1.0)
@param L              lightness  (0.0-1.0)
@param fromH		  starting hue value (optional) 
@param toH			  stopping hue value (optional) 
@return               table of RGB color as hex rep. (e.g., "#ff6a8b")
-------------------------------------------------------------------------------]]
function gnuplot.RGBbyEquiHUE(n, S, L, fromH, toH)
	if n < 1 then print("Requirement: n>=2"); exit() end
	if not fromH then fromH = 0 end
	if not toH 	 then toH = 360*(1-1/n) end
	if toH < fromH then print("Requirement:  toH > fromH"); exit() end

	local step = (toH - fromH) / (n-1)
	local colors = {}
	for i=1,n do
		colors[i] = gnuplot.HSLtoRGBhex( (fromH + (i-1)*step)%360, S, L)
	end 
	return colors
end
