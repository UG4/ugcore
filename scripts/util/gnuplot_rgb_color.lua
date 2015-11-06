-- Copyright (c) 2014:  G-CSC, Goethe University Frankfurt
-- Author: Andreas Vogel
-- 
-- This file is part of UG4.
-- 
-- UG4 is free software: you can redistribute it and/or modify it under the
-- terms of the GNU Lesser General Public License version 3 (as published by the
-- Free Software Foundation) with the following additional attribution
-- requirements (according to LGPL/GPL v3 §7):
-- 
-- (1) The following notice must be displayed in the Appropriate Legal Notices
-- of covered and combined works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (3) The following bibliography is recommended for citation and must be
-- preserved in all covered files:
-- "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
--   parallel geometric multigrid solver on hierarchically distributed grids.
--   Computing and visualization in science 16, 4 (2013), 151-164"
-- "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
--   flexible software system for simulating pde based models on high performance
--   computers. Computing and visualization in science 16, 4 (2013), 165-179"
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU Lesser General Public License for more details.


gnuplot = gnuplot or {}


--[[! --------------------------------------------------------------------------
hey color codes for color names as defined in gnuplot
-------------------------------------------------------------------------------]]

gnuplot.colorname = gnuplot.colorname or 
{
["white"]              ="#ffffff",
["black"]              ="#000000",
["dark-grey"]          ="#a0a0a0",
["red"]                ="#ff0000",
["web-green"]          ="#00c000",
["web-blue"]           ="#0080ff",
["dark-magenta"]       ="#c000ff",
["dark-cyan"]          ="#00eeee",
["dark-orange"]        ="#c04000",
["dark-yellow"]        ="#c8c800",
["royalblue"]          ="#4169e1",
["goldenrod"]          ="#ffc020",
["dark-spring-green"]  ="#008040",
["purple"]             ="#c080ff",
["steelblue"]          ="#306080",
["dark-red"]           ="#8b0000",
["dark-chartreuse"]    ="#408000",
["orchid"]             ="#ff80ff",
["aquamarine"]         ="#7fffd4",
["brown"]              ="#a52a2a",
["yellow"]             ="#ffff00",
["turquoise"]          ="#40e0d0",
["grey0"]              ="#000000",
["grey10"]             ="#1a1a1a",
["grey20"]             ="#333333",
["grey30"]             ="#4d4d4d",
["grey40"]             ="#666666",
["grey50"]             ="#7f7f7f",
["grey60"]             ="#999999",
["grey70"]             ="#b3b3b3",
["grey"]               ="#c0c0c0",
["grey80"]             ="#cccccc",
["grey90"]             ="#e5e5e5",
["grey100"]            ="#ffffff",
["light-red"]          ="#f03232",
["light-green"]        ="#90ee90",
["light-blue"]         ="#add8e6",
["light-magenta"]      ="#f055f0",
["light-cyan"]         ="#e0ffff",
["light-goldenrod"]    ="#eedd82",
["light-pink"]         ="#ffb6c1",
["light-turquoise"]    ="#afeeee",
["gold"]               ="#ffd700",
["green"]              ="#00ff00",
["dark-green"]         ="#006400",
["spring-green"]       ="#00ff7f",
["forest-green"]       ="#228b22",
["sea-green"]          ="#2e8b57",
["blue"]               ="#0000ff",
["dark-blue"]          ="#00008b",
["midnight-blue"]      ="#191970",
["navy"]               ="#000080",
["medium-blue"]        ="#0000cd",
["skyblue"]            ="#87ceeb",
["cyan"]               ="#00ffff",
["magenta"]            ="#ff00ff",
["dark-turquoise"]     ="#00ced1",
["dark-pink"]          ="#ff1493",
["coral"]              ="#ff7f50",
["light-coral"]        ="#f08080",
["orange-red"]         ="#ff4500",
["salmon"]             ="#fa8072",
["dark-salmon"]        ="#e9967a",
["khaki"]              ="#f0e68c",
["dark-khaki"]         ="#bdb76b",
["dark-goldenrod"]     ="#b8860b",
["beige"]              ="#f5f5dc",
["olive"]              ="#a08020",
["orange"]             ="#ffa500",
["violet"]             ="#ee82ee",
["dark-violet"]        ="#9400d3",
["plum"]               ="#dda0dd",
["dark-plum"]          ="#905040",
["dark-olivegreen"]    ="#556b2f",
["orangered4"]         ="#801400",
["brown4"]             ="#801414",
["sienna4"]            ="#804014",
["orchid4"]            ="#804080",
["mediumpurple3"]      ="#8060c0",
["slateblue1"]         ="#8060ff",
["yellow4"]            ="#808000",
["sienna1"]            ="#ff8040",
["tan1"]               ="#ffa040",
["sandybrown"]         ="#ffa060",
["light-salmon"]       ="#ffa070",
["pink"]               ="#ffc0c0",
["khaki1"]             ="#ffff80",
["lemonchiffon"]       ="#ffffc0",
["bisque"]             ="#cdb79e",
["honeydew"]           ="#f0fff0",
["slategrey"]          ="#a0b6cd",
["seagreen"]           ="#c1ffc1",
["antiquewhite"]       ="#cdc0b0",
["chartreuse"]         ="#7cff40",
["greenyellow"]        ="#a0ff20",
["gray"]               ="#bebebe",
["light-gray"]         ="#d3d3d3",
["light-grey"]         ="#d3d3d3",
["dark-gray"]          ="#a0a0a0",
["slategray"]          ="#a0b6cd",
["gray0"]              ="#000000",
["gray10"]             ="#1a1a1a",
["gray20"]             ="#333333",
["gray30"]             ="#4d4d4d",
["gray40"]             ="#666666",
["gray50"]             ="#7f7f7f",
["gray60"]             ="#999999",
["gray70"]             ="#b3b3b3",
["gray80"]             ="#cccccc",
["gray90"]             ="#e5e5e5",
["gray100"]            ="#ffffff",
}


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
function gnuplot.RGBbyLinearHUE(n, S, L, fromH, toH)
	if n < 1 then print("Requirement: n>=2"); exit() end
	if not fromH then fromH = 0 end
	if not toH 	 then toH = 360*(1-1/n) end

	local step = (toH - fromH) / (n-1)
	local colors = {}
	for i=1,n do
		colors[i] = gnuplot.HSLtoRGBhex( (fromH + (i-1)*step + 720)%360, S, L)
	end 
	return colors
end


--[[! --------------------------------------------------------------------------
Returns by lightness (in HSL color space) equally spaced colors 

@param n              number of desired colors
@param H              hue (0-360)
@param S              saturation  (0.0-1.0)
@param fromL		  starting lightnesse value (optional) 
@param toL			  stopping lightnesse value (optional) 
@return               table of RGB color as hex rep. (e.g., "#ff6a8b")
-------------------------------------------------------------------------------]]
function gnuplot.RGBbyLinearLightness(n, H, S, fromL, toL)
	if n < 1 then print("Requirement: n>=2"); exit() end
	if not fromL then fromL = 1/(n+1) end
	if not toL 	 then toL = (1-1/(n+1)) end

	local step = (toL - fromL) / (n-1)
	local colors = {}
	for i=1,n do
		colors[i] = gnuplot.HSLtoRGBhex( H, S, fromL + (i-1)*step)
	end 
	return colors
end

--[[! --------------------------------------------------------------------------
Returns by heu and lighness (in HSL color space) equally spaced color 

@param n              number of desired colors
@param S              saturation  (0.0-1.0)
@param fromH		  starting hue value (optional) 
@param toH			  stopping hue value (optional) 
@param fromL		  starting lightnesse value (optional) 
@param toL			  stopping lightnesse value (optional) 
@return               table of RGB color as hex rep. (e.g., "#ff6a8b")
-------------------------------------------------------------------------------]]
function gnuplot.RGBbyLinearHUEandLightness(n, S, fromH, toH, fromL, toL)
	if n < 1 then print("Requirement: n>=2"); exit() end
	if not fromH then fromH = 0 end
	if not toH 	 then toH = 360*(1-1/n) end
	if not fromL then fromL = 1/(n+1) end
	if not toL 	 then toL = (1-1/(n+1)) end

	local stepL = (toL - fromL) / (n-1)
	local stepH = (toH - fromH) / (n-1)
	local colors = {}
	for i=1,n do
		colors[i] = gnuplot.HSLtoRGBhex( (fromH + (i-1)*stepH + 720)%360, S, fromL + (i-1)*stepL)
	end 
	return colors
end
