-- Copyright 2016 Stanford University
--
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
--
--     http://www.apache.org/licenses/LICENSE-2.0
--
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

import "regent"

local c     = regentlib.c
local std   = terralib.includec("stdlib.h")
local cmath = terralib.includec("math.h")

-- Some math definitions

local min = regentlib.fmin
local max = regentlib.fmax
local pi  = 2.0*cmath.acos(0.0)

-- Quadrature file name

local quad_file = "radiation_solver/S8.dat"

-- Grid size (x cells, y cells)

local Nx = 1000
local Ny = 1000

--todo: Read from file in Lua

terra get_number_angles()
	var filename : rawstring = "radiation_solver/S8.dat"
  	var f = c.fopen(filename, "rb")
  	var N   : int64[1]
  	c.fscanf(f, "%d\n", N)
  	c.fclose(f)
  	return N[0]
end

local N_angles = get_number_angles()

-- Domain size

local Lx = 1.0
local Ly = 1.0

-- Grid spacing

local dx = Lx/Nx
local dy = Ly/Ny

-- Albedo

local omega = .7

-- Wall emissivity

local emiss_east  = 1
local emiss_west  = 1
local emiss_south = 1
local emiss_north = 1

-- Wall temperatures

local SB = 5.67e-8

local T_west  = 2725
local T_east  = 300
local T_south = 300
local T_north = 300

-- Procedure parameters

local tol   = 1e-6   -- solution tolerance
local res   = 1      -- initial residual
local gamma = 0.5    -- 1 for step differencing, 0.5 for diamond differencing


terra read_val(f : &c.FILE, val : &double)
  return c.fscanf(f, "%lf\n", &val[0])
end

-- Fieldspaces

-- Internal cell values are essentially private,
-- face values are what need to be passed to downstream neighbor
-- Update cell value, then update downstream face values

-- quadrature information

fspace angle_value {
	xi : double,
	eta : double,
	w : double
}

fspace point {

	-- Intensities per angle quadrant
	I_1 : double[N_angles/4],		-- cell center intensity per angle
	I_2 : double[N_angles/4],		-- cell center intensity per angle
	I_3 : double[N_angles/4],		-- cell center intensity per angle
	I_4 : double[N_angles/4],		-- cell center intensity per angle

	Iiter_1 : double[N_angles/4],   -- iterative intensity per angle
	Iiter_2 : double[N_angles/4],   -- iterative intensity per angle
	Iiter_3 : double[N_angles/4],   -- iterative intensity per angle
	Iiter_4 : double[N_angles/4],   -- iterative intensity per angle

	G : double,						-- intensity summation over all angles

	-- source term
	S : double,  	

	-- constants
	T : double,		-- blackbody source
	Ib : double,
	sigma : double,	-- Extinction coefficient

	-- x,y coordinates for visualization
	x : double,
	y : double
}

fspace x_face {
	Ifx_1 : double[N_angles/4],		-- x face intensity per angle
	Ifx_2 : double[N_angles/4],		-- x face intensity per angle
	Ifx_3 : double[N_angles/4],		-- x face intensity per angle
	Ifx_4 : double[N_angles/4]		-- x face intensity per angle
}

fspace y_face {
	Ify_1 : double[N_angles/4], 	-- y face intensity per angle
	Ify_2 : double[N_angles/4], 	-- y face intensity per angle
	Ify_3 : double[N_angles/4], 	-- y face intensity per angle
	Ify_4 : double[N_angles/4] 		-- y face intensity per angle
}

task initialize_x_faces(x_faces : region(ispace(int2d), x_face))
where reads writes(x_faces.Ifx_1, x_faces.Ifx_2,
				x_faces.Ifx_3, x_faces.Ifx_4)
do

	for i in x_faces do
		for m = 0, N_angles/4 do
			x_faces[i].Ifx_1[m] = 0.0
			x_faces[i].Ifx_2[m] = 0.0
			x_faces[i].Ifx_3[m] = 0.0
			x_faces[i].Ifx_4[m] = 0.0
		end
	end

end

task initialize_y_faces(y_faces : region(ispace(int2d), y_face))
where reads writes(y_faces.Ify_1, y_faces.Ify_2,
				y_faces.Ify_3, y_faces.Ify_4)
do

	for i in y_faces do
		for m = 0, N_angles/4 do
			y_faces[i].Ify_1[m] = 0.0
			y_faces[i].Ify_2[m] = 0.0
			y_faces[i].Ify_3[m] = 0.0
			y_faces[i].Ify_4[m] = 0.0
		end
	end

end

task initialize_angle_values(angle_values : region(ispace(int1d), angle_value))
where writes (angle_values.xi, angle_values.eta, angle_values.w)

do

	-- Now set the angle_value information, read in from file.

  	var val : double[1]

  	var f = c.fopen("radiation_solver/S8.dat", "rb")

  	read_val(f, val) -- gets rid of num angles

  	for i in angle_values do
  		read_val(f, val)
  		var index : int = (i % 4) * (N_angles/4) + (i / 4)
		angle_values[index].xi = val[0]
	end

	for i in angle_values do
  		read_val(f, val)
  		var index : int = (i % 4) * (N_angles/4) + (i / 4)
		angle_values[index].eta = val[0]
	end

	for i in angle_values do
  		read_val(f, val)
  		var index : int = (i % 4) * (N_angles/4) + (i / 4)
		angle_values[index].w = val[0]
	end
  

 	c.fclose(f)

end

task initialize(points : region(ispace(int2d), point))
where
  reads writes(points.T,
  	points.I_1,
  	points.I_2,
  	points.I_3,
  	points.I_4,
  	points.Iiter_1,
  	points.Iiter_2,
  	points.Iiter_3,
  	points.Iiter_4), 
  writes(points.Ib, points.sigma, points.G, points.S)
do

  	-- First loop over all points to set the constant values.

  	for i in points do

	    -- Blackbody source

	    points[i].T  = 300.0
	    points[i].Ib = (SB/pi)*cmath.pow(points[i].T,4.0)

	    -- Extinction coefficient

	    points[i].sigma = 3.0

	    -- Intensities (cell-based)

	    for m = 0, N_angles/4 do
			points[i].I_1[m]     = 0.0
	    	points[i].Iiter_1[m] = 0.0

	    	points[i].I_2[m]     = 0.0
	    	points[i].Iiter_2[m] = 0.0

	    	points[i].I_3[m]     = 0.0
	    	points[i].Iiter_3[m] = 0.0

	    	points[i].I_4[m]     = 0.0
	    	points[i].Iiter_4[m] = 0.0
		end
	    
    	points[i].G = 0.0

    	-- Source term

    	points[i].S = 0.0

  	end

end

-- Update source term
task source_term(points : region(ispace(int2d), point),
				angle_values : region(ispace(int1d), angle_value))
where
  reads (points.Iiter_1, points.Iiter_2, points.Iiter_3, points.Iiter_4,
  	points.Ib, points.sigma, angle_values.w),
  reads writes (points.S)
do

  	-- Get array bounds

  	var limits = points.bounds

  	-- Loop over all angles and grid cells to compute the source term
  	-- for the current iteration.

  	for i = limits.lo.x, limits.hi.x + 1 do
    	for j = limits.lo.y, limits.hi.y + 1 do
      		points[{i,j}].S = (1.0-omega)*SB*points[{i,j}].sigma*points[{i,j}].Ib
     	 	for m = 0, N_angles/4 do
        		points[{i,j}].S = points[{i,j}].S + omega*points[{i,j}].sigma/(4.0*pi)*angle_values[m].w*points[{i,j}].Iiter_1[m]
      		end
      		for m = 0, N_angles/4 do
        		points[{i,j}].S = points[{i,j}].S + omega*points[{i,j}].sigma/(4.0*pi)*angle_values[m+N_angles/4].w*points[{i,j}].Iiter_2[m]
      		end
      		for m = 0, N_angles/4 do
        		points[{i,j}].S = points[{i,j}].S + omega*points[{i,j}].sigma/(4.0*pi)*angle_values[m+(N_angles/4)*2].w*points[{i,j}].Iiter_3[m]
      		end
      		for m = 0, N_angles/4 do
        		points[{i,j}].S = points[{i,j}].S + omega*points[{i,j}].sigma/(4.0*pi)*angle_values[m+(N_angles/4)*3].w*points[{i,j}].Iiter_4[m]
      		end
      		-- c.printf("source term %d %d = %lf ", i, j, points[{i,j}].S)
    	end
    	-- c.printf("\n")
  	end

end

-- Partition x,y grid into tiles
task make_interior_partition(points : region(ispace(int2d), point),
                            tiles : ispace(int2d),
                            ntiles : int64, nx : int64, ny : int64)
	var coloring = c.legion_domain_point_coloring_create()
	for tile in tiles do
	    var lo = int2d { x = tile.x * nx / ntiles, y = tile.y * ny / ntiles}
		var hi = int2d { x = (tile.x + 1) * nx / ntiles - 1, y = (tile.y + 1) * ny / ntiles - 1}
		var rect = rect2d {lo = lo, hi = hi}
		c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
	end
	var p = partition(disjoint, points, coloring, tiles)
	c.legion_domain_point_coloring_destroy(coloring)
	return p
end

-- x
task make_interior_partition_x_hi(faces : region(ispace(int2d), x_face),
							tiles : ispace(int2d),
							ntiles : int64, nx : int64, ny : int64) 
	
	var coloring = c.legion_domain_point_coloring_create()
	for tile in tiles do

		var val : int = -1
		if tile.x == tiles.bounds.hi.x-1 then val = 0 end -- include extra face in last partition
	    var lo = int2d { x = tile.x * (nx-1) / ntiles, y = tile.y * ny / ntiles} 
		var hi = int2d { x = (tile.x+1) * (nx-1) / ntiles + val, y = (tile.y+1) * ny / ntiles - 1}

		-- Create an empty partition
		if hi.x >= nx then
			lo.x = 1
			hi.x = 0
		end

		var rect = rect2d {lo = lo, hi = hi}
		c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
	end
	var p = partition(disjoint, faces, coloring, tiles)
	c.legion_domain_point_coloring_destroy(coloring)
	return p
end

task make_interior_partition_x_lo(faces : region(ispace(int2d), x_face),
							tiles : ispace(int2d),
							ntiles : int64, nx : int64, ny : int64) 
	
	var coloring = c.legion_domain_point_coloring_create()
	for tile in tiles do

		-- c.printf("tile.x=%d tile.y=%d \n", tile.x, tile.y)

		var val : int = 1
		if tile.x == 1 then val = 0 end -- include extra face in first partition
	    var lo = int2d { x = (tile.x-1) * (nx-1) / ntiles + val, y = tile.y * ny / ntiles} 
		var hi = int2d { x = tile.x * (nx-1) / ntiles, y = (tile.y+1) * ny / ntiles - 1}

		-- Create an empty partition
		if lo.x < 0 then
			lo.x = 1
			hi.x = 0
		end

		-- c.printf("x lo bounds x=%d-%d y=%d-%d \n", lo.x, hi.x, lo.y, hi.y)

		var rect = rect2d {lo = lo, hi = hi}
		c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
	end
	var p = partition(disjoint, faces, coloring, tiles)
	c.legion_domain_point_coloring_destroy(coloring)
	return p
end


task make_interior_partition_y_hi(faces : region(ispace(int2d), y_face),
							tiles : ispace(int2d),
							ntiles : int64, nx : int64, ny : int64)
	
	var coloring = c.legion_domain_point_coloring_create()
	for tile in tiles do

		var val : int = -1
		if tile.y == tiles.bounds.hi.y-1 then val = 0 end -- include extra face in last partition
	   	var lo = int2d { x = tile.x * nx / ntiles, y = tile.y * (ny-1) / ntiles}
		var hi = int2d { x = (tile.x+1) * nx / ntiles - 1, y = (tile.y+1) * (ny-1) / ntiles + val}

		-- Create an empty partition

		if hi.y >= ny then
			lo.y = 1
			hi.y = 0
		end

		var rect = rect2d {lo = lo, hi = hi}
		c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
	end
	var p = partition(disjoint, faces, coloring, tiles)
	c.legion_domain_point_coloring_destroy(coloring)
	return p
end

task make_interior_partition_y_lo(faces : region(ispace(int2d), y_face),
							tiles : ispace(int2d),
							ntiles : int64, nx : int64, ny : int64)
	
	var coloring = c.legion_domain_point_coloring_create()
	for tile in tiles do

		var val : int = 1
		if tile.y == 1 then val = 0 end -- include extra face in first partition
	   	var lo = int2d { x = tile.x * nx / ntiles, y = (tile.y-1) * (ny-1) / ntiles + val}
		var hi = int2d { x = (tile.x+1) * nx / ntiles - 1, y = tile.y * (ny-1) / ntiles}

		-- Create an empty partition

		if lo.y < 0 then
			lo.y = 1
			hi.y = 0
		end

		var rect = rect2d {lo = lo, hi = hi}
		c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
	end
	var p = partition(disjoint, faces, coloring, tiles)
	c.legion_domain_point_coloring_destroy(coloring)
	return p
end


task west_bound(faces : region(ispace(int2d), x_face),
				angle_values : region(ispace(int1d), angle_value))
where
	reads (angle_values.w, angle_values.xi),
	reads writes (faces.Ifx_1, faces.Ifx_2, faces.Ifx_3, faces.Ifx_4) 
do

	-- Get array bounds

  	var limits = faces.bounds

  	-- c.printf("west bound x = %d ", limits.lo.x)

 	 -- Temporary variables for the west bound

  	var reflect : double = 0.0
  	var epsw    : double = emiss_west
  	var Tw      : double = T_west

  	for j = limits.lo.y, limits.hi.y + 1 do

  		-- Calculate reflect

  		reflect = 0
  		for m = 0, N_angles do
  			var face_value : double = 0.0
  			if m < N_angles/4 then
  				face_value = faces[{limits.lo.x,j}].Ifx_1[m]
  			elseif m < (N_angles/4)*2 then
  				face_value = faces[{limits.lo.x,j}].Ifx_2[m - (N_angles/4)]
  			elseif m < (N_angles/4)*3 then
  				face_value = faces[{limits.lo.x,j}].Ifx_3[m - (N_angles/4)*2]
  			else
  				face_value = faces[{limits.lo.x,j}].Ifx_4[m - (N_angles/4)*3]
  			end
  			if angle_values[m].xi < 0 then
  				reflect = reflect + (1.0-epsw)/pi*angle_values[m].w*cmath.fabs(angle_values[m].xi)*face_value
  			end
  		end


  		-- Set Ifx values using reflect

  		for m = 0, N_angles/4 do
  			var value : double = epsw*SB*cmath.pow(Tw,4.0)/pi + reflect

  			var angle : int = m
  			if angle_values[angle].xi > 0 then
  				faces[{limits.lo.x,j}].Ifx_1[m] = value
  				-- c.printf("west bound=%lf angle=%d x=%d y=%d\n", faces[{limits.lo.x,j}].Ifx_1[m], angle, limits.lo.x, j)
  			end

  			angle = m + (N_angles/4)
  			if angle_values[angle].xi > 0 then
  				faces[{limits.lo.x,j}].Ifx_2[m] = value
  				-- c.printf("west bound=%lf angle=%d x=%d y=%d\n", faces[{limits.lo.x,j}].Ifx_2[m], angle, limits.lo.x, j)
  			end

  			angle = m + (N_angles/4) * 2
  			if angle_values[angle].xi > 0 then
  				faces[{limits.lo.x,j}].Ifx_3[m] = value
  				-- c.printf("west bound=%lf angle=%d x=%d y=%d\n", faces[{limits.lo.x,j}].Ifx_3[m], angle, limits.lo.x, j)
  			end

  			angle = m + (N_angles/4) * 3
  			if angle_values[angle].xi > 0 then
  				faces[{limits.lo.x,j}].Ifx_4[m] = value
  				-- c.printf("west bound=%lf angle=%d x=%d y=%d\n", faces[{limits.lo.x,j}].Ifx_4[m], angle, limits.lo.x, j)
  			end
  		end

  	end
end

task east_bound(faces : region(ispace(int2d), x_face),
				angle_values : region(ispace(int1d), angle_value))
where
	reads (angle_values.w, angle_values.xi),
	reads writes (faces.Ifx_1, faces.Ifx_2, faces.Ifx_3, faces.Ifx_4) 
do

	-- Get array bounds

  	var limits = faces.bounds

 	 -- Temporary variables for the east bound

  	var reflect : double = 0.0
  	var epsw    : double = emiss_east
  	var Tw      : double = T_east

  	for j = limits.lo.y, limits.hi.y+1 do

  		-- Calculate reflect

  		reflect = 0
  		for m = 0, N_angles do
  			var face_value : double = 0.0
  			if m < N_angles/4 then
  				face_value = faces[{limits.hi.x,j}].Ifx_1[m]
  			elseif m < (N_angles/4)*2 then
  				face_value = faces[{limits.hi.x,j}].Ifx_2[m - (N_angles/4)]
  			elseif m < (N_angles/4)*3 then
  				face_value = faces[{limits.hi.x,j}].Ifx_3[m - (N_angles/4)*2]
  			else
  				face_value = faces[{limits.hi.x,j}].Ifx_4[m - (N_angles/4)*3]
  			end
  			if angle_values[m].xi > 0 then
  				reflect = reflect + (1.0-epsw)/pi*angle_values[m].w*angle_values[m].xi*face_value
  			end
  		end


  		-- Set Ifx values using reflect

  		for m = 0, N_angles/4 do
  			var value : double = epsw*SB*cmath.pow(Tw,4.0)/pi + reflect

  			var angle : int = m
  			if angle_values[angle].xi < 0 then
  				faces[{limits.hi.x,j}].Ifx_1[m] = value
  			end

  			angle = m + (N_angles/4)
  			if angle_values[angle].xi < 0 then
  				faces[{limits.hi.x,j}].Ifx_2[m] = value
  			end

  			angle = m + (N_angles/4) * 2
  			if angle_values[angle].xi < 0 then
  				faces[{limits.hi.x,j}].Ifx_3[m] = value
  			end

  			angle = m + (N_angles/4) * 3
  			if angle_values[angle].xi < 0 then
  				faces[{limits.hi.x,j}].Ifx_4[m] = value
  			end
  		end

  	end
end

task north_bound(faces : region(ispace(int2d), y_face),
				angle_values : region(ispace(int1d), angle_value))
where
	reads (angle_values.w, angle_values.eta),
	reads writes (faces.Ify_1, faces.Ify_2, faces.Ify_3, faces.Ify_4)
do

	-- Get array bounds

  	var limits = faces.bounds

 	 -- Temporary variables for the west bound

  	var reflect : double = 0.0
  	var epsw    : double = emiss_north
  	var Tw      : double = T_north

  	for i = limits.lo.x, limits.hi.x+1 do

  		-- Calculate reflect

  		reflect = 0
  		for m = 0, N_angles do
  			var face_value : double = 0.0
  			if m < N_angles/4 then
  				face_value = faces[{i,limits.hi.y}].Ify_1[m]
  			elseif m < (N_angles/4)*2 then
  				face_value = faces[{i,limits.hi.y}].Ify_2[m - (N_angles/4)]
  			elseif m < (N_angles/4)*3 then
  				face_value = faces[{i,limits.hi.y}].Ify_3[m - (N_angles/4)*2]
  			else
  				face_value = faces[{i,limits.hi.y}].Ify_4[m - (N_angles/4)*3]
  			end
  			if angle_values[m].eta > 0 then
  				reflect = reflect + (1.0-epsw)/pi*angle_values[m].w*angle_values[m].eta*face_value
  			end
  		end


  		-- Set Ify values using reflect

  		for m = 0, N_angles/4 do
  			var value : double = epsw*SB*cmath.pow(Tw,4.0)/pi + reflect

  			var angle : int = m
  			if angle_values[angle].eta < 0 then
  				faces[{i,limits.hi.y}].Ify_1[m] = value
  			end

  			angle = m + (N_angles/4)
  			if angle_values[angle].eta < 0 then
  				faces[{i,limits.hi.y}].Ify_2[m] = value
  			end

  			angle = m + (N_angles/4) * 2
  			if angle_values[angle].eta < 0 then
  				faces[{i,limits.hi.y}].Ify_3[m] = value
  			end

  			angle = m + (N_angles/4) * 3
  			if angle_values[angle].eta < 0 then
  				faces[{i,limits.hi.y}].Ify_4[m] = value
  			end
  		end
  	end
end

task south_bound(faces : region(ispace(int2d), y_face),
				angle_values : region(ispace(int1d), angle_value))
where
	reads (angle_values.w, angle_values.eta),
	reads writes (faces.Ify_1, faces.Ify_2, faces.Ify_3, faces.Ify_4)
do

	-- Get array bounds

  	var limits = faces.bounds

 	 -- Temporary variables for the west bound

  	var reflect : double = 0.0
  	var epsw    : double = emiss_south
  	var Tw      : double = T_south

  	for i = limits.lo.x, limits.hi.x+1 do

  		-- Calculate reflect

  		reflect = 0
  		for m = 0, N_angles do
  			var face_value : double = 0.0
  			if m < N_angles/4 then
  				face_value = faces[{i,limits.lo.y}].Ify_1[m]
  			elseif m < (N_angles/4)*2 then
  				face_value = faces[{i,limits.lo.y}].Ify_2[m - (N_angles/4)]
  			elseif m < (N_angles/4)*3 then
  				face_value = faces[{i,limits.lo.y}].Ify_3[m - (N_angles/4)*2]
  			else
  				face_value = faces[{i,limits.lo.y}].Ify_4[m - (N_angles/4)*3]
  			end
  			if angle_values[m].eta < 0 then
  				reflect = reflect + (1.0-epsw)/pi*angle_values[m].w*cmath.fabs(angle_values[m].eta)*face_value
  			end
  		end


  		-- Set Ify values using reflect

  		for m = 0, N_angles/4 do
  			var value : double = epsw*SB*cmath.pow(Tw,4.0)/pi + reflect

  			var angle : int = m
  			if angle_values[m].eta > 0 then
  				faces[{i,limits.lo.y}].Ify_1[m] = value
  			end

  			angle = m + (N_angles/4)
  			if angle_values[angle].eta > 0 then
  				faces[{i,limits.lo.y}].Ify_2[m] = value
  			end

  			angle = m + (N_angles/4) * 2
  			if angle_values[angle].eta > 0 then
  				faces[{i,limits.lo.y}].Ify_3[m] = value
  			end

  			angle = m + (N_angles/4) * 3
  			if angle_values[angle].eta > 0 then
  				faces[{i,limits.lo.y}].Ify_4[m] = value
  			end
  		end

  	end
end

-- +x, +y
task sweep_1(points : region(ispace(int2d), point),
		   x_faces : region(ispace(int2d), x_face),
		   y_faces : region(ispace(int2d), y_face),
		   ghost_x_faces : region(ispace(int2d), x_face),
		   ghost_y_faces : region(ispace(int2d), y_face),
		   angle_values : region(ispace(int1d), angle_value))
where
  reads (angle_values.xi, angle_values.eta, points.S, points.sigma, ghost_x_faces.Ifx_1, ghost_y_faces.Ify_1),
  reads writes(points.I_1, x_faces.Ifx_1, y_faces.Ify_1)
do

  -- Get array bounds and some temporary index variables for sweeping

  var limits = points.bounds

  var indx   : int64 = 0
  var indy   : int64 = 0

  var dindx  : int64 = 1
  var startx : int64 = limits.lo.x
  var endx   : int64 = limits.hi.x + 1

  var dindy  : int64 = 1
  var starty : int64 = limits.lo.y
  var endy   : int64 = limits.hi.y + 1

  
  -- Outer loop over all angles.
  for m = 0, N_angles/4 do

  	var angle : int64 = m

  	-- c.printf("angle value x = %lf, y = %lf\n", angle_values[angle].xi, angle_values[angle].eta)

    -- Use our direction and increments for the sweep.

    for j = starty,endy,dindy do
      for i = startx,endx,dindx do

      	-- indx and indy are the upwind indices
      	indx = i - min(dindx,0)
        indy = j - min(dindy,0)

        var upwind_x_value : double = 0.0
        if indx < x_faces.bounds.lo.x then
        	var ghost_x_limits = ghost_x_faces.bounds
        	upwind_x_value = ghost_x_faces[{ghost_x_limits.hi.x,j}].Ifx_1[m]
        else
        	upwind_x_value = x_faces[{indx, j}].Ifx_1[m]
       	end

       	var upwind_y_value : double = 0.0
       	if indy < y_faces.bounds.lo.y then
       		var ghost_y_limits = ghost_y_faces.bounds
        	upwind_y_value = ghost_y_faces[{i,ghost_y_limits.hi.y}].Ify_1[m]
       	else
       		upwind_y_value = y_faces[{i,indy}].Ify_1[m]
       	end

        -- Integrate to compute cell-centered value of I.

        points[{i,j}].I_1[m] = (points[{i,j}].S * dx * dy 
        	+ cmath.fabs(angle_values[angle].xi) * dy * upwind_x_value/gamma 
        	+ cmath.fabs(angle_values[angle].eta) * dx * upwind_y_value/gamma)
        	/(points[{i,j}].sigma * dx * dy 
        	+ cmath.fabs(angle_values[angle].xi) * dy/gamma 
        	+ cmath.fabs(angle_values[angle].eta) * dx/gamma)

        -- c.printf("x=%d,y=%d,angle=%d I = %lf \n", i, j, angle, points[{i,j}].I_1[m])

        -- Compute intensities on downwind faces

        x_faces[{indx+dindx, j}].Ifx_1[m] = (points[{i,j}].I_1[m] - (1-gamma)*upwind_x_value)/gamma
        y_faces[{i, indy+dindy}].Ify_1[m] = (points[{i,j}].I_1[m] - (1-gamma)*upwind_y_value)/gamma

      end
    end
  end
end

-- +x, -y
task sweep_2(points : region(ispace(int2d), point),
		   x_faces : region(ispace(int2d), x_face),
		   y_faces : region(ispace(int2d), y_face),
		   ghost_x_faces : region(ispace(int2d), x_face),
		   ghost_y_faces : region(ispace(int2d), y_face),
		   angle_values : region(ispace(int1d), angle_value))
where
  reads (angle_values.xi, angle_values.eta, points.S, points.sigma, ghost_x_faces.Ifx_2, ghost_y_faces.Ify_2),
  reads writes(points.I_2, x_faces.Ifx_2, y_faces.Ify_2)
do

  -- Get array bounds and some temporary index variables for sweeping

  var limits = points.bounds

  var indx   : int64 = 0
  var indy   : int64 = 0

  var dindx  : int64 = 1
  var startx : int64 = limits.lo.x
  var endx   : int64 = limits.hi.x + 1

  var dindy  : int64 = -1
  var starty : int64 = limits.hi.y
  var endy   : int64 = limits.lo.y - 1
  
  -- Outer loop over all angles.
  for m = 0, N_angles/4 do

  	var angle : int64 = m + N_angles/4

  	-- c.printf("angle value x = %lf, y = %lf\n", angle_values[angle].xi, angle_values[angle].eta)

    -- Use our direction and increments for the sweep.

    for j = starty,endy,dindy do
      for i = startx,endx,dindx do

      	indx = i - min(dindx,0)
        indy = j - min(dindy,0)

        var upwind_x_value : double = 0.0
        if indx < x_faces.bounds.lo.x then
        	var ghost_x_limits = ghost_x_faces.bounds
        	upwind_x_value = ghost_x_faces[{ghost_x_limits.hi.x,j}].Ifx_2[m]
        else
        	upwind_x_value = x_faces[{indx,j}].Ifx_2[m]
       	end

       	var upwind_y_value : double = 0.0
       	if indy > y_faces.bounds.hi.y then
       		var ghost_y_limits = ghost_y_faces.bounds
        	upwind_y_value = ghost_y_faces[{i,ghost_y_limits.lo.y}].Ify_2[m]
       	else
       		upwind_y_value = y_faces[{i,indy}].Ify_2[m]
       	end

        -- Integrate to compute cell-centered value of I.

        points[{i,j}].I_2[m] = (points[{i,j}].S * dx * dy 
        	+ cmath.fabs(angle_values[angle].xi) * dy * upwind_x_value/gamma 
        	+ cmath.fabs(angle_values[angle].eta) * dx * upwind_y_value/gamma)
        	/(points[{i,j}].sigma * dx * dy + cmath.fabs(angle_values[angle].xi) * dy/gamma + cmath.fabs(angle_values[angle].eta) * dx/gamma)

        -- c.printf("x=%d,y=%d,angle=%d I = %lf \n", i, j, angle, points[{i,j}].I_2[m])

        x_faces[{indx+dindx, j}].Ifx_2[m] = (points[{i,j}].I_2[m] - (1-gamma)*upwind_x_value)/gamma
        y_faces[{i, indy+dindy}].Ify_2[m] = (points[{i,j}].I_2[m] - (1-gamma)*upwind_y_value)/gamma

      end
    end
  end
end

-- -x, +y
task sweep_3(points : region(ispace(int2d), point),
		   x_faces : region(ispace(int2d), x_face),
		   y_faces : region(ispace(int2d), y_face),
		   ghost_x_faces : region(ispace(int2d), x_face),
		   ghost_y_faces : region(ispace(int2d), y_face),
		   angle_values : region(ispace(int1d), angle_value))
where
  reads (angle_values.xi, angle_values.eta, points.S, points.sigma, ghost_x_faces.Ifx_3, ghost_y_faces.Ify_3),
  reads writes(points.I_3, x_faces.Ifx_3, y_faces.Ify_3)
do

  -- Get array bounds and some temporary index variables for sweeping

  var limits = points.bounds

  var indx   : int64 = 0
  var indy   : int64 = 0

  var dindx  : int64 = -1
  var startx : int64 = limits.hi.x
  var endx   : int64 = limits.lo.x - 1

  var dindy  : int64 = 1
  var starty : int64 = limits.lo.y
  var endy   : int64 = limits.hi.y + 1

  
  -- Outer loop over all angles.
  for m = 0, N_angles/4 do

  	var angle : int64 = m + (N_angles/4)*2

  	-- c.printf("angle value x = %lf, y = %lf\n", angle_values[angle].xi, angle_values[angle].eta)

    -- Use our direction and increments for the sweep.

    for j = starty,endy,dindy do
      for i = startx,endx,dindx do

      	indx = i - min(dindx,0)
        indy = j - min(dindy,0)

        var upwind_x_value : double = 0.0
        if indx > x_faces.bounds.hi.x then
        	var ghost_x_limits = ghost_x_faces.bounds
        	upwind_x_value = ghost_x_faces[{ghost_x_limits.lo.x,j}].Ifx_3[m]
        else
        	upwind_x_value = x_faces[{indx,j}].Ifx_3[m]
       	end

       	var upwind_y_value : double = 0.0
       	if indy < y_faces.bounds.lo.y then
       		var ghost_y_limits = ghost_y_faces.bounds
        	upwind_y_value = ghost_y_faces[{i,ghost_y_limits.hi.y}].Ify_3[m]
       	else
       		upwind_y_value = y_faces[{i,indy}].Ify_3[m]
       	end

        -- Integrate to compute cell-centered value of I.

        points[{i,j}].I_3[m] = (points[{i,j}].S * dx * dy 
        	+ cmath.fabs(angle_values[angle].xi) * dy * upwind_x_value/gamma 
        	+ cmath.fabs(angle_values[angle].eta) * dx * upwind_y_value/gamma)
        	/(points[{i,j}].sigma * dx * dy + cmath.fabs(angle_values[angle].xi) * dy/gamma + cmath.fabs(angle_values[angle].eta) * dx/gamma)

        -- c.printf("x=%d,y=%d,angle=%d I = %lf \n", i, j, angle, points[{i,j}].I_3[m])


        x_faces[{indx+dindx, j}].Ifx_3[m] = (points[{i,j}].I_3[m] - (1-gamma)*upwind_x_value)/gamma
        y_faces[{i, indy+dindy}].Ify_3[m] = (points[{i,j}].I_3[m] - (1-gamma)*upwind_y_value)/gamma

      end
      -- c.printf("\n")
    end
  end
end

-- -x, -y
task sweep_4(points : region(ispace(int2d), point),
		   x_faces : region(ispace(int2d), x_face),
		   y_faces : region(ispace(int2d), y_face),
		   ghost_x_faces : region(ispace(int2d), x_face),
		   ghost_y_faces : region(ispace(int2d), y_face),
		   angle_values : region(ispace(int1d), angle_value))
where
  reads (angle_values.xi, angle_values.eta, points.S, points.sigma, ghost_x_faces.Ifx_4, ghost_y_faces.Ify_4),
  reads writes(points.I_4, x_faces.Ifx_4, y_faces.Ify_4)
do

	-- c.printf("x face bounds x=%d-%d y=%d-%d", x_faces.bounds.lo.x, x_faces.bounds.hi.x, x_faces.bounds.lo.y, x_faces.bounds.hi.y)
	-- c.printf("y face bounds x=%d-%d y=%d-%d", y_faces.bounds.lo.x, y_faces.bounds.hi.x, y_faces.bounds.lo.y, y_faces.bounds.hi.y)
  -- Get array bounds and some temporary index variables for sweeping

  var limits = points.bounds

  var indx   : int64 = 0
  var indy   : int64 = 0

  var dindx  : int64 = -1
  var startx : int64 = limits.hi.x
  var endx   : int64 = limits.lo.x - 1

  var dindy  : int64 = -1
  var starty : int64 = limits.hi.y
  var endy   : int64 = limits.lo.y - 1

  
  -- Outer loop over all angles.
  for m = 0, N_angles/4 do

  	var angle : int64 = m + (N_angles/4) *3

  	-- c.printf("angle value x = %lf, y = %lf\n", angle_values[angle].xi, angle_values[angle].eta)

    -- Use our direction and increments for the sweep.

    for j = starty,endy,dindy do
      for i = startx,endx,dindx do

      	indx = i - min(dindx,0)
        indy = j - min(dindy,0)

        var upwind_x_value : double = 0.0
        if indx > x_faces.bounds.hi.x then
        	var ghost_x_limits = ghost_x_faces.bounds
        	upwind_x_value = ghost_x_faces[{ghost_x_limits.lo.x,j}].Ifx_4[m]
        	-- c.printf("x=%d,y=%d,angle=%d upwind x face ghost = %lf bounds = %d \n", indx, j, angle, upwind_x_value, x_faces.bounds.hi.x)
        else
        	upwind_x_value = x_faces[{indx,j}].Ifx_4[m]
        	-- c.printf("x=%d,y=%d,angle=%d upwind x face = %lf \n", indx, j, angle, upwind_x_value)
       	end

       	var upwind_y_value : double = 0.0
       	if indy > y_faces.bounds.hi.y then
       		var ghost_y_limits = ghost_y_faces.bounds
        	upwind_y_value = ghost_y_faces[{i,ghost_y_limits.lo.y}].Ify_4[m]
        	-- c.printf("x=%d,y=%d,angle=%d upwind y face ghost = %lf \n", i, indy, angle, upwind_y_value)
       	else
       		upwind_y_value = y_faces[{i,indy}].Ify_4[m]
       		-- c.printf("x=%d,y=%d,angle=%d upwind y face = %lf \n", i, indy, angle, upwind_y_value)
       	end

        -- Integrate to compute cell-centered value of I.

        points[{i,j}].I_4[m] = (points[{i,j}].S * dx * dy 
        	+ cmath.fabs(angle_values[angle].xi) * dy * upwind_x_value/gamma 
        	+ cmath.fabs(angle_values[angle].eta) * dx * upwind_y_value/gamma)
        	/(points[{i,j}].sigma * dx * dy + cmath.fabs(angle_values[angle].xi) * dy/gamma + cmath.fabs(angle_values[angle].eta) * dx/gamma)

        -- c.printf("x=%d,y=%d,angle=%d I = %lf \n", i, j, angle, points[{i,j}].I_4[m])

        x_faces[{indx+dindx, j}].Ifx_4[m] = (points[{i,j}].I_4[m] - (1-gamma)*upwind_x_value)/gamma
        y_faces[{i, indy+dindy}].Ify_4[m] = (points[{i,j}].I_4[m] - (1-gamma)*upwind_y_value)/gamma

      end
    end
  end
end

task residual(points : region(ispace(int2d), point))
where
  reads (points.I_1, points.I_2, points.I_3, points.I_4, 
  	points.Iiter_1, points.Iiter_2, points.Iiter_3, points.Iiter_4)
do

   -- Compute the residual after each iteration and return the value.

 	var res : double = 0.0
  	var limits = points.bounds

  	for m = 0, N_angles/4 do
    	for i = limits.lo.x, limits.hi.x + 1 do
      		for j = limits.lo.y, limits.hi.y + 1 do
        		res = res + (1.0/(Nx*Ny*N_angles+1))
          			*cmath.pow((points[{i,j}].I_1[m]-points[{i,j}].Iiter_1[m]),2.0)/cmath.pow((points[{i,j}].I_1[m]),2.0)

          		res = res + (1.0/(Nx*Ny*N_angles+1))
          			*cmath.pow((points[{i,j}].I_2[m]-points[{i,j}].Iiter_2[m]),2.0)/cmath.pow((points[{i,j}].I_2[m]),2.0)

          		res = res + (1.0/(Nx*Ny*N_angles+1))
          			*cmath.pow((points[{i,j}].I_3[m]-points[{i,j}].Iiter_3[m]),2.0)/cmath.pow((points[{i,j}].I_3[m]),2.0)

          		res = res + (1.0/(Nx*Ny*N_angles+1))
          			*cmath.pow((points[{i,j}].I_4[m]-points[{i,j}].Iiter_4[m]),2.0)/cmath.pow((points[{i,j}].I_4[m]),2.0)
      		end
    	end
  	end

  	return res
end

task update(points : region(ispace(int2d), point))
where
  reads (points.I_1, points.I_2, points.I_3, points.I_4), 
  reads writes(points.Iiter_1, points.Iiter_2, points.Iiter_3, points.Iiter_4)
do

  	-- Update the intensity before moving to the next iteration.

  	for i in points do 
  		for m = 0, N_angles/4 do
  			points[i].Iiter_1[m] = points[i].I_1[m]
  			points[i].Iiter_2[m] = points[i].I_2[m]
  			points[i].Iiter_3[m] = points[i].I_3[m]
  			points[i].Iiter_4[m] = points[i].I_4[m]
  		end
  	end
end

task reduce_intensity(points : region(ispace(int2d), point),
					angle_values : region(ispace(int1d), angle_value))
where
  reads (points.I_1, points.I_2, points.I_3, points.I_4, angle_values.w),
  reads writes (points.G)
do

   -- Reduce the intensity to summation over all angles

  	var limits = points.bounds
  	for i = limits.lo.x, limits.hi.x+1 do
    	for j = limits.lo.y, limits.hi.y+1 do
      		for m = 0, N_angles/4 do
      			var angle : int = m
        		points[{i,j}].G = points[{i,j}].G + angle_values[angle].w*points[{i,j}].I_1[m]
        		angle = m + N_angles/4
        		points[{i,j}].G = points[{i,j}].G + angle_values[angle].w*points[{i,j}].I_2[m]
        		angle = m + (N_angles/4)*2
        		points[{i,j}].G = points[{i,j}].G + angle_values[angle].w*points[{i,j}].I_3[m]
        		angle = m + (N_angles/4)*3
        		points[{i,j}].G = points[{i,j}].G + angle_values[angle].w*points[{i,j}].I_4[m]        		
      		end
    	end
  	end
end

task create_tecplot_file(points : region(ispace(int2d), point))
where
	reads (points.G),
	reads writes (points.x, points.y)

do

	-- Compute the x- and y-coordinates for vizualization
	var limits = points.bounds
  	for i = limits.lo.x, limits.hi.x+1 do
    	for j = limits.lo.y, limits.hi.y+1 do
        	points[{i,j}].x = (dx * [double](i))
        	points[{i,j}].y = (dy * [double](j))
    	end
  	end

 	-- Tecplot ASCII format (cell-centered)

  	var f = c.fopen("radiation_solver/intensity_parallel.dat", "w")

  	-- Write header

  	c.fprintf(f,'\n\n')
  	c.fprintf(f,'TITLE = "DOM Intensity"\n')
  	c.fprintf(f,'VARIABLES = "X", "Y", "Intensity"\n')
  	c.fprintf(f,'ZONE I= %d J= %d DATAPACKING=BLOCK VARLOCATION=([3]=CELLCENTERED)\n', Nx+1,Ny+1)

  	-- Write the x & y coords, then cell-centered intensity.

  	for i = limits.lo.x, limits.hi.x+1 do
    	for j = limits.lo.y, limits.hi.y+1 do
      		c.fprintf(f,' %.15e ', points[{i,j}].x)
    	end
    	c.fprintf(f,'\n')
 	end

  	for i = limits.lo.x, limits.hi.x+1 do
    	for j = limits.lo.y, limits.hi.y+1 do
      		c.fprintf(f,' %.15e ', points[{i,j}].y)
    	end
    	c.fprintf(f,'\n')
  	end

  	for i = limits.lo.x, limits.hi.x+1 do
    	for j = limits.lo.y, limits.hi.y+1 do
      		c.fprintf(f,' %.15e ', points[{i,j}].G)
    	end
    	c.fprintf(f,'\n')
  	end

  	-- Close the Tecplot file.

  	c.fclose(f)
end

task main()

	-- Some local variables needed for the iterative algorithm.

	var t   : int64  = 1
	var res : double = 1.0
	var N   : int64[1]

	var nt : int64 = 4 -- # tiles per direction

	-- var filename : rawstring = quad_file

	-- Create a region from our grid index space (angles + 2D grid in space)
	-- and our cell field space defined above.

	var grid = ispace(int2d, {x = Nx, y = Ny})
	var points = region(grid, point)

	-- 2D Region from grid index space (+1 in x direction) and x_face field space

	var grid_x = ispace(int2d, {x = Nx+1, y = Ny})

	var x_faces = region(grid_x, x_face)

	-- 2D Region from grid index space (+1 in y direction) and y_face field space

	var grid_y = ispace(int2d, {x = Nx, y = Ny+1})

	var y_faces = region(grid_y, y_face)

	-- 1D Region from angle values

	c.printf(' Number of DOM angles: %d\n', N_angles)
	var angle_indices = ispace(int1d, N_angles)
	var angle_values = region(angle_indices, angle_value)

	-- Initialize all arrays in our field space on the grid. 

	initialize_angle_values(angle_values)

	-- Tile partition cells
	var tiles = ispace(int2d, {x = nt, y = nt})

	var private_cells = make_interior_partition(points, tiles, nt, Nx, Ny)

	-- Extra tile required for ghost
	var nt_face : int64 = nt+1 -- # tiles per direction
	var x_faces_tiles = ispace(int2d, {x = nt_face, y = nt})
	var y_faces_tiles = ispace(int2d, {x = nt, y = nt_face})

	-- Partition faces
	var private_x_faces_hi = make_interior_partition_x_hi(x_faces, x_faces_tiles, nt, Nx+1, Ny)

	var private_y_faces_hi = make_interior_partition_y_hi(y_faces, y_faces_tiles, nt, Nx, Ny+1)

	var private_x_faces_lo = make_interior_partition_x_lo(x_faces, x_faces_tiles, nt, Nx+1, Ny)

	var private_y_faces_lo = make_interior_partition_y_lo(y_faces, y_faces_tiles, nt, Nx, Ny+1)

	for color in tiles do
		initialize(private_cells[color])
	end

	-- for color in x_faces_tiles do
	-- 	initialize_x_faces(private_x_faces_hi[color])
	-- end

	-- for color in y_faces_tiles do
	-- 	initialize_y_faces(private_y_faces_hi[color])
	-- end


	while (res > tol) do
    
	    -- Update the source term (in this problem, isotropic).

	    for color in tiles do
	   		source_term(private_cells[color], angle_values)
	   	end

	   	-- Update the grid boundary intensities.
	
	  	-- Update x faces (west bound/east bound)
	  	for j = [int](tiles.bounds.lo.y), [int](tiles.bounds.hi.y) + 1 do
	  		west_bound(private_x_faces_hi[{x_faces_tiles.bounds.lo.x,j}], angle_values)
	  		east_bound(private_x_faces_lo[{x_faces_tiles.bounds.hi.x,j}], angle_values)
	  	end
	  	
	  	-- Update y faces (north bound/south bound)
	  	for i = [int](tiles.bounds.lo.x), [int](tiles.bounds.hi.x) + 1 do
	  		south_bound(private_y_faces_hi[{i,y_faces_tiles.bounds.lo.y}], angle_values)
	  		north_bound(private_y_faces_lo[{i,y_faces_tiles.bounds.hi.y}], angle_values)
	  	end

	  	-- Perform the sweep for computing new intensities.

	  	-- Quadrant 1 - +x, +y
		for i = tiles.bounds.lo.x, tiles.bounds.hi.x + 1 do
			for j = tiles.bounds.lo.y, tiles.bounds.hi.y + 1 do
			
				sweep_1(private_cells[{i,j}], private_x_faces_lo[{i+1,j}], private_y_faces_lo[{i,j+1}], 
					private_x_faces_lo[{i,j}], private_y_faces_lo[{i,j}], angle_values)
			end
		end 

		-- Quadrant 2 - +x, -y
		for i = tiles.bounds.lo.x, tiles.bounds.hi.x + 1 do
			for j = tiles.bounds.hi.y, tiles.bounds.lo.y - 1, -1 do 

				sweep_2(private_cells[{i,j}], private_x_faces_lo[{i+1,j}], private_y_faces_hi[{i,j}], 
					private_x_faces_lo[{i,j}], private_y_faces_hi[{i,j+1}], angle_values)
			end
		end

		-- Quadrant 3 - -x, +y
		for i = tiles.bounds.hi.x, tiles.bounds.lo.x - 1, -1 do 
			for j = tiles.bounds.lo.y, tiles.bounds.hi.y + 1 do

				sweep_3(private_cells[{i,j}], private_x_faces_hi[{i,j}], private_y_faces_lo[{i,j+1}], 
					private_x_faces_hi[{i+1,j}], private_y_faces_lo[{i,j}], angle_values)
			end
		end

		-- Quadrant 4 - -x, -y
		for i = tiles.bounds.hi.x, tiles.bounds.lo.x - 1, -1 do 
			for j = tiles.bounds.hi.y, tiles.bounds.lo.y - 1, -1 do 

				sweep_4(private_cells[{i,j}], private_x_faces_hi[{i,j}], private_y_faces_hi[{i,j}], 
					private_x_faces_hi[{i+1,j}], private_y_faces_hi[{i,j+1}], angle_values)
			end
		end

  		-- Compute the residual and output to the screen.
  		res = 0.0
  		for color in tiles do
  			res += residual(private_cells[color])
  		end
  		res = cmath.sqrt(res)


  		if (t == 1) then
    		c.printf("\n")
    		c.printf(" Iteration     Residual         \n")
    		c.printf(" ------------------------------ \n")
  		end
  		c.printf( "   %3d    %.15e \n", t, res)

  		-- Update the intensities and the iteration number.

        for color in tiles do
        	update(private_cells[color])
        end

		t = t + 1

		if t > 10 then
			break
		end

	end

  	-- Reduce intensity
    -- reduce_intensity(points, angle_values)

    -- Write a Tecplot file to vizualize the intensity.
    --todo: divide by tile?
    -- create_tecplot_file(points)

end

regentlib.start(main)








