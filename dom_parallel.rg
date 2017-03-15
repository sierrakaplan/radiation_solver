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

local quad_file = "S8.dat"

-- Grid size (x cells, y cells)

local Nx = 100
local Ny = 100

--todo: Read from file in Lua

local N_angles = 1

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


-- Fieldspaces

-- Internal cell values are essentially private,
-- face values are what need to be passed to downstream neighbor
-- Update cell value, then update downstream face values

-- quadrature information

-- todo: create 4 regions of angle values
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
	y : double,
}

fspace x_face {
	Ifx_1 : double[N_angles/4]		-- x face intensity per angle
	Ifx_2 : double[N_angles/4]		-- x face intensity per angle
	Ifx_3 : double[N_angles/4]		-- x face intensity per angle
	Ifx_4 : double[N_angles/4]		-- x face intensity per angle
}

fspace y_face {
	Ify_1 : double[N_angles/4] 		-- y face intensity per angle
	Ify_2 : double[N_angles/4] 		-- y face intensity per angle
	Ify_3 : double[N_angles/4] 		-- y face intensity per angle
	Ify_4 : double[N_angles/4] 		-- y face intensity per angle
}

task initialize_faces(x_faces : region(ispace(int2d), x_face),
					  y_faces : region(ispace(int2d), y_face))
where writes(x_faces.Ifx, y_faces,Ify)
do

	for i in x_faces do
		for m = 1, N_angles do
			x_faces[i].Ifx[m] = 0.0
		end
	end

	for i in y_faces do
		for m = 1, N_angles do
			y_faces[i].Ify[m] = 0.0
		end
	end

end

-- todo: generate angle values rather than read from file
task initialize_angle_values(angle_values : region(ispace(int1d), angle_value),
						   filename : rawstring)
where writes (angle_values.xi, angle_values.eta, angle_values.w)

do

	-- Now set the angle_value information, read in from file.

  	-- todo: num angles
  	var N   : int64[1]
  	var val : double[1]

  	var f = c.fopen(filename, "rb")
  	get_number_angles(f, N)

  	var limits = points.bounds

  	for i in angle_values do
  		read_val(f, val)
		angle_values[i].xi = val[0]
	end

	for i in angle_values do
  		read_val(f, val)
		angle_values[i].eta = val[0]
	end

	for i in angle_values do
  		read_val(f, val)
		angle_values[i].w = val[0]
	end
  

 	c.fclose(f)

end

task initialize(points : region(ispace(int2d), point))
where
  reads writes(points.T), writes(points.Ib, points.sigma, points.I, points.G,
                                 points.Iiter, points.S)
do

  	-- First loop over all points to set the constant values.

  	for i in points do

	    -- Blackbody source

	    points[i].T  = 300.0
	    points[i].Ib = (SB/pi)*cmath.pow(points[i].T,4.0)

	    -- Extinction coefficient

	    points[i].sigma = 3.0

	    -- Intensities (cell-based)

	    for m = 1, N_angles do
			points[i].I[m]     = 0.0
	    	points[i].Iiter[m] = 0.0
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
  reads (points.Iiter, points.w, points.Ib, points.sigma, angle_values.w),
  reads writes (points.S)
do

  	-- Get array bounds

  	var  limits = points.bounds

  	-- Loop over all angles and grid cells to compute the source term
  	-- for the current iteration.

  	for i = limits.lo.x, limits.hi.x do
    	for j = limits.lo.y, limits.hi.y do
      		points[{i,j}].S = (1.0-omega)*SB*points[{i,j}].sigma*points[{i,j}].Ib
     	 	for m = 1, N_angles do
        		points[{i,j}].S = points[{i,j}].S + omega*points[{i,j}].sigma/(4.0*pi)*angle_values[m].w*points[{i,j}].Iiter[m]
      		end
    	end
  	end

end

-- Partition x,y grid into tiles
task make_interior_partition(points : region(ispace(int2d), point),
                            tiles : ispace(int2d),
                            ntiles : int64, nx : int64, ny : int64) do
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

task make_interior_partition_x(faces : region(ispace(int2d), x_face),
							tiles : ispace(int2d),
							ntiles : int64, nx : int64, ny : int64) do
	
	var coloring = c.legion_domain_point_coloring_create()
	for tile in tiles do
	    var lo = int2d { x = tile.x * nx / ntiles, y = tile.y * ny / ntiles}
		var hi = int2d { x = (tile.x + 1) * nx / ntiles - 1, y = (tile.y + 1) * ny / ntiles - 1}
		var rect = rect2d {lo = lo, hi = hi}
		c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
	end
	var p = partition(disjoint, faces, coloring, tiles)
	c.legion_domain_point_coloring_destroy(coloring)
	return p

end

-- todo: repetition for all interior partitions, way to avoid?
task make_interior_partition_y(faces : region(ispace(int2d), y_face),
							tiles : ispace(int2d),
							ntiles : int64, nx : int64, ny : int64) do
	
	var coloring = c.legion_domain_point_coloring_create()
	for tile in tiles do
	    var lo = int2d { x = tile.x * nx / ntiles, y = tile.y * ny / ntiles}
		var hi = int2d { x = (tile.x + 1) * nx / ntiles - 1, y = (tile.y + 1) * ny / ntiles - 1}
		var rect = rect2d {lo = lo, hi = hi}
		c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
	end
	var p = partition(disjoint, faces, coloring, tiles)
	c.legion_domain_point_coloring_destroy(coloring)
	return p

end


-- The ghost region is the right most column of each tile 
task make_ghost_partition_x(faces : region(ispace(int2d), x_face),
							tiles : ispace(int2d),
							ntiles : int64, nx : int64, ny : int64) do
	var coloring = c.legion_domain_point_coloring_create()

	for tile in tiles do
		var lo = int2d { x = (tile.x + 1) * n / ntiles - 1, y = tile.y * n / ntiles}
		var hi = int2d { x = (tile.x + 1) * n / ntiles - 1, y = (tile.y + 1) * n / ntiles - 1}
		var rect = rect2d {lo = lo, hi = hi}
		c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
	end

	var p = partition(disjoint, faces, coloring, tiles)
	c.legion_domain_point_coloring_destroy(coloring)
	return p

end

-- The ghost region is the bottom row of each tile
task make_ghost_partition_y(faces : region(ispace(int2d), y_face),
							tiles : ispace(int2d),
							ntiles : int64, nx : int64, ny : int64) do
	var coloring = c.legion_domain_point_coloring_create()

	for tile in tiles do
		var lo = int2d { x = tile.x * n / ntiles, y = (tile.y + 1) * n / ntiles - 1}
		var hi = int2d { x = (tile.x + 1) * n / ntiles - 1, y = (tile.y + 1) * n / ntiles - 1}
		var rect = rect2d {lo = lo, hi = hi}
		c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
	end

	var p = partition(disjoint, faces, coloring, tiles)
	c.legion_domain_point_coloring_destroy(coloring)
	return p

end

-- pass all 4 quadrants of faces and all 4 angle_values
task west_bound(faces : region(ispace(int2d), x_face),
				angle_values : region(ispace(int1d), angle_value))
where
	reads (angle_values.w, angle_values.xi),
	reads writes (faces.Ifx) do

	-- Get array bounds

  	var limits = points.bounds

 	 -- Temporary variables for the west bound

  	var reflect : double = 0.0
  	var epsw    : double = emiss_west
  	var Tw      : double = T_west

  	for j = limits.lo.y, limits.hi.y do
  		reflect = 0
  		for m = 1, N_angles do
  			if angle_values[m].xi < 0 then
  				reflect = reflect + (1.0-epsw)/pi*angle_values[m].w*cmath.fabs(angle_values[m].xi)*faces[{limits.lo.x,j}].Ifx[m]
  			end
  		end
  		for m = 1, N_angles do
  			if angle_values[m].xi > 0 then
  				faces[{0,j}].Ifx = epsw*SB*cmath.pow(Tw,4.0)/pi + reflect
  			end
  		end
  	end
end

task east_bound(faces : region(ispace(int2d), x_face),
				angle_values : region(ispace(int1d), angle_value))
where
	reads (angle_values.w, angle_values.xi),
	reads writes (faces.Ifx) do

	-- Get array bounds

  	var limits = points.bounds

 	 -- Temporary variables for the west bound

  	var reflect : double = 0.0
  	var epsw    : double = emiss_east
  	var Tw      : double = T_east

  	for j = limits.lo.y, limits.hi.y do
  		reflect = 0
  		for m = 1, N_angles do
  			if angle_values[m].xi < 0 then
  				reflect = reflect + (1.0-epsw)/pi*angle_values[m].w*cmath.fabs(angle_values[m].xi)*faces[{limits.hi.x,j}].Ifx
  			end
  		end
  		for m = 1, N_angles do
  			if angle_values[m].xi > 0 then
  				faces[{0,j}].Ifx = epsw*SB*cmath.pow(Tw,4.0)/pi + reflect
  			end
  		end
  	end
end

task north_bound(faces : region(ispace(int2d), y_face),
				angle_values : region(ispace(int1d), angle_value))
where
	reads (angle_values.w, angle_values.eta),
	reads writes (faces.Ify) do

	-- Get array bounds

  	var limits = points.bounds

 	 -- Temporary variables for the west bound

  	var reflect : double = 0.0
  	var epsw    : double = emiss_north
  	var Tw      : double = T_north

  	for i = limits.lo.x, limits.hi.x do
  		reflect = 0
  		for m = 1, N_angles do
  			if angle_values[m].eta < 0 then
  				reflect = reflect + (1.0-epsw)/pi*angle_values[m].w*cmath.fabs(angle_values[m].eta)*faces[{i,limits.hi.y}].Ify
  			end
  		end
  		for m = 1, N_angles do
  			if angle_values[m].eta > 0 then
  				faces[{0,j}].Ify = epsw*SB*cmath.pow(Tw,4.0)/pi + reflect
  			end
  		end
  	end
end

task south_bound(faces : region(ispace(int2d), y_face),
				angle_values : region(ispace(int1d), angle_value))
where
	reads (angle_values.w, angle_values.eta)
	reads writes (faces.Ify)

	-- Get array bounds

  	var limits = points.bounds

 	 -- Temporary variables for the west bound

  	var reflect : double = 0.0
  	var epsw    : double = emiss_south
  	var Tw      : double = T_south

  	for i = limits.lo.x, limits.hi.x do
  		reflect = 0
  		for m = 1, N_angles do
  			if angle_values[m].eta < 0 then
  				reflect = reflect + (1.0-epsw)/pi*angle_values[m].w*cmath.fabs(angle_values[m].eta)*faces[{i, limits.lo.y}].Ify
  			end
  		end
  		for m = 1, N_angles do
  			if angle_values[m].eta > 0 then
  				faces[{0,j}].Ify = epsw*SB*cmath.pow(Tw,4.0)/pi + reflect
  			end
  		end
  	end
end

task sweep(points : region(ispace(int2d), point),
		   x_faces : region(ispace(int2d), x_face),
		   y_faces : region(ispace(int2d), y_face),
		   ghost_x_faces : region(ispace(int2d), x_face),
		   ghost_y_faces : region(ispace(int2d), y_face),
		   angle_values : region(ispace(int1d), angle_value))
where
  reads (angle_values.xi, angle_values.eta, points.S, ghost_x_faces.Ifx, ghost_y_faces.Ify),
  reads writes(points.I, x_faces.Ifx, y_faces.Ify)
do

  -- Get array bounds and some temporary index variables for sweeping

  var limits = points.bounds

  var indx   : int64 = 0
  var indy   : int64 = 0

  var dindx  : int64 = 0
  var startx : int64 = 0
  var endx   : int64 = 0

  var dindy  : int64 = 0
  var starty : int64 = 0
  var endy   : int64 = 0

  -- Outer loop over all angles.

  for m = 1, N_angles do

    -- Determine our sweeping direction. If xi > 0, angle points in +x,
    -- so we sweep from left to right. Otherwise, angle points in -x, 
    -- so sweep from right to left.

    if (angle_values[m].xi > 0) then
      dindx  = 1
      startx = limits.lo.x
      endx   = limits.hi.x
    else
      dindx  = -1
      startx = limits.hi.x - 1
      endx   = limits.lo.x - 1
    end

    -- If eta > 0, angle points in +y, sp sweep bottom to top.
    -- Otherwise, angle points in -y, so sweep from top to bottom.

    if (angle_values[m].eta > 0) then
      dindy  = 1
      starty = 0
      endy   = limits.hi.y
    else
      dindy  = -1
      starty = limits.hi.y - 1
      endy   = -1
    end

    -- Use our direction and increments for the sweep.

    for j = starty,endy,dindy do
      for i = startx,endx,dindx do

        -- Index of upwind x-face & y-face intensity

        indx = i - min(dindx,0)
        indy = j - min(dindy,0)

        -- Integrate to compute cell-centered value of I.

        points[{i,j}].I = (points[{i,j}].S*dx*dy + cmath.fabs(angle_values[m].xi)*dy*x_faces[{indx,j}].Ifx/gamma + cmath.fabs(angle_values[m].eta)*dx*y_faces[{i,indy}].Ify/gamma)
        /(sigma*dx*dy + cmath.fabs(angle_values[m].xi)*dy/gamma + cmath.fabs(angle_values[m].eta)*dx/gamma)

        -- Compute downwind intensities on cell faces.

        x_faces[{indx+dindx,j}].Ifx = (points[{i,j}].I - (1-gamma)*x_faces[{indx,j}].Ifx)/gamma
        y_faces[{i,indy+dindy}].Ify = (points[{i,j}].I - (1-gamma)*y_faces[{i,indy}].Ify)/gamma

      end
    end
  end

end

-- partition rect with hi < lo
task empty_region() do

end

task main()

	-- Some local variables needed for the iterative algorithm.

	var t   : int64  = 1
	var res : double = 1.0
	var N   : int64[1]

	var nt : int64 = 4 -- # tiles per direction

	-- Check the file containing the angle_values

	-- todo: using constant number of angles

	var filename : rawstring = quad_file
	var f = c.fopen(filename, "rb")
	get_number_angles(f, N)
	c.fclose(f)
	c.printf(' Number of DOM angles: %d\n', N[0])
	N_angles = N[0]

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

	var angle_indices = ispace(int1d, N_angles)

	var angle_values = region(angle_indices, angle_value)

	-- Initialize all arrays in our field space on the grid. 

	initialize(points, filename)

	initialize_faces(x_faces, y_faces)

	-- Tile partition cells
	var tiles = ispace(int2d, {x = nt, y = nt})
	var private_cells = make_interior_partition(points, tiles, nt, Nx, Ny)

	-- Partition faces
	var private_x_faces = make_interior_partition_x(x_faces, tiles, nt, Nx+1, Ny)

	var private_y_faces = make_interior_partition_y(y_faces, tiles, nt, Nx, Ny+1)

	var ghost_x_faces = make_ghost_partition_x(x_faces, tiles, nt, Nx+1, Ny)

	var ghost_y_faces = make_ghost_partition_y(y_faces, tiles, nt, Nx, Ny+1)


	while (res > tol) do
    
    -- for all quadrants?

	    -- Update the source term (in this problem, isotropic).
	    for color in tiles do
	   		source_term(private_cells[color])
	   	end
	-- end

	   	-- Update the grid boundary intensities.
	   	-- pass in all quadrants for both faces and angle values
	  	-- Update x faces (west bound/east bound)
	  	for j = [int](tiles.bounds.lo.y), [int](tiles.bounds.hi.y) + 1 do
	  		west_bound(private_x_faces[0][j], angle_values)
	  		east_bound(private_x_faces[Nx][j], angle_values)
	  	end
	  	
	  	-- Update y faces (north bound/south bound)
	  	for i = [int](tiles.bounds.lo.x), [int](tiles.bounds.hi.x) + 1 do
	  		north_bound(private_y_faces[i][0], angle_values)
	  		south_bound(private_y_faces[i][Ny], angle_values)
	  	end

	-- for all quadrants
		for color in tiles do

			-- Perform the sweep for computing new intensities.
			-- todo: empty regions?

			var ghost_x_region = empty_region()
			if color.y > 0 do
				ghost_x_region = ghost_x_faces[color.x][color.y-1]
			end
			var ghost_y_region = empty_region()
			if color.x > 0 do 
				ghost_y_region = ghost_y_faces[color.x-1][color.y]
			end
			sweep(private_cells[color], private_x_faces[color], private_y_faces[color], 
				ghost_x_region, ghost_y_region, angle_values)
		end  

  
  	-- end all quadrants 

  	-- for all quadrants
  		-- residual & update 
  	-- end
  		-- 

    
		t = t + 1

	end

  -- Write a Tecplot file to vizualize the intensity.

  -- todo: reduce_intensity(points)


  -- read n angles from file in lua
  -- change ghost regions to be 1 over
  -- update & residual for each tile
  -- loop over quadrants before source term, boundaries, sweep, and residual
  -- change sweep code to take which range of angles (copy it 4 times)
  -- create 4 regions of points (each direction of angles) and faces  and angle values
  	-- BEFORE partitioning (partition each separately)

  -- create 4 different partitioning code so that 1,2,3,4 is the right order (start with corner 
  	-- with angle direction always)




  -- 1 quadrangle
  	-- Start with 1 angle
  		-- * Create region of angle weights 
  		-- * Partition cells into tiles
  		-- * Partition each face into tiles
  		-- * Partition each face into ghost regions per tile
  		-- Pass empty regions for boundary cases ghost regions
  		-- * West, East, South, North, Sweep - per tile
  	-- Group 12 angles

  -- (work on 2nd quadrangle, then try meta programming to do all 4)

  -- create two slides on work
  -- create design doc on partitions
  -- poster by May

  -- Compute the residual and output to the screen.
    	--todo:
  		-- res = residual(points, t)
    


    -- Update the intensities and the iteration number.
        
	--  todo: update(points)

end

regentlib.start(main)










