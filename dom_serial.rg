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

local quad_file = "radiation_solver/S4.dat"

-- Grid size (x cells, y cells)

local Nx = 100
local Ny = 100

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

-- Create our fieldspace

-- Internal cell values are essentially private,
-- face values are what need to be passed to downstream neighbor
-- Update cell value, then update downstream face values

-- Another option for launching tasks is index space launch, which you can't
-- do if you have alias partitions

-- fieldspace for x and y values (grid) n+1 x n+1
-- fieldspace for cell values (grid) n x n

-- boundary conditions faces
-- cell update taking upstream faces
-- face update downstream taking upstream cells


fspace point {
  x     : double,
  y     : double,
  xi    : double,
  eta   : double,
  w     : double,
  T     : double,
  Ib    : double,
  sigma : double,
  I     : double, --  cell center
  Ifx   : double, --  face values (could switch to be downstream rather than upstream)
  Ify   : double, --  face values
  Iiter : double,
  S     : double,
  G     : double
}

terra get_number_angles(f : &c.FILE, N : &int64)
  return c.fscanf(f, "%d\n", &N[0])
end

terra read_val(f : &c.FILE, val : &double)
  return c.fscanf(f, "%lf\n", &val[0])
end

task initialize(points : region(ispace(int3d), point),
                filename : rawstring)
where
  reads writes(points.T), writes(points.Ib, points.sigma, points.I, points.G,
                                 points.Ifx, points.Ify, points.Iiter,
                                 points.S, points.xi, points.eta, points.w)
do

  -- First loop over all points to set the constant values.

  for i in points do

    -- Blackbody source

    points[i].T  = 300.0
    points[i].Ib = (SB/pi)*cmath.pow(points[i].T,4.0)

    -- Extinction coefficient

    points[i].sigma = 3.0

    -- Intensities (cell- and face-based)

    points[i].I     = 0.0
    points[i].Ifx   = 0.0
    points[i].Ify   = 0.0
    points[i].Iiter = 0.0
    points[i].G     = 0.0

    -- Source term

    points[i].S = 0.0

  end

  -- Now set the quadrature information. Note that this is hard-coded for
  -- now but could be read in from a file instead.

  var N   : int64[1]
  var val : double[1]

  var f = c.fopen(filename, "rb")
  get_number_angles(f, N)

  var limits = points.bounds

  for m = limits.lo.x, limits.hi.x + 1 do
    read_val(f, val)
    for i = limits.lo.y, limits.hi.y + 1 do
      for j = limits.lo.z, limits.hi.z + 1 do
        points[{m,i,j}].xi = val[0]
      end
    end
  end

  for m = limits.lo.x, limits.hi.x + 1 do
    read_val(f, val)
    for i = limits.lo.y, limits.hi.y + 1 do
      for j = limits.lo.z, limits.hi.z + 1 do
        points[{m,i,j}].eta = val[0]
      end
    end
  end

  for m = limits.lo.x, limits.hi.x + 1 do
    read_val(f, val)
    for i = limits.lo.y, limits.hi.y + 1 do
      for j = limits.lo.z, limits.hi.z + 1 do
        points[{m,i,j}].w = val[0]
      end
    end
  end

  c.fclose(f)

end

task source_term(points : region(ispace(int3d), point))
where
  reads (points.Iiter, points.w, points.Ib, points.sigma),
  reads writes (points.S)
do

  -- Get array bounds

  var limits = points.bounds

  -- Loop over all angles and grid cells to compute the source term
  -- for the current iteration.

  for i = limits.lo.y, limits.hi.y do
    for j = limits.lo.z, limits.hi.z do
      points[{0,i,j}].S = (1.0-omega)*SB*points[{0,i,j}].sigma*points[{0,i,j}].Ib
      for m = limits.lo.x, limits.hi.x + 1 do
        points[{0,i,j}].S = points[{0,i,j}].S + omega*points[{0,i,j}].sigma/(4.0*pi)*points[{m,0,0}].w*points[{m,i,j}].Iiter
      end
    end
  end

end

task west_bound(points : region(ispace(int3d), point))
where
  reads (points.w, points.xi, points.sigma),
  reads writes (points.Ifx)
do

  -- Get array bounds

  var limits = points.bounds

  -- Temporary variables for the west bound

  var reflect : double = 0.0
  var epsw    : double = emiss_west
  var Tw      : double = T_west

  -- Loop over the west boundary

  for j = limits.lo.z, limits.hi.z do
    reflect = 0
    for m = limits.lo.x, limits.hi.x + 1 do
      if points[{m,0,0}].xi < 0 then
        reflect = reflect + (1.0-epsw)/pi*points[{m,0,0}].w*cmath.fabs(points[{m,0,0}].xi)*points[{m,0,j}].Ifx
      end
    end
    for m = limits.lo.x, limits.hi.x + 1 do
      if points[{m,0,0}].xi > 0 then
        points[{m,0,j}].Ifx = epsw*SB*cmath.pow(Tw,4.0)/pi + reflect
      end
    end
  end

end

task east_bound(points : region(ispace(int3d), point))
where
  reads (points.w, points.xi, points.sigma),
  reads writes(points.Ifx)
do

  -- Get array bounds

  var limits = points.bounds

  -- Temporary variables for the east bound

  var reflect : double = 0.0
  var epsw    : double = emiss_east
  var Tw      : double = T_east

  -- Loop over the east boundary

  for j = limits.lo.z, limits.hi.z do
    reflect = 0
    for m = limits.lo.x, limits.hi.x + 1 do
      if points[{m,0,0}].xi > 0 then
        reflect = reflect + (1.0-epsw)/pi*points[{m,0,0}].w*points[{m,0,0}].xi*points[{m,Nx,j}].Ifx
      end
    end
    for m = limits.lo.x, limits.hi.x + 1 do
      if points[{m,0,0}].xi < 0 then
        points[{m,Nx,j}].Ifx = epsw*SB*cmath.pow(Tw,4.0)/pi + reflect
      end
    end
  end

end

task south_bound(points : region(ispace(int3d), point))
where
  reads (points.w, points.eta, points.sigma),
  reads writes(points.Ify)
do

  -- Get array bounds

  var limits = points.bounds

  -- Temporary variables for the south bound

  var reflect : double = 0.0
  var epsw    : double = emiss_south
  var Tw      : double = T_south

  -- Loop over the south boundary

  for i = limits.lo.y, limits.hi.y do
    reflect = 0
    for m = limits.lo.x, limits.hi.x + 1 do
      if points[{m,0,0}].eta < 0 then
        reflect = reflect + (1.0-epsw)/pi*points[{m,0,0}].w*cmath.fabs(points[{m,0,0}].eta)*points[{m,i,0}].Ify
      end
    end
    for m = limits.lo.x, limits.hi.x + 1 do
      if points[{m,0,0}].eta > 0 then
        points[{m,i,0}].Ify = epsw*SB*cmath.pow(Tw,4.0)/pi + reflect
      end
    end
  end

end

task north_bound(points : region(ispace(int3d), point))
where
  reads (points.w, points.eta, points.sigma),
  reads writes(points.Ify)
do

  -- Get array bounds

  var limits = points.bounds

  -- Temporary variables for the north bound

  var reflect : double = 0.0
  var epsw    : double = emiss_north
  var Tw      : double = T_north

  -- Loop over the north boundary

  for i = limits.lo.y, limits.hi.y do
    reflect = 0
    for m = limits.lo.x, limits.hi.x + 1 do
      if points[{m,0,0}].eta > 0 then
        reflect = reflect + (1.0-epsw)/pi*points[{m,0,0}].w*points[{m,0,0}].eta*points[{m,i,Ny}].Ify
      end
    end
    for m = limits.lo.x, limits.hi.x + 1 do
      if points[{m,0,0}].eta < 0 then
        points[{m,i,Ny}].Ify = epsw*SB*cmath.pow(Tw,4.0)/pi + reflect
      end
    end
  end

end


task sweep(points : region(ispace(int3d), point))
where
  reads (points.xi, points.eta, points.sigma, points.S),
  reads writes(points.I, points.Ifx, points.Ify)
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

  for m = limits.lo.x, limits.hi.x + 1 do


    -- Determine our sweeping direction. If xi > 0, angle points in +x,
    -- so we sweep from left to right. Otherwise, angle points in -x, 
    -- so sweep from right to left.

    if (points[{m,0,0}].xi > 0) then
      dindx  = 1
      startx = 0
      endx   = Nx
      -- c.printf("+x")
    else
      dindx  = -1
      startx = Nx-1
      endx   = -1
      -- c.printf("-x")
    end

    -- If eta > 0, angle points in +y, sp sweep bottom to top.
    -- Otherwise, angle points in -y, so sweep from top to bottom.

    if (points[{m,0,0}].eta > 0) then
      dindy  = 1
      starty = 0
      endy   = Ny
      -- c.printf("+y\n")
    else
      dindy  = -1
      starty = Ny-1
      endy   = -1
      -- c.printf("-y\n")
    end

    -- Use our direction and increments for the sweep.

    for j = starty,endy,dindy do
      for i = startx,endx,dindx do

        -- Index of upwind x-face & y-face intensity

        indx = i - min(dindx,0)
        indy = j - min(dindy,0)

        -- Integrate to compute cell-centered value of I.

        points[{m,i,j}].I = (points[{0,i,j}].S*dx*dy 
          + cmath.fabs(points[{m,0,0}].xi)*dy*points[{m,indx,j}].Ifx/gamma 
          + cmath.fabs(points[{m,0,0}].eta)*dx*points[{m,i,indy}].Ify/gamma)
        /(points[{0,i,j}].sigma*dx*dy 
          + cmath.fabs(points[{m,0,0}].xi)*dy/gamma 
          + cmath.fabs(points[{m,0,0}].eta)*dx/gamma)


        -- if dindy > 0 and dindx < 0 then
          -- c.printf("x=%d,y=%d,angle=%d I = %lf \n", i, j, m, points[{m,i,j}].I)
        -- end


        -- Compute downwind intensities on cell faces.

        points[{m,indx+dindx,j}].Ifx = (points[{m,i,j}].I - (1-gamma)*points[{m,indx,j}].Ifx)/gamma
        points[{m,i,indy+dindy}].Ify = (points[{m,i,j}].I - (1-gamma)*points[{m,i,indy}].Ify)/gamma

      end
    end
  end

end

task residual(points : region(ispace(int3d), point),
              t : int64)
where
  reads (points.I, points.Iiter)
do

   -- Compute the residual after each iteration and return the value.

  var res : double = 0.0
  var limits = points.bounds
  for m = limits.lo.x, limits.hi.x + 1 do
    for i = limits.lo.y, limits.hi.y do
      for j = limits.lo.z, limits.hi.z do
        res = res + (1.0/(Nx*Ny*(limits.hi.x+1)))
          *cmath.pow((points[{m,i,j}].I-points[{m,i,j}].Iiter),2.0)/cmath.pow((points[{m,i,j}].I),2.0)
      end
    end
  end
  res = cmath.sqrt(res)

  if (t == 1) then
    c.printf("\n")
    c.printf(" Iteration     Residual         \n")
    c.printf(" ------------------------------ \n")
  end
  c.printf( "   %3d    %.15e \n", t, res )

  return res

end

task update(points : region(ispace(int3d), point))
where
  reads (points.I), writes(points.Iiter)
do

  -- Update the intensity before moving to the next iteration.

  for i in points do points[i].Iiter = points[i].I end

end

task reduce_intensity(points : region(ispace(int3d), point))
where
  reads (points.I, points.w),
  reads writes (points.x, points.y, points.G)
do

   -- Reduce the intensity to summation over all angles

  var limits = points.bounds
  for i = limits.lo.y, limits.hi.y do
    for j = limits.lo.z, limits.hi.z do
      for m = limits.lo.x, limits.hi.x + 1 do
        points[{0,i,j}].G = points[{0,i,j}].G + points[{m,i,j}].w*points[{m,i,j}].I
        -- c.printf("x=%d,y=%d,angle=%d I = %lf \n", i, j, m, points[{m,i,j}].I)
      end
    end
  end

  -- Compute the x- and y-coordinates for vizualization

  for i = limits.lo.y, limits.hi.y+1 do
    for j = limits.lo.z, limits.hi.z+1 do
        points[{0,i,j}].x = (dx * [double](i))
        points[{0,i,j}].y = (dy * [double](j))
    end
  end

  -- Tecplot ASCII format (cell-centered)

  var f = c.fopen("radiation_solver/intensity_serial.dat", "w")

  -- Write header

  c.fprintf(f,'\n\n')
  c.fprintf(f,'TITLE = "DOM Intensity"\n')
  c.fprintf(f,'VARIABLES = "X", "Y", "Intensity"\n')
  c.fprintf(f,'ZONE I= %d J= %d DATAPACKING=BLOCK VARLOCATION=([3]=CELLCENTERED)\n', Nx,Ny)

  -- Write the x & y coords, then cell-centered intensity.

  for i = limits.lo.y, limits.hi.y do
    for j = limits.lo.z, limits.hi.z do
      c.fprintf(f,' %.15e ', points[{0,i,j}].x)
    end
    c.fprintf(f,'\n')
  end

  for i = limits.lo.y, limits.hi.y do
    for j = limits.lo.z, limits.hi.z do
      c.fprintf(f,' %.15e ', points[{0,i,j}].y)
    end
    c.fprintf(f,'\n')
  end

  for i = limits.lo.y, limits.hi.y do
    for j = limits.lo.z, limits.hi.z do
      c.fprintf(f,' %.15e ', points[{0,i,j}].G)
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

  -- Check the file containing the quadrature info for the length of
  -- the array of angles. We will use this to create the 3D index space.

  var filename : rawstring = quad_file
  var f = c.fopen(filename, "rb")
  get_number_angles(f, N)
  c.fclose(f)
  c.printf(' Number of DOM angles: %d\n', N[0])

  -- Create our grid index space. Here, we are using a 3D index space to
  -- defined a set of angles and a 2D grid in space. Maybe there is a better
  -- way to do this with separate regions? Note the +1 for the grid directions, 
  -- these are due to the face intensity arrays being face-based and needing 
  -- an additional slot.

  var grid = ispace(int3d, { x = N[0], y = Nx+1, z = Ny+1 })

  -- Create a region from our grid index space (angles + 2D grid in space)
  -- and our point field space defined above.

  var points = region(grid, point)

  -- Initialize all arrays in our field space on the grid. Note that some
  -- arrays only truly need to be along all angles or on the grid. In these
  -- cases, we will have extra unused data, and there is likely a better way
  -- to take this into account.

  initialize(points, filename)

  while (res > tol) do
    
    -- Update the source term (in this problem, isotropic).
    
    source_term(points)
    
    -- Update the grid boundary intensities.
    
    west_bound(points)
    east_bound(points)

    south_bound(points)
    north_bound(points)
        
    -- Perform the sweep for computing new intensities.
    
    sweep(points)
    
    -- Compute the residual and output to the screen.
    
    res = residual(points, t)

    -- Update the intensities and the iteration number.
        
    update(points)
    t = t + 1

    -- if t > 2 then
    --   break
    -- end

  end

  -- Write a Tecplot file to vizualize the intensity.

  reduce_intensity(points)

end

regentlib.start(main)
