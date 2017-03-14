--
--        |  1  2  3  4  5 |
--        |  6  7  8  9  10|
--        |  11 12 13 14 15 |
--        |  16 17 18 19 20 |
--        |  21 22 23 24 25 |
--	A(i,j) = (A(i-1,j) + A(i,j-1) + A(i-1,j-1))/3 (average of three to the upper left)

import "regent"

local c = regentlib.c

fspace point  {
	input : double,
	output : double,
}

task check(points : region(ispace(int2d), point))
where reads(points.{input, output}) do
	var rect = points.bounds
	-- input
	c.printf("Input: \n")
	for row = [int](rect.lo.y), [int](rect.hi.y) + 1 do
		for col = [int](rect.lo.x), [int](rect.hi.x) + 1 do
			c.printf("{%.2f} ", points[{row,col}].input)
		end
		c.printf("\n")
	end
	-- output
	c.printf("\nOutput: \n")
	for row = [int](rect.lo.y), [int](rect.hi.y) + 1 do
		for col = [int](rect.lo.x), [int](rect.hi.x) + 1 do
			c.printf("{%.2f} ", points[{row,col}].output)
		end
		c.printf("\n")
	end
end

task make_tile_partition(points : region(ispace(int2d), point),
                            tiles : ispace(int2d),
                            n : int64, ntiles : int64)
  var coloring = c.legion_domain_point_coloring_create()
  for tile in tiles do
    var lo = int2d { x = tile.x * n / ntiles, y = tile.y * n / ntiles}
	var hi = int2d { x = (tile.x + 1) * n / ntiles - 1, y = (tile.y + 1) * n / ntiles - 1}
	var rect = rect2d {lo = lo, hi = hi}
	c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
  end
  var p = partition(disjoint, points, coloring, tiles)
  c.legion_domain_point_coloring_destroy(coloring)
  return p
end

-- Partitions right column of each tile (left ghost region)
task make_ghost_left_partition(points : region(ispace(int2d), point),
                            tiles : ispace(int2d),
                            n : int64, ntiles : int64)

	var coloring = c.legion_domain_point_coloring_create()
	for tile in tiles do
	    var lo = int2d { x = (tile.x + 1) * n / ntiles - 1, y = tile.y * n / ntiles}
		var hi = int2d { x = (tile.x + 1) * n / ntiles - 1, y = (tile.y + 1) * n / ntiles - 1}
		var rect = rect2d {lo = lo, hi = hi}
		c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
	end
  	var p = partition(disjoint, points, coloring, tiles)
  	c.legion_domain_point_coloring_destroy(coloring)
  	return p
end

-- Partitions bottom row of each tile (upper ghost region)
task make_ghost_upper_partition(points : region(ispace(int2d), point),
                            tiles : ispace(int2d),
                            n : int64, ntiles : int64)

	var coloring = c.legion_domain_point_coloring_create()
	for tile in tiles do
	    var lo = int2d { x = tile.x * n / ntiles, y = (tile.y + 1) * n / ntiles - 1}
		var hi = int2d { x = (tile.x + 1) * n / ntiles - 1, y = (tile.y + 1) * n / ntiles - 1}
		var rect = rect2d {lo = lo, hi = hi}
		c.legion_domain_point_coloring_color_domain(coloring, tile, rect)
	end
  	var p = partition(disjoint, points, coloring, tiles)
  	c.legion_domain_point_coloring_destroy(coloring)
  	return p
end

task average_interior(tile: region(ispace(int2d), point))
where reads (tile.input), reads writes (tile.output) do

	var rect = tile.bounds

	-- row, col = y, x
	-- Interior calculations
	for row = [int](rect.lo.y) + 1, [int](rect.hi.y) + 1 do
		for col = [int](rect.lo.x) + 1, [int](rect.hi.x) + 1 do
			tile[{row,col}].output = (tile[{row,col}].input + tile[{row-1,col}].output + 
				tile[{row,col-1}].output + tile[{row-1,col-1}].output)/4;
		end
	end


end

task average_upper_row(tile: region(ispace(int2d), point),
			tile_upper: region(ispace(int2d), point))
where reads (tile.input), reads (tile_upper.output), reads writes (tile.output) do

	-- row, col = y, x

	var rect = tile.bounds

	var upper_row_index = tile_upper.bounds.hi.y
	var upper_row_index_col = tile_upper.bounds.lo.x + 1

	var row = rect.lo.y
	for col = [int](rect.lo.x) + 1, [int](rect.hi.x) + 1 do
		tile[{row, col}].output = (tile[{row, col}].input 
			+ tile[{row, col-1}].output 
			+ tile_upper[{upper_row_index, upper_row_index_col}].output
			+ tile_upper[{upper_row_index, upper_row_index_col-1}].output)/4
		upper_row_index_col += 1
	end
end

task average_left_col(tile: region(ispace(int2d), point),
			tile_left: region(ispace(int2d), point))
where reads (tile.input), reads (tile_left.output), reads writes (tile.output) do

	-- row, col = y, x

	var rect = tile.bounds

	var left_col_index = tile_left.bounds.hi.x
	var left_col_index_row = tile_left.bounds.lo.y + 1

	var col = rect.lo.x
	for row = [int](rect.lo.y) + 1, [int](rect.hi.y) + 1 do
		tile[{row, col}].output = (tile[{row, col}].input 
			+ tile[{row-1, col}].output 
			+ tile_left[{left_col_index_row, left_col_index}].output
			+ tile_left[{left_col_index_row-1, left_col_index}].output)/4
		left_col_index_row += 1
	end
end

task average_upper_corner(tile: region(ispace(int2d), point),
			tile_left: region(ispace(int2d), point),
			tile_upper: region(ispace(int2d), point),
			tile_corner: region(ispace(int2d), point))
where reads (tile.input), reads (tile_left.output), reads (tile_upper.output), 
reads (tile_corner.output), reads writes (tile.output) do

	-- row, col = y, x
	var row = tile.bounds.lo.y
	var col = tile.bounds.lo.x

	tile[{row, col}].output = (tile[{row, col}].input 
		+ tile_left[{row, col-1}].output 
		+ tile_upper[{row-1, col}].output
		+ tile_corner[{row-1, col-1}].output)/4
end

-- todo: index space launch
task main()
	var n : int64 = 10 -- grid size along each dimension
	var ntiles : int64 = 2 -- number of regions in each dimension
	var tsteps : int64 = 1

	-- Grid must be larger than # of tiles
	regentlib.assert(n >= ntiles, "grid too small")

	var grid = ispace(int2d, {x = n, y = n})
	var tiles = ispace(int2d, {x = ntiles, y = ntiles})


	var points = region(grid, point)
	var tile_partition = make_tile_partition(points, tiles, n, ntiles)
	var ghost_left = make_ghost_left_partition(points, tiles, n, ntiles)
	var ghost_upper = make_ghost_upper_partition(points, tiles, n, ntiles) 

	--	Loop over matrix, initialize grid
	for row = [int] (0), [int] (n) do
		for col = [int] (0), [int] (n) do
			var value : double = row*n + col
			points[{row,col}].{input, output} = value
		end
	end

	for t = 0, tsteps do
		for i = [int](tiles.bounds.lo.x), [int](tiles.bounds.hi.x) + 1 do
			for j = [int](tiles.bounds.lo.y), [int](tiles.bounds.hi.y) + 1 do

				if i > tiles.bounds.lo.x and j > tiles.bounds.lo.y then
					average_upper_corner(tile_partition[{i,j}], ghost_left[{i-1,j}],ghost_upper[{i,j-1}],ghost_left[{i-1,j-1}])
				end
				if i > tiles.bounds.lo.x then 
					average_left_col(tile_partition[{i,j}], ghost_left[{i-1,j}])
				end 
				if j > tiles.bounds.lo.y then
					average_upper_row(tile_partition[{i,j}], ghost_upper[{i,j-1}])
				end
				average_interior(tile_partition[{i,j}])
			end
		end
	end

	-- Print/check result
	check(points)

end

regentlib.start(main)