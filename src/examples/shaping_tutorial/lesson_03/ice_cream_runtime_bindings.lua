local dim = 2

local scoop_radius = 1.1
local scoop_lift = 2.0
local sprinkle_lift = 3.0
local cone_lift = -2.0

local function point(x, y, z)
  if dim == 2 then
    return {x, y}
  end
  return {x, y, z or 0.0}
end

local function offset(y)
  return function()
    return point(0.0, y)
  end
end

local function scoop_scale()
  return {scoop_radius}
end

return {
  dimensions = dim,
  scoop_lift = scoop_lift,
  sprinkle_lift = sprinkle_lift,
  cone_lift = cone_lift,
  offset = offset,
  scoop_scale = scoop_scale
}
