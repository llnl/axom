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

dimensions = dim

shapes = {
  {
    name = "background",
    material = "background",
    geometry = {
      format = "none"
    }
  },
  {
    name = "vanilla_scoop",
    material = "ice_cream",
    geometry = {
      format = "mfem",
      path = "ice_cream_scoop.mesh",
      units = "cm",
      operators = {
        { scale = scoop_scale },
        { rotate = 5 },
        { translate = offset(scoop_lift) }
      }
    }
  },
  {
    name = "colorful_sprinkles",
    material = "sprinkles",
    geometry = {
      format = "mfem",
      path = "ice_cream_sprinkles.mesh",
      units = "cm",
      operators = {
        { rotate = 15 },
        { translate = offset(sprinkle_lift) }
      }
    },
    replaces = {"ice_cream"}
  },
  {
    name = "cone",
    material = "batter",
    geometry = {
      format = "mfem",
      path = "ice_cream_cone.mesh",
      units = "cm",
      operators = {
        { rotate = -5 },
        { translate = offset(cone_lift) }
      }
    },
    does_not_replace = {"ice_cream", "sprinkles"}
  }
}
