local dim = 2

local scoop_radius = 1.1
local scoop_lift = 2.0
local sprinkle_lift = 3.0
local cone_lift = -2.0

local scoop_scale = {scoop_radius}
local scoop_offset = {0.0, scoop_lift}
local sprinkle_offset = {0.0, sprinkle_lift}
local cone_offset = {0.0, cone_lift}

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
        { translate = scoop_offset }
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
        { translate = sprinkle_offset }
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
        { translate = cone_offset }
      }
    },
    does_not_replace = {"ice_cream", "sprinkles"}
  }
}
