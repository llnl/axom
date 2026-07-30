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
