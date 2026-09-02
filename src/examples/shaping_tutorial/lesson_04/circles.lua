local outer_radius = 5.0
local inner_radius_mm = 25.0

dimensions = 2

shapes = {
  {
    name = "background",
    material = "void",
    geometry = {
      format = "none"
    }
  },
  {
    name = "outer_shell",
    material = "steel",
    geometry = {
      format = "mfem",
      path = "unit_circle.mesh",
      units = "cm",
      operators = {
        {
          scale = function()
            return {outer_radius}
          end
        }
      }
    }
  },
  {
    name = "inner_ball",
    material = "void",
    geometry = {
      format = "mfem",
      path = "unit_circle.mesh",
      start_units = "mm",
      end_units = "cm",
      operators = {
        {
          scale = function()
            return {inner_radius_mm}
          end
        },
        { convert_units_to = "cm" }
      }
    }
  }
}
