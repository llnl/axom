namespace axom::bump::utilities
{
namespace detail
{
struct PathPattern
{
  static constexpr std::size_t MaxSegments = 5;

  std::array<std::string_view, MaxSegments> m_segments {};
  std::size_t m_size {0};

  constexpr PathPattern() = default;

  template <std::size_t N>
  constexpr PathPattern(const std::array<std::string_view, N>& segments) : m_size(N)
  {
    static_assert(N <= MaxSegments, "PathPattern has too many segments.");
    for(std::size_t i = 0; i < N; i++)
    {
      m_segments[i] = segments[i];
    }
  }
};

static std::string_view s_path_wildcard {"*"};

static std::array<PathPattern, 43> s_blueprint_array_patterns {
  PathPattern {
    std::array<std::string_view, 4> {"coordsets", s_path_wildcard, "values", s_path_wildcard}},
  PathPattern {std::array<std::string_view, 3> {"fields", s_path_wildcard, "values"}},
  PathPattern {std::array<std::string_view, 4> {"fields", s_path_wildcard, "values", s_path_wildcard}},
  PathPattern {std::array<std::string_view, 3> {"fields", s_path_wildcard, "matset_values"}},
  PathPattern {
    std::array<std::string_view, 4> {"fields", s_path_wildcard, "matset_values", s_path_wildcard}},
  PathPattern {std::array<std::string_view, 5> {"fields",
                                                s_path_wildcard,
                                                "matset_values",
                                                s_path_wildcard,
                                                "values"}},
  PathPattern {
    std::array<std::string_view, 4> {"topologies", s_path_wildcard, "elements", "connectivity"}},
  PathPattern {std::array<std::string_view, 4> {"topologies", s_path_wildcard, "elements", "sizes"}},
  PathPattern {
    std::array<std::string_view, 4> {"topologies", s_path_wildcard, "elements", "offsets"}},
  PathPattern {
    std::array<std::string_view, 4> {"topologies", s_path_wildcard, "elements", "shapes"}},
  PathPattern {
    std::array<std::string_view, 4> {"topologies", s_path_wildcard, "subelements", "connectivity"}},
  PathPattern {
    std::array<std::string_view, 4> {"topologies", s_path_wildcard, "subelements", "sizes"}},
  PathPattern {
    std::array<std::string_view, 4> {"topologies", s_path_wildcard, "subelements", "offsets"}},
  PathPattern {std::array<std::string_view, 3> {"matsets", s_path_wildcard, "material_ids"}},
  PathPattern {std::array<std::string_view, 3> {"matsets", s_path_wildcard, "volume_fractions"}},
  PathPattern {
    std::array<std::string_view, 4> {"matsets", s_path_wildcard, "volume_fractions", s_path_wildcard}},
  PathPattern {std::array<std::string_view, 5> {"matsets",
                                                s_path_wildcard,
                                                "volume_fractions",
                                                s_path_wildcard,
                                                "values"}},
  PathPattern {std::array<std::string_view, 5> {"matsets",
                                                s_path_wildcard,
                                                "volume_fractions",
                                                s_path_wildcard,
                                                "indices"}},
  PathPattern {std::array<std::string_view, 3> {"matsets", s_path_wildcard, "sizes"}},
  PathPattern {std::array<std::string_view, 3> {"matsets", s_path_wildcard, "offsets"}},
  PathPattern {std::array<std::string_view, 3> {"matsets", s_path_wildcard, "indices"}},
  PathPattern {
    std::array<std::string_view, 4> {"matsets", s_path_wildcard, "element_ids", s_path_wildcard}},
};

bool pathMatchesPattern(std::string_view path, const PathPattern& pattern)
{
  std::size_t cursor = 0;
  std::size_t segment_index = 0;
  std::string_view segment;

  while(nextPathSegment(path, cursor, segment))
  {
    if(segment_index >= pattern.m_size)
    {
      return false;
    }

    const auto expected = pattern.m_segments[segment_index];
    if(expected != s_path_wildcard && expected != segment)
    {
      return false;
    }

    segment_index++;
  }

  return segment_index == pattern.m_size;
}
} // end namespace detail

bool hasArrayValuedBlueprintPath(std::string_view path)
{
  for(const auto& pattern : s_blueprint_array_patterns)
  {
    if(detail::pathMatchesPattern(path, pattern))
    {
      return true;
    }
  }

  return false;
}

} // end namespace axom::bump::utilities
