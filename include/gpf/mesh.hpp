#pragma once

#include "gpf/handles.hpp"
#include "gpf/ids.hpp"
#include "gpf/manifold_mesh.hpp"
#include "gpf/property.hpp"
#include "gpf/surface_mesh.hpp"
#include "gpf/utils.hpp"

namespace gpf {
struct Empty {
  [[nodiscard]] constexpr bool operator==(const Empty&) const = default;
};
}
