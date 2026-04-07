#ifndef ANALYSIS_IMAGE_PRODUCTS_COMPAT_H
#define ANALYSIS_IMAGE_PRODUCTS_COMPAT_H

#include "Products/ImageFeatures.h"

namespace image {

// Backward-compatible alias for stale dictionary dependencies and includes.
using ImageProduct = ImageFeatures;

} // namespace image

#endif // ANALYSIS_IMAGE_PRODUCTS_COMPAT_H
