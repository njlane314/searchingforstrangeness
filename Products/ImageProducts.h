#ifndef ANALYSIS_IMAGE_PRODUCTS_COMPAT_H
#define ANALYSIS_IMAGE_PRODUCTS_COMPAT_H

#include "Products/SparsePlaneImage.h"

namespace image {

// Backward-compatible alias for stale dictionary dependencies and includes.
using ImageProduct = SparsePlaneImage;

} // namespace image

#endif // ANALYSIS_IMAGE_PRODUCTS_COMPAT_H
