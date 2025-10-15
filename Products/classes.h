#include "canvas/Persistency/Common/Wrapper.h"
#include "Products/ImageProducts.h"
#include "Products/SegmentationProducts.h"
#include <vector>

template class std::vector<image::PixelImage>;
template class std::vector<image::Segmentation>;
