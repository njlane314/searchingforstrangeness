#include <vector>

#include "canvas/Persistency/Common/Wrapper.h"

#include "Products/SparsePlaneImage.h"
#include "Products/InferenceMetrics.h"
#include "Products/InferencePredictions.h"

// Sparse plane image products
template class std::vector<image::SparsePlaneImage>;
template class art::Wrapper<image::SparsePlaneImage>;
template class art::Wrapper<std::vector<image::SparsePlaneImage>>;

// Inference metrics products
template class std::vector<image::ModelMetrics>;
template class art::Wrapper<std::vector<image::ModelMetrics>>;
template class art::Wrapper<image::InferenceMetrics>;

// Inference prediction products
template class std::vector<image::ModelPrediction>;
template class art::Wrapper<std::vector<image::ModelPrediction>>;
template class art::Wrapper<image::InferencePredictions>;
