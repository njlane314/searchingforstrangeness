#ifndef SELECTIONAGGREGATOR_H
#define SELECTIONAGGREGATOR_H

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

class FilterAggregatorProducer;

class FilterAggregatorProducer : public art::EDProducer {
public:
    explicit FilterAggregatorProducer(fhicl::ParameterSet const& pset);
    void produce(art::Event& e) override;

private:
    std::vector<std::string> filter_labels_;
    bool use_and_;
};

FilterAggregatorProducer::FilterAggregatorProducer(fhicl::ParameterSet const& pset)
    : EDProducer{pset},
      filter_labels_(pset.get<std::vector<std::string>>("filter_labels")),
      use_and_(pset.get<bool>("use_and", true)) {
    produces<bool>();
}

void FilterAggregatorProducer::produce(art::Event& e) {
    bool result = use_and_ ? true : false;
    for (const auto& label : filter_labels_) {
        auto handle = e.getHandle<bool>(label);
        if (!handle.isValid()) {
            throw cet::exception("FilterAggregatorProducer") << "Invalid handle for " << label;
        }
        bool pass = *handle;
        if (use_and_) {
            result &= pass;
        } else {
            result |= pass;
        }
    }
    e.put(std::make_unique<bool>(result));
}

DEFINE_ART_MODULE(FilterAggregatorProducer)
#endif