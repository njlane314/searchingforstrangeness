#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include <vector>
#include <map>
#include <memory>

class SignatureProducer : public art::EDProducer {
public:
    explicit SignatureProducer(fhicl::ParameterSet const& pset);

    void produce(art::Event& evt) override;

private:
    recob::Track buildTrack(const std::vector<art::Ptr<recob::Hit>>& hits);
    recob::SpacePoint buildSpacePointFromHit(const art::Ptr<recob::Hit>& hit, size_t id);
    recob::Cluster buildClusterFromHits(const std::vector<art::Ptr<recob::Hit>>& hits, size_t id);
    recob::PFParticle buildPFParticle(int pdgCode, size_t id);

    std::unique_ptr<ConvolutionNetworkAlgo> _cnn_algo;
};

SignatureProducer::SignatureProducer(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
{
    _cnn_algo = std::make_unique<ConvolutionNetworkAlgo>(pset);

    produces<std::vector<recob::Track>>();
    produces<std::vector<recob::SpacePoint>>();
    produces<std::vector<recob::Cluster>>();
    produces<std::vector<recob::PFParticle>>();

    produces<art::Assns<recob::SpacePoint, recob::Hit>>();
    produces<art::Assns<recob::Cluster, recob::Hit>>();
    produces<art::Assns<recob::PFParticle, recob::Cluster>>();
}

void SignatureProducer::produce(art::Event& evt) 
{
    auto clus_coll = std::make_unique<std::vector<recob::Cluster>>();
    auto sp_coll = std::make_unique<std::vector<recob::SpacePoint>>();
    auto pfp_coll = std::make_unique<std::vector<recob::PFParticle>>();
    auto trck_coll = std::make_unique<std::vector<recob::Track>>();

    auto sp_hit_assns = std::make_unique<art::Assns<recob::SpacePoint, recob::Hit>>();
    auto clus_hit_assns = std::make_unique<art::Assns<recob::Cluster, recob::Hit>>();
    auto pfp_clus_assns = std::make_unique<art::Assns<recob::PFParticle, recob::Cluster>>();

    std::map<int, std::vector<art::Ptr<recob::Hit>>> class_hits;
    _cnn_algo->infer(evt, class_hits);  

    this->buildSpacePoints(evt, class_hits, *sp_coll, *sp_hit_assns);
    this->buildClusters(evt, class_hits, *clus_coll, *clus_hit_assns);
    this->buildPFParticles(evt, *pfp_coll, *pfp_clus_assns, *clus_coll);

    evt.put(std::move(clus_coll));
    evt.put(std::move(sp_coll));
    evt.put(std::move(pfp_coll));
    evt.put(std::move(trck_coll));

    evt.put(std::move(sp_hit_assns));
    evt.put(std::move(clus_hit_assns));
    evt.put(std::move(pfp_clus_assns));
}

DEFINE_ART_MODULE(SignatureProducer)