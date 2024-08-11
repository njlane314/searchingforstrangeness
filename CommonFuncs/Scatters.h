#ifndef SCATTERSFUNCS_H
#define SCATTERSFUNCS_H

namespace common
{
    std::vector< art::Ptr<simb::MCParticle>> GetDaughters(const art::Ptr<simb::MCParticle> &particle, const std::map<int, art::Ptr<simb::MCParticle> > &mcParticleMap)
    {
        std::vector< art::Ptr<simb::MCParticle>> daughters;
        for (int i = 0; i < particle->NumberDaughters(); ++i)
        {
            const auto daughterIter = mcParticleMap.find(particle->Daughter(i));
            if (daughterIter != mcParticleMap.end()) daughters.push_back(daughterIter->second);
        }
        return daughters;
    }

    void GetNScatters(const art::ValidHandle<std::vector<simb::MCParticle>> &mcp_h, const art::Ptr<simb::MCParticle> &mcParticle, art::Ptr<simb::MCParticle> &mcScatteredParticle, unsigned int &nElastic, unsigned int &nInelastic)
    {
        mcScatteredParticle = mcParticle;

        std::map<int, art::Ptr<simb::MCParticle> > mcParticleMap;
        for (size_t d = 0; d < mcp_h->size(); d++)
        {
            const art::Ptr<simb::MCParticle> mcParticle(mcp_h, d);
            if (!mcParticleMap.emplace(mcParticle->TrackId(), mcParticle).second)
                throw cet::exception("::GetNScatters") << " - Found repeated MCParticle with TrackId = " << mcParticle->TrackId() << "." << std::endl;
        }

        art::Ptr<simb::MCParticle> finalStateParticle;
        bool foundInelasticScatter = false;
        for (const auto &daughter : GetDaughters(mcParticle, mcParticleMap))
        {
            const auto& process = daughter->Process();

            if (process == "hadElastic") 
            {
                nElastic++;
            }
            else if (process.find("Inelastic") != std::string::npos)
            {
                if (daughter->PdgCode() != mcParticle->PdgCode()) continue;

                if (foundInelasticScatter) 
                {
                    foundInelasticScatter = false;

                    break;
                }
                
                finalStateParticle = daughter;
                foundInelasticScatter = true;
            }
        }

        if (foundInelasticScatter)
        {
            nInelastic++;
            GetNScatters(mcp_h, finalStateParticle, mcScatteredParticle, nElastic, nInelastic);
        }
    }

    std::string GetEndState(const art::Ptr<simb::MCParticle> &particle, const art::ValidHandle<std::vector<simb::MCParticle>> &mcp_h)
    {
        std::string type = "Other";
        bool hasPi0 = false;
        bool hasDecayMuon = false;
        bool hasDecayMuonNeutrino = false;

        std::map<int, art::Ptr<simb::MCParticle>> mcParticleMap;
        for (size_t d = 0; d < mcp_h->size(); d++)
        {
            const art::Ptr<simb::MCParticle> mcParticle(mcp_h, d);
            mcParticleMap[mcParticle->TrackId()] = mcParticle;
        }

        art::Ptr<simb::MCParticle> scatteredParticle = particle;
        unsigned int nElastic = 0;
        unsigned int nInelastic = 0;
        GetNScatters(mcp_h, particle, scatteredParticle, nElastic, nInelastic);

        std::vector<art::Ptr<simb::MCParticle>> products;
        for (const auto &daughter : GetDaughters(scatteredParticle, mcParticleMap))
        {
            const auto process = daughter->Process();

            if (daughter->PdgCode() == 11 && process == "hIoni")
                continue;

            if (process == "hadElastic")
                continue;

            products.push_back(daughter);

            if (daughter->PdgCode() == 111 && (process == "pi+Inelastic" || process == "pi-Inelastic"))
                hasPi0 = true;

            if (daughter->PdgCode() == -13 && process == "Decay")
                hasDecayMuon = true;

            if (daughter->PdgCode() == 14 && process == "Decay")
                hasDecayMuonNeutrino = true;
        }

        if (products.empty())
        {
            type = "None";
        }
        else if (hasDecayMuon && hasDecayMuonNeutrino && products.size() == 2)
        {
            type = "DecayToMuon";
        }
        else if (scatteredParticle->EndProcess() == "pi+Inelastic" || scatteredParticle->EndProcess() == "pi-Inelastic")
        {
            type = hasPi0 ? "Pi0ChargeExchange" : "InelasticAbsorption";
        }
        else
        {
            type = "Other";
        }

        return type;
    }
}

#endif