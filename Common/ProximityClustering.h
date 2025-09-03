



#ifndef PROXIMITYCLUSTERING_H
#define PROXIMITYCLUSTERING_H

#include <map>

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace common {
  
  
  void MakeHitMap(const std::vector< art::Ptr<recob::Hit> >& hitlist,
		  int plane,
		  const float& _time2cm, const float& _wire2cm,
		  const float& _cellSize,
		  std::map<std::pair<int,int>, std::vector<size_t> >& _hitMap) {

    _hitMap.clear();
    
    std::pair<int,int> tmpPair;
    
    
    for (size_t h=0; h < hitlist.size(); h++){
      
      auto const& hit = hitlist.at(h);
      
      if (hit->View() != plane)
	continue;

      double t = hit->PeakTime()*_time2cm;
      double w = hit->WireID().Wire*_wire2cm;

      
      
      
      int i = int(w/_cellSize);
      int j = int(t/_cellSize);
      tmpPair = std::make_pair(i,j);
      
      
      
      if (_hitMap.find(tmpPair) == _hitMap.end()){
	std::vector<size_t> aaa = {h};
	_hitMap[tmpPair] = aaa;
      }
      else
	_hitMap[tmpPair].push_back(h);
    }

    return;
  }

  
  bool TimeOverlap(const art::Ptr<recob::Hit>& h1, const art::Ptr<recob::Hit>& h2,
		   const float& _time2cm, double& dmin) {

    auto T1 = h1->PeakTime() * _time2cm; 
    auto T2 = h2->PeakTime() * _time2cm;
    auto W1 = h1->RMS() * _time2cm;
    auto W2 = h2->RMS() * _time2cm;
    
    double d = dmin;
    
    if (T1 > T2) {
      
      if ( (T2+W2) > (T1-W1) ) return true;
      
      d = (T1-W1) - (T2+W2);
      if (d < dmin) dmin = d;
      
    }
    
    else {
      
      if ( (T1+W1) > (T2-W2) ) return true;
      
      d = (T2-W2) - (T1+W1);
      if (d < dmin) dmin = d;
      
    }
    
    return false;
  }



  
  
  bool HitsCompatible(const art::Ptr<recob::Hit>& h1, const art::Ptr<recob::Hit>& h2,
		      const float& _time2cm, const float& _wire2cm, const float& _radius) {

    if (h1->WireID().Plane != h2->WireID().Plane)
      return false;
    
    double dt = ( h1->PeakTime() - h2->PeakTime() ) * _time2cm;
    
    if (TimeOverlap(h1,h2,_time2cm,dt) == true)
      dt = 0;
    double dw = fabs(((double)h1->Channel()-(double)h2->Channel())*_wire2cm);
    if (dw >  0.3) dw -= 0.3;
    
    double d = dt*dt + dw*dw;

    if (d > (_radius*_radius))
      return false;

    return true;
  }

  
  
  void getNeighboringHits(const std::pair<int,int>& pair, std::vector<size_t>& hitIndices,
			  std::map<std::pair<int,int>, std::vector<size_t> >& _hitMap) {

    auto const& i       = pair.first;
    
    auto const& j       = pair.second;

    
    
    
    
    if (_hitMap.find(std::make_pair(i,j)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i,j)])
	hitIndices.push_back(h);
    }

    
    
    
    
    
    if (_hitMap.find(std::make_pair(i-1,j)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i-1,j)])
	hitIndices.push_back(h);
    }
    
    
    
    
    if (_hitMap.find(std::make_pair(i,j-1)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i,j-1)])
	hitIndices.push_back(h);
    }
    
    
    
    
    if ( _hitMap.find(std::make_pair(i-1,j-1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i-1,j-1)])
	hitIndices.push_back(h);
    }
    
    
    
    
    if ( _hitMap.find(std::make_pair(i,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i,j+1)])
	hitIndices.push_back(h);
    }
    
    
    
    
    if ( _hitMap.find(std::make_pair(i+1,j)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j)])
	hitIndices.push_back(h);
    }
    
    
    
    
    if ( _hitMap.find(std::make_pair(i+1,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j+1)])
	hitIndices.push_back(h);
    }
    
    
    
    
    if ( _hitMap.find(std::make_pair(i-1,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i-1,j+1)])
	hitIndices.push_back(h);
    }
    
    
    
    
    if ( _hitMap.find(std::make_pair(i+1,j-1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j-1)])
	hitIndices.push_back(h);
    }
  }

  
  
  
  
  bool cluster(const std::vector< art::Ptr<recob::Hit> >& hit_ptr_v,
	       std::vector<std::vector<unsigned int> >& _out_cluster_vector,
	       const float& cellSize, const float& radius) {

    if (hit_ptr_v.size() == 0)
      return false;
    
    
    double _cellSize = cellSize;
    
    
    double _radius = radius;
    
    
    double _wire2cm, _time2cm;
    
    
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    _wire2cm = geom->WirePitch(0,0,0);
    _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
    
    
    
    std::map<std::pair<int,int>, std::vector<size_t> > _hitMap;

    
    
    
    std::map<size_t, size_t> _clusterMap;
    
    
    std::map<size_t,std::vector<size_t> > _clusters;

    
    size_t maxClusterID = 0;

    for (int pl=0; pl < 3; pl++){
      
      
      MakeHitMap(hit_ptr_v,pl,_time2cm, _wire2cm, _cellSize, _hitMap);
      
      
      std::map<std::pair<int,int>, std::vector<size_t> >::iterator it;
      
      
      for (it = _hitMap.begin(); it != _hitMap.end(); it++){
	
	
	auto const& pair = it->first;
	
	
	
	
	
	
	
	std::vector<size_t> cellhits = it->second;

	std::vector<size_t> neighborhits;
	getNeighboringHits(pair,neighborhits, _hitMap);

	for (size_t h1=0; h1 < cellhits.size(); h1++){

	  
	  
	  auto const& hit1 = cellhits[h1];
	  
	  bool matched = false;
	  
	  for (size_t h2=0; h2 < neighborhits.size(); h2++){
	    auto const& hit2 = neighborhits[h2];
	    if (hit1 == hit2) continue;
	    
	    bool compat = HitsCompatible(hit_ptr_v.at(hit1),
					 hit_ptr_v.at(hit2),
					 _time2cm, _wire2cm, _radius);
	    
	    if (compat){
	      matched = true;
	      
	      if ( (_clusterMap.find(hit1) != _clusterMap.end()) and
		   (_clusterMap.find(hit2) != _clusterMap.end()) ){
		
		
		if (_clusterMap[hit1] != _clusterMap[hit2]){
		  auto idx1 = _clusterMap[hit1];
		  auto idx2 = _clusterMap[hit2];
		  
		  auto hits1 = _clusters[idx1];
		  auto hits2 = _clusters[idx2];
		  
		  for (auto h : hits2){
		    hits1.push_back(h);
		    
		    _clusterMap[h] = idx1;
		  }
		  _clusters[idx1] = hits1;
		  
		  _clusters.erase(idx2);
		}
	      }
	      
	      
	      else if ( (_clusterMap.find(hit2) != _clusterMap.end()) and
			(_clusterMap.find(hit1) == _clusterMap.end()) ){
		auto clusIdx = _clusterMap[hit2];
		_clusterMap[hit1] = clusIdx;
		_clusters[clusIdx].push_back(hit1);
	      }
	      
	      else if ( (_clusterMap.find(hit1) != _clusterMap.end()) and
			(_clusterMap.find(hit2) == _clusterMap.end()) ){
		auto clusIdx = _clusterMap[hit1];
		_clusterMap[hit2] = clusIdx;
		_clusters[clusIdx].push_back(hit2);
	      }
	      
	      else{
		
		_clusterMap[hit1] = maxClusterID;
		_clusterMap[hit2] = maxClusterID;
		std::vector<size_t> cl = {hit1,hit2};
		_clusters[maxClusterID] = cl;
		maxClusterID += 1;
	      }
	    }
	  }
	  
	  if (matched == false){
	    _clusterMap[hit1] = maxClusterID;
	    _clusters[maxClusterID] = {hit1};
	    maxClusterID += 1;
	  }
	}
      }

    }

    
    for (auto it = _clusters.begin(); it != _clusters.end(); it++){
      auto indices = it->second;
      
      if (indices.size() >= 1){
	std::vector<unsigned int> clus;
	for (auto idx : indices)
	  clus.push_back(idx);
	_out_cluster_vector.push_back(clus);
      }
    }
    
    return true;
  }

      
}
#endif

/ 