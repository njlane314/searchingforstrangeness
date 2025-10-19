#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "Products/ImageProducts.h"        // PlaneImage (your existing definition)
#include "Products/VertexDisplacement.h"   // our single product

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace {

struct IAOKHeader {
  char     magic[4];
  uint32_t version;
  uint32_t K;
  uint32_t segW;
  uint32_t segH;
  uint32_t has_conf;
  uint64_t cls_bytes;
  uint64_t seg_bytes;
  uint64_t conf_bytes;
};

inline std::string jpath(std::string a, const std::string& b) {
  if (!a.empty() && a.back() != '/') a.push_back('/');
  return a + b;
}

inline void write_chw_f32(const std::string& path,
                          const std::vector<float>& U,
                          const std::vector<float>& V,
                          const std::vector<float>& W) {
  std::ofstream ofs(path, std::ios::binary | std::ios::trunc);
  if (!ofs) throw cet::exception("VDI") << "Cannot open " << path;
  ofs.write(reinterpret_cast<const char*>(U.data()), sizeof(float)*U.size());
  ofs.write(reinterpret_cast<const char*>(V.data()), sizeof(float)*V.size());
  ofs.write(reinterpret_cast<const char*>(W.data()), sizeof(float)*W.size());
  if (!ofs) throw cet::exception("VDI") << "Short write to " << path;
}

inline bool exists(const std::string& p) { struct stat st{}; return ::stat(p.c_str(), &st) == 0; }

// ---- Beam utilities (fixed NuMI @ MicroBooNE) ----
struct Beam3D { double bx, by, bz; };

inline Beam3D beam_from_angles_deg(double theta_xz_deg, double theta_yz_deg) {
  const double tx = std::tan(theta_xz_deg * M_PI/180.0);
  const double ty = std::tan(theta_yz_deg * M_PI/180.0);
  const double nx = tx, ny = ty, nz = 1.0;
  const double n = std::sqrt(nx*nx + ny*ny + nz*nz);
  return {nx/n, ny/n, nz/n};
}

// Column mapping: U: z*c - y*s; V: z*c + y*s; W: z ; Row is x
inline std::pair<double,double>
beam2d_pixels_for_view(int view, double pixel_w, double pixel_h, const Beam3D& b) {
  const double c60 = 0.5, s60 = std::sqrt(3.0)/2.0;
  double col_cm = 0.0;
  if      (view == static_cast<int>(geo::kU)) col_cm =  b.bz*c60 - b.by*s60;
  else if (view == static_cast<int>(geo::kV)) col_cm =  b.bz*c60 + b.by*s60;
  else                                        col_cm =  b.bz; // W or fallback
  const double row_cm = b.bx;
  return { col_cm / (pixel_w == 0 ? 1.0 : pixel_w),
           row_cm / (pixel_h == 0 ? 1.0 : pixel_h) };
}

// Read flat float32 file
inline std::vector<float> read_f32(const std::string& path, size_t n) {
  std::vector<float> v(n, 0.f);
  std::ifstream f(path, std::ios::binary);
  if (f) f.read(reinterpret_cast<char*>(v.data()),
                static_cast<std::streamsize>(n*sizeof(float)));
  return v;
}

} // namespace

namespace image {

class VertexDisplacementInference : public art::EDProducer {
public:
  explicit VertexDisplacementInference(fhicl::ParameterSet const& p)
  : planesTag_{p.get<std::string>("PlanesTag")}
  , assetsBaseDir_{p.get<std::string>("AssetsBaseDir")}
  , wrapper_{p.get<std::string>("InferenceWrapper", jpath(assetsBaseDir_,"scripts/inference_wrapper.sh"))}
  , container_{p.get<std::string>("ContainerImage",
      "/cvmfs/uboone.opensciencegrid.org/containers/lantern_v2_me_06_03_prod")}
  , apptainer_{p.get<std::string>("ApptainerPath","apptainer")} {
    produces<image::VertexDisplacement>();
  }

  void produce(art::Event& e) override {
    auto planes = e.getValidHandle<std::vector<image::PlaneImage>>(planesTag_);
    if (planes->size() < 3)
      throw cet::exception("VDI") << "Need 3 planes (U,V,W), got " << planes->size();

    // Identify U,V,W by view (fallback to first 3)
    auto const* U=&planes->at(0);
    auto const* V=&planes->at(1);
    auto const* W=&planes->at(2);
    for (auto const& P : *planes) {
      if      (P.view == static_cast<int>(geo::kU)) U = &P;
      else if (P.view == static_cast<int>(geo::kV)) V = &P;
      else if (P.view == static_cast<int>(geo::kW)) W = &P;
    }

    const uint32_t Wpix = U->width, Hpix = U->height;
    auto chk = [&](auto const& P){
      if (P.width!=Wpix || P.height!=Hpix) throw cet::exception("VDI")<<"Plane size mismatch";
      if (P.adc.size() != size_t(Wpix)*Hpix) throw cet::exception("VDI")<<"ADC size mismatch";
    };
    chk(*V); chk(*W);

    // Clamp negatives to zero (model uses occupancy-like features)
    auto clamp0 = [](std::vector<float> v){ for(auto& x: v) if (x<0) x=0; return v; };
    auto u = clamp0(U->adc), v = clamp0(V->adc), w = clamp0(W->adc);

    // Fixed NuMI-at-MicroBooNE beam; not overridable
    constexpr double kThetaXZ = +8.0;   // degrees
    constexpr double kThetaYZ = -3.0;   // degrees
    auto b3d = beam_from_angles_deg(kThetaXZ, kThetaYZ);

    auto bu = beam2d_pixels_for_view(U->view, U->pixel_w, U->pixel_h, b3d);
    auto bv = beam2d_pixels_for_view(V->view, V->pixel_w, V->pixel_h, b3d);
    auto bw = beam2d_pixels_for_view(W->view, W->pixel_w, W->pixel_h, b3d);

    // Scratch files
    std::string scratch = ".";
    if (const char* s = std::getenv("_CONDOR_SCRATCH_DIR")) scratch = s;
    char buf[4096]; std::string abs = realpath(scratch.c_str(), buf) ? std::string(buf) : scratch;
    auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    int pid = ::getpid();
    std::string base = "vd_" + std::to_string(pid) + "_" + std::to_string(now);

    std::string req   = jpath(abs, base + "_planes.bin");
    std::string out   = jpath(abs, base + "_out.bin");
    std::string meta  = out + ".meta";
    std::string sout  = jpath(abs, base + "_stdout.txt");
    std::string serr  = jpath(abs, base + "_stderr.txt");

    auto t0 = std::chrono::steady_clock::now();
    write_chw_f32(req, u, v, w);
    auto t1 = std::chrono::steady_clock::now();

    // Build apptainer command (no extra flags; only what we use)
    std::ostringstream cmd;
    if (!container_.empty()) {
      std::ostringstream binds; binds << "/cvmfs," << abs;
      if (!assetsBaseDir_.empty()) binds << "," << assetsBaseDir_;
      cmd << apptainer_ << " exec --cleanenv --bind " << binds.str()
          << " " << container_ << " /bin/bash " << wrapper_;
    } else {
      cmd << "/bin/bash " << wrapper_;
    }
    cmd << " --in " << req
        << " --out " << out
        << " --metrics " << meta
        << " --W " << Wpix
        << " --H " << Hpix
        << " --arch pvdv2d_me"
        << " --beam2d-u " << bu.first << " " << bu.second
        << " --beam2d-v " << bv.first << " " << bv.second
        << " --beam2d-w " << bw.first << " " << bw.second
        << " > " << sout << " 2> " << serr;

    mf::LogInfo("VDI") << "Executing: " << cmd.str();

    auto start = std::chrono::steady_clock::now();
    int code = std::system(cmd.str().c_str());
    auto end = std::chrono::steady_clock::now();
    auto t2 = std::chrono::steady_clock::now();

    if (code != 0 || !exists(out)) {
      std::ifstream es(serr);
      std::string err((std::istreambuf_iterator<char>(es)), {});
      throw cet::exception("VDI") << "External run failed, code=" << code << "\n" << err;
    }

    // Read IAOK (classification only)
    std::ifstream ifs(out, std::ios::binary);
    IAOKHeader h{};
    ifs.read(reinterpret_cast<char*>(&h), sizeof(h));
    if (!ifs || std::string(h.magic, h.magic+4)!="IAOK" || h.version!=1)
      throw cet::exception("VDI") << "Bad IAOK header";
    float response = 0.f;
    if (h.K == 1 && h.cls_bytes == sizeof(float)) {
      ifs.read(reinterpret_cast<char*>(&response), sizeof(float));
    }

    auto t3 = std::chrono::steady_clock::now();

    // Parse metrics (timings + features sidecar pointer)
    VDPerf perf{};
    std::string features_path; int feat_dim = 0;
    {
      std::ifstream m(meta);
      std::string line;
      auto num = [](const std::string& s, double& d)->bool{
        char* ep=nullptr; d = std::strtod(s.c_str(), &ep); return ep!=s.c_str();
      };
      while (std::getline(m,line)) {
        auto eq = line.find('='); if (eq==std::string::npos) continue;
        std::string k = line.substr(0,eq), v = line.substr(eq+1);
        double d=0.0;
        if      (k=="t_total_ms" && num(v,d)) perf.t_child_total_ms = d;
        else if (k=="t_setup_ms" && num(v,d)) perf.t_child_setup_ms = d;
        else if (k=="t_infer_ms" && num(v,d)) perf.t_child_infer_ms = d;
        else if (k=="t_post_ms"  && num(v,d)) perf.t_child_post_ms  = d;
        else if (k=="max_rss_mb" && num(v,d)) perf.child_max_rss_mb = d;
        else if (k=="cuda_mem_mb"&& num(v,d)) perf.child_cuda_mem_mb= d;
        else if (k=="features_path")          features_path = v;
        else if (k=="feat_dim" && num(v,d))   feat_dim = static_cast<int>(d);
      }
    }

    perf.t_write_req_ms  = std::chrono::duration<float,std::milli>(t1 - t0).count();
    perf.t_exec_total_ms = std::chrono::duration<float,std::milli>(end - start).count();
    perf.t_read_resp_ms  = std::chrono::duration<float,std::milli>(t3 - t2).count();

    // Build product
    auto outprod = std::make_unique<image::VertexDisplacement>();
    outprod->algorithm = "pvdv2d_me";
    outprod->response  = response;
    outprod->theta_xz_deg = static_cast<float>(kThetaXZ);
    outprod->theta_yz_deg = static_cast<float>(kThetaYZ);
    outprod->perf = perf;

    auto fill_plane = [&](const image::PlaneImage* P,
                          const std::pair<double,double>& b2d) -> image::VDPlaneSummary {
      image::VDPlaneSummary s;
      s.view   = P->view;
      s.width  = P->width;
      s.height = P->height;
      s.beam_col_pix = static_cast<float>(b2d.first);
      s.beam_row_pix = static_cast<float>(b2d.second);
      return s;
    };
    outprod->U = fill_plane(U, bu);
    outprod->V = fill_plane(V, bv);
    outprod->W = fill_plane(W, bw);

    // Parse features sidecar (PV/R_pix/peaks per plane)
    // Layout per plane (34 floats):
    // [ pv_col, pv_row, R_pix, K,
    //   dv0_col, dv0_row, dv0_s, ..., up to 10 triplets ]
    if (!features_path.empty() && feat_dim == 34*3) {
      std::vector<float> F = read_f32(features_path, static_cast<size_t>(feat_dim));
      auto parse_plane = [&](int i, image::VDPlaneSummary& S){
        const int off = 34*i;
        const float pv_c = F[off+0], pv_r = F[off+1];
        const float Rpix = F[off+2];
        int K = static_cast<int>(std::round(F[off+3]));
        if (K < 0) K = 0; if (K > 10) K = 10;
        S.pv_col = static_cast<std::uint16_t>(std::lround(pv_c));
        S.pv_row = static_cast<std::uint16_t>(std::lround(pv_r));
        S.R_pix  = Rpix;
        S.peaks.clear(); S.peaks.reserve(static_cast<size_t>(K));
        float top = 0.f;
        for (int j=0; j<K; ++j) {
          const int b = off + 4 + 3*j;
          VDPeak p;
          p.col   = static_cast<std::uint16_t>(std::lround(F[b+0]));
          p.row   = static_cast<std::uint16_t>(std::lround(F[b+1]));
          p.score = F[b+2];
          if (j==0) top = p.score; // already sorted desc in python
          S.peaks.push_back(p);
        }
        S.top_score = top;
      };
      parse_plane(0, outprod->U);
      parse_plane(1, outprod->V);
      parse_plane(2, outprod->W);
    }

    e.put(std::move(outprod));

    // Cleanup everything (including sidecar)
    std::remove(req.c_str());
    std::remove(out.c_str());
    std::remove(meta.c_str());
    std::remove(sout.c_str());
    std::remove(serr.c_str());
    if (!features_path.empty()) std::remove(features_path.c_str());
  }

private:
  art::InputTag planesTag_;
  std::string   assetsBaseDir_;
  std::string   wrapper_;
  std::string   container_;
  std::string   apptainer_;
};

DEFINE_ART_MODULE(image::VertexDisplacementInference)

} // namespace image
