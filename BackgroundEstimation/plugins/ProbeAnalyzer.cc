#include "DisappTrks/BackgroundEstimation/plugins/ProbeAnalyzer.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

template<class T>
ProbeAnalyzer<T>::ProbeAnalyzer(const edm::ParameterSet &cfg):
 tokenProbes_      (cfg.getParameter<edm::InputTag> ("probeTracks")),
 EBRecHitsTag_     (cfg.getParameter<edm::InputTag> ("EBRecHits")),
 EERecHitsTag_     (cfg.getParameter<edm::InputTag> ("EERecHits")),
 HBHERecHitsTag_   (cfg.getParameter<edm::InputTag> ("HBHERecHits")),

 x_lo              (cfg.getParameter<double>        ("x_lo")),
 x_hi              (cfg.getParameter<double>        ("x_hi")),
 y_lo              (cfg.getParameter<double>        ("y_lo")),
 y_hi              (cfg.getParameter<double>        ("y_hi")),
 n_x               (cfg.getParameter<int>           ("n_x")),
 n_y               (cfg.getParameter<int>           ("n_y"))
{  
  probeTrackToken_      = consumes<vector<T>>            (tokenTags_);
 
  EBRecHitsToken_       = consumes<EBRecHitCollection>   (EBRecHitsTag_);
  EERecHitsToken_       = consumes<EERecHitCollection>   (EERecHitsTag_);
  HBHERecHitsToken_     = consumes<HBHERecHitCollection> (HBHERecHitsTag_);

  for(int i = 0; i < n_x; i++) {
    image_ecal.push_back(vector<double>(n_y, 0.0));
    image_hcal.push_back(vector<double>(n_y, 0.0));
  }
 
  tree_ = fs_->make<TTree>("tree", "tree");
  tree_->Branch("ecal", &image_ecal);
  tree_->Branch("hcal", &image_hcal);

}

template<class T>
ProbeAnalyzer<T>::~ProbeAnalyzer()
{
}

template<class T> void
ProbeAnalyzer<T>::analyze(const edm::Event &event, const edm::EventSetup &setup)
{
  edm::Handle<vector<T> > probes;
  event.getByToken (tokenTags_, probes);

  edm::Handle<EBRecHitCollection> EBRecHits;
  event.getByToken(EBRecHitsToken_, EBRecHits);
 
  edm::Handle<EERecHitCollection> EERecHits;
  event.getByToken(EERecHitsToken_, EERecHits);
 
  edm::Handle<HBHERecHitCollection> HBHERecHits;
  event.getByToken(HBHERecHitsToken_, HBHERecHits);
 
  setup.get<CaloGeometryRecord>().get(caloGeometry_);
  if (!caloGeometry_.isValid())
    throw cms::Exception("FatalError") << "Unable to find CaloGeometryRecord in event!\n"; 

  for(const auto &probe : *probes) { 
    getImage(probe, *EBRecHits, *EERecHits, *HBHERecHits);
    tree_->Fill();
  }

}

template<class T> void 
ProbeAnalyzer<T>::getImage(
  const T &probe,
  const EBRecHitCollection &EBRecHits,
  const EERecHitCollection &EERecHits,
  const HBHERecHitCollection &HBHERecHits)
{

  for(int i = 0; i < n_x; i++) {
    for(int j = 0; j < n_y; j++) {
      image_ecal[i][j] = 0.0;
      image_hcal[i][j] = 0.0;
    }
  }

  for(const auto &hit : EBRecHits) {
    math::XYZVector pos = getPosition(hit.detid());
    double dEta = probe.eta() - pos.eta();
    if(dEta < x_lo || dEta >= x_hi) continue;
    double dPhi = deltaPhi(probe, pos);
    if(dPhi < y_lo || dPhi >= y_hi) continue;

    const unsigned int ix = (dEta - x_lo) * n_x / (x_hi - x_lo);
    const unsigned int iy = (dPhi - y_lo) * n_y / (y_hi - y_lo);
    image_ecal[ix][iy] += hit.energy();
  }

  for(const auto &hit : EERecHits) {
    math::XYZVector pos = getPosition(hit.detid());
    double dEta = probe.eta() - pos.eta();
    if(dEta < x_lo || dEta >= x_hi) continue;
    double dPhi = deltaPhi(probe, pos);
    if(dPhi < y_lo || dPhi >= y_hi) continue;

    const unsigned int ix = (dEta - x_lo) * n_x / (x_hi - x_lo);
    const unsigned int iy = (dPhi - y_lo) * n_y / (y_hi - y_lo);
    image_ecal[ix][iy] += hit.energy();
  }

  for(const auto &hit : HBHERecHits) {
    math::XYZVector pos = getPosition(hit.detid());
    double dEta = probe.eta() - pos.eta();
    if(dEta < x_lo || dEta >= x_hi) continue;
    double dPhi = deltaPhi(probe, pos);
    if(dPhi < y_lo || dPhi >= y_hi) continue;

    const unsigned int ix = (dEta - x_lo) * n_x / (x_hi - x_lo);
    const unsigned int iy = (dPhi - y_lo) * n_y / (y_hi - y_lo);
    image_hcal[ix][iy] += hit.energy();
  }

}

template<class T>const math::XYZVector 
ProbeAnalyzer<T>::getPosition(const DetId& id) const
{
   if ( ! caloGeometry_.isValid() ||
        ! caloGeometry_->getSubdetectorGeometry(id) ||
        ! caloGeometry_->getSubdetectorGeometry(id)->getGeometry(id) ) {
      throw cms::Exception("FatalError") << "Failed to access geometry for DetId: " << id.rawId();
      return math::XYZVector(0,0,0);
   }
   const GlobalPoint idPosition = caloGeometry_->getSubdetectorGeometry(id)->getGeometry(id)->getPosition();
   math::XYZVector idPositionRoot(idPosition.x(), idPosition.y(), idPosition.z());
   return idPositionRoot;
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TrackProbeAnalyzer);
DEFINE_FWK_MODULE(ElectronProbeAnalyzer);
DEFINE_FWK_MODULE(MuonProbeAnalyzer);
