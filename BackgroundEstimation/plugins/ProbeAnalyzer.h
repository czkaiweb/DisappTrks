#ifndef PROBE_ANALYZER

#define PROBE_ANALYZER

#include <map>
#include <string>

#include "TTree.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"


#include "DataFormats/Math/interface/deltaPhi.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"

using namespace std;

template<class T, class... Args> class ProbeAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ProbeAnalyzer (const edm::ParameterSet &);
      ~ProbeAnalyzer ();

   private:
      void analyze (const edm::Event &, const edm::EventSetup &);

      void getImage(const T &, const EBRecHitCollection &, const EERecHitCollection &, const HBHERecHitCollection &);
      const math::XYZVector getPosition(const DetId &) const;

      edm::InputTag tokenProbes_;
      edm::InputTag EBRecHitsTag_;
      edm::InputTag EERecHitsTag_;
      edm::InputTag HBHERecHitsTag_;

      edm::EDGetTokenT<vector<T> >                 probeTrackToken_;

      edm::EDGetTokenT<EBRecHitCollection>         EBRecHitsToken_;
      edm::EDGetTokenT<EERecHitCollection>         EERecHitsToken_;
      edm::EDGetTokenT<HBHERecHitCollection>       HBHERecHitsToken_;

      edm::ESHandle<CaloGeometry> caloGeometry_;

      edm::Service<TFileService> fs_;
      TTree * tree_;
      TH2D * averageImage_;

      vector<vector<double> > image_ecal, image_hcal;

      const double x_lo, x_hi, y_lo, y_hi;
      const int n_x, n_y;
 
      const int pdgId;
};

typedef ProbeAnalyzer<reco::Track> TrackProbeAnalyzer;
typedef ProbeAnalyzer<reco::GsfElectron> ElectronProbeAnlyzer;
typedef ProbeAnalyzer<reco::Muon> MuonProbeAnlyzer;

#endif
