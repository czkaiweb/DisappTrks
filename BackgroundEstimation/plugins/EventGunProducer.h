#ifndef EVENT_TAG_AND_PROBE_PRODUCER
#define EVENT_TAG_AND_PROBE_PRODUCER

#include <sstream>
#include "TLorentzVector.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "OSUT3Analysis/AnaTools/interface/EventVariableProducer.h"
#include "OSUT3Analysis/AnaTools/interface/DataFormat.h"
#include "OSUT3Analysis/AnaTools/interface/CommonUtils.h"
#include "OSUT3Analysis/Collections/interface/Electron.h"
#include "OSUT3Analysis/Collections/interface/Muon.h"
#include "OSUT3Analysis/Collections/interface/Tau.h"
#include "OSUT3Analysis/Collections/interface/Track.h"
#include "OSUT3Analysis/Collections/interface/GenMatchable.h"

template<class T, class... Args> class EventGunProducer : public EventVariableProducer
{
  public:
    EventGunProducer (const edm::ParameterSet &);
    ~EventGunProducer ();

  private:
    void AddVariables(const edm::Event &);
    edm::EDGetTokenT<vector<T> > tokenProbes_;
    edm::EDGetTokenT<vector<TYPE(muons)> > tokenMuons_;
    edm::EDGetTokenT<vector<reco::GenParticle> > genParticlesToken_;
    bool doFilter_;
    bool doSSFilter_;
    bool doJetFilter_;

    const string tagCollectionParameter () const;
    bool passesVeto (const osu::Track &) const;
    bool passesOuterHits (const osu::Track &) const;
    template<class T0> const double getTrackIsolation (const T0 &, const vector<T0> &, const double, const double = 1.0e-12) const;
    bool jetMatchedToMuon (const pat::Jet &, const vector<pat::PackedCandidate> &) const;
};

typedef EventGunProducer<osu::Track,TYPE(muons) > EventMuonPGunProducer;
#endif