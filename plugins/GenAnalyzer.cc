#include "GenAnalyzer.h"


GenAnalyzer::GenAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    GenToken(CColl.consumes<GenEventInfoProduct>(PSet.getParameter<edm::InputTag>("genProduct"))),
    LheToken(CColl.consumes<LHEEventProduct>(PSet.getParameter<edm::InputTag>("lheProduct"))),
    GenParticlesToken(CColl.consumes<std::vector<reco::GenParticle> >(PSet.getParameter<edm::InputTag>("genParticles")))
{
    /*
    Sample=sample;
    isDYFile=false;
    Sum=Num=NULL;
    // Read root file
    DYFile=new TFile("data/DYWeight.root", "READ");
    DYFile->cd();
    if(!DYFile->IsZombie()) {
      Num=(TH3F*)DYFile->Get(Sample.c_str());
      Sum=(TH3F*)DYFile->Get("Sum");
      if(Sum && Num && Sum->GetEntries()>0 && Num->GetEntries()>0) {
        npMax=Sum->GetXaxis()->GetBinCenter(Sum->GetNbinsX());
        ptMax=Sum->GetYaxis()->GetBinCenter(Sum->GetNbinsY());
        htMax=Sum->GetZaxis()->GetBinCenter(Sum->GetNbinsZ());
        isDYFile=true;
      }
      else std::cout << " - GenAnalyzer Warning: Drell-Yan initialization failed, check rootfile" << std::endl;
      
      
    }
    else std::cout << " - GenAnalyzer Warning: No Drell-Yan File" << std::endl;
    */
    
    // PU reweighting
//    LumiWeights=new edm::LumiReWeighting("data/MC_True.root", "data/Prod6.root", "S10", "pileup");
    
    std::cout << " --- GenAnalyzer initialization ---" << std::endl;
    std::cout << std::endl;
}

GenAnalyzer::~GenAnalyzer() {
//    delete LumiWeights;
    /*
    DYFile->Close();
    */
}

// ---------- GEN WEIGHTS ----------

std::map<std::string, float> GenAnalyzer::FillWeightsMap(const edm::Event& iEvent) {
    std::map<std::string, float> Weights;
    Weights["event"] = 1.;
    if(iEvent.isRealData()) return Weights;
    // Declare and open collection
    edm::Handle<GenEventInfoProduct> GenEventCollection;
    iEvent.getByToken(GenToken, GenEventCollection);
    // Declare and open collection
    edm::Handle<LHEEventProduct> LheEventCollection;
    iEvent.getByToken(LheToken, LheEventCollection);
    
    Weights["event"] = fabs(LheEventCollection.product()->originalXWGTUP()) / LheEventCollection.product()->originalXWGTUP();
    return Weights;
}


// ---------- GEN PARTICLES ----------

std::vector<reco::GenParticle> GenAnalyzer::FillGenVector(const edm::Event& iEvent) {
    std::vector<reco::GenParticle> Vect;
    if(iEvent.isRealData()) return Vect;
    // Declare and open collection
    edm::Handle<std::vector<reco::GenParticle> > GenCollection;
    iEvent.getByToken(GenParticlesToken, GenCollection);
    // Loop on Gen Particles collection
    for(std::vector<reco::GenParticle>::const_iterator it = GenCollection->begin(); it != GenCollection->end(); ++it) {Vect.push_back(*it); // Fill vector
//        std::cout << it->pdgId() << "  " << it->status() << "  " << it->pt() << "  " << it->eta() << "  " << it->phi() << "  " << it->mass() << std::endl;
//        if(it->numberOfDaughters()>0) std::cout << "  " << it->daughter(0)->pdgId() << "  " << it->daughter(0)->status() << "  " << it->daughter(0)->pt() << std::endl;
//        if(it->numberOfDaughters()>1) std::cout << "  " << it->daughter(1)->pdgId() << "  " << it->daughter(1)->status() << "  " << it->daughter(1)->pt() << std::endl;
    }
//    std::cout << "\n\n\n" << std::endl;
    return Vect;
}

reco::Candidate* GenAnalyzer::FindGenParticle(std::vector<reco::GenParticle>& Vect, int pdg) {
    for(unsigned int i = 0; i < Vect.size(); i++) {
        if(Vect[i].pdgId() == pdg) return FindLastDaughter(dynamic_cast<reco::Candidate*>(&Vect[i]));
    }
    return NULL;
}

// Recursive function to find the last particle in the chain before decay: e.g. 23 -> 23 -> *23* -> 13 -13
reco::Candidate* GenAnalyzer::FindLastDaughter(reco::Candidate* p) {
    if(p->numberOfDaughters() <= 0 || !p->daughter(0)) return p;
    if(p->daughter(0)->pdgId() != p->pdgId()) return p;
    return FindLastDaughter(p->daughter(0));
}





// ---------- PILEUP ----------

float GenAnalyzer::GetPUWeight(const edm::Event& iEvent) {
  //  int nPT(0);
  //  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchool2012PileupReweighting
  //  if(!iEvent.isRealData()) {
  //    edm::Handle<std::vector<PileupSummaryInfo> > PUInfo;
  //    iEvent.getByLabel(edm::InputTag("PileupSummaryInfo"), PUInfo);
  //    for(std::vector<PileupSummaryInfo>::const_iterator pvi=PUInfo->begin(), pvn=PUInfo->end(); pvi!=pvn; ++pvi) {
  //      if(pvi->getBunchCrossing()==0) nPT=pvi->getTrueNumInteractions(); // getPU_NumInteractions();
  //    }
  //    return LumiWeights->weight( nPT );
  //  }
    return 1.;
}


// ---------- PDF ----------

std::pair<float, float> GenAnalyzer::GetQ2Weight(const edm::Event& iEvent) {
  //  float Q, id1, id2, x1, x2;
  //  // Open LHE event
  //  edm::Handle<LHEEventProduct> lheProduct;
  //  iEvent.getByLabel(edm::InputTag("source"), lheProduct);
  //  // Access event info from LHE
  //  if(lheProduct.isValid()) {
  //    const lhef::HEPEUP hepeup=lheProduct->hepeup();
  //    // PDF
  //    Q   = hepeup.SCALUP;
  //    // id of the particle: 0 is for gluons
  //    id1 = hepeup.IDUP[0]==21 ? 0 : hepeup.IDUP[0];
  //    id2 = hepeup.IDUP[1]==21 ? 0 : hepeup.IDUP[1];
  //    x1  = fabs(hepeup.PUP[0][2]/6500.);
  //    x2  = fabs(hepeup.PUP[1][2]/6500.);
  //    //xfx1 = hepeup.XPDWUP.first;
  //    //xfx2 = hepeup.XPDWUP.second;
  //  }
  //  // Access event info from Gen if LHE info are not available
  //  else {
  //    edm::Handle<GenEventInfoProduct> genProduct;
  //    iEvent.getByLabel(edm::InputTag("generator"), genProduct);
  //    Q   = genProduct->pdf()->scalePDF;
  //    id1 = genProduct->pdf()->id.first;
  //    id2 = genProduct->pdf()->id.second;
  //    x1  = genProduct->pdf()->x.first;
  //    x2  = genProduct->pdf()->x.second;
  //    //pdf1 = genProduct->pdf()->xPDF.first;
  //    //pdf2 = genProduct->pdf()->xPDF.second;
  //  }
  //  //LHAPDF::usePDFMember(1, 0);
  //  float pdf1 = LHAPDF::xfx(1, x1, Q, id1)/x1;
  //  float pdf2 = LHAPDF::xfx(1, x2, Q, id2)/x2;
  //  // Q2 scale
  //  float newqcd1_up = LHAPDF::xfx(1, x1, Q*2, id1)/x1;
  //  float newqcd2_up = LHAPDF::xfx(1, x2, Q*2, id2)/x2;
  //  float newqcd1_down = LHAPDF::xfx(1, x1, Q/2, id1)/x1;
  //  float newqcd2_down = LHAPDF::xfx(1, x2, Q/2, id2)/x2;
  //  float Q2ScaleWeightUp = (newqcd1_up/pdf1)*(newqcd2_up/pdf2);
  //  float Q2ScaleWeightDown = (newqcd1_down/pdf1)*(newqcd2_down/pdf2);
    float Q2ScaleWeightUp = 1.;
    float Q2ScaleWeightDown = 1.;
    return std::pair<float, float>(Q2ScaleWeightUp, Q2ScaleWeightDown);
}

//float GenAnalyzer::GetPDFWeight(const edm::Event& iEvent) {
//  float Q, id1, id2, x1, x2;
//  // Open LHE event
//  edm::Handle<LHEEventProduct> lheProduct;
//  iEvent.getByLabel(edm::InputTag("source"), lheProduct);
//  // Access event info from LHE
//  if(lheProduct.isValid()) {
//    const lhef::HEPEUP hepeup=lheProduct->hepeup();
//    // PDF
//    Q   = hepeup.SCALUP;
//    // id of the particle: 0 is for gluons
//    id1 = hepeup.IDUP[0]==21 ? 0 : hepeup.IDUP[0];
//    id2 = hepeup.IDUP[1]==21 ? 0 : hepeup.IDUP[1];
//    x1  = fabs(hepeup.PUP[0][2]/6500.);
//    x2  = fabs(hepeup.PUP[1][2]/6500.);
//    //xfx1 = hepeup.XPDWUP.first;
//    //xfx2 = hepeup.XPDWUP.second;
//  }
//  // Access event info from Gen if LHE info are not available
//  else {
//    edm::Handle<GenEventInfoProduct> genProduct;
//    iEvent.getByLabel(edm::InputTag("generator"), genProduct);
//    Q   = genProduct->pdf()->scalePDF;
//    id1 = genProduct->pdf()->id.first;
//    id2 = genProduct->pdf()->id.second;
//    x1  = genProduct->pdf()->x.first;
//    x2  = genProduct->pdf()->x.second;
//    //pdf1 = genProduct->pdf()->xPDF.first;
//    //pdf2 = genProduct->pdf()->xPDF.second;
//  }
//  //LHAPDF::usePDFMember(1, 0);
//  float pdf1 = LHAPDF::xfx(1, x1, Q, id1)/x1;
//  float pdf2 = LHAPDF::xfx(1, x2, Q, id2)/x2;
//  // New PDF, if not working type <scramv1 setup lhapdffull> to enable more than one LHAPDF set
//  //LHAPDF::usePDFMember(2, 0);
//  float newpdf1 = LHAPDF::xfx(2, x1, Q, id1)/x1;
//  float newpdf2 = LHAPDF::xfx(2, x2, Q, id2)/x2;
//  return (newpdf1/pdf1)*(newpdf2/pdf2);
//}

// ---------- LHE Event ----------
/*
float GenAnalyzer::GetDYWeight(const edm::Event& iEvent) {
  if(!isDYFile) return 1.;
  int lhePartons(0), lheBPartons(0);
  float lheHT(0.), lhePtZ(0.);
  // Open LHE event if possible
  edm::Handle<LHEEventProduct> lheProduct;
  iEvent.getByLabel(edm::InputTag("source"), lheProduct);
  if(lheProduct.isValid()) return 1.;
  const lhef::HEPEUP hepeup=lheProduct->hepeup();
  
  // Find LHE Z and partons
  for(int i=0; i<hepeup.NUP; ++i) {
    int id=abs(hepeup.IDUP[i]);
    // Lab frame momentum (Px, Py, Pz, E and M in GeV) for the particle entries in this event
    reco::Candidate::LorentzVector P4(hepeup.PUP[i][0], hepeup.PUP[i][1], hepeup.PUP[i][2], hepeup.PUP[i][3]);
    if(hepeup.ISTUP[i]==1 && (id<6 || id==21)) {
      lheHT+=P4.pt();
      lhePartons++;
      if(id==5) lheBPartons++;
    }
    if(hepeup.ISTUP[i]==2 && abs(hepeup.IDUP[i])==23) lhePtZ=P4.pt();
  }
  // Check ranges
  if(lhePartons>npMax) lhePartons=npMax;
  if(lhePtZ>ptMax) lhePtZ=ptMax;
  if(lheHT>htMax) lheHT=htMax;
  // Get numbers
  int bin=Sum->FindBin(lhePartons, lhePtZ, lheHT);
  float num=Num->GetBinContent(bin);
  float den=Sum->GetBinContent(bin);
  if(den!=den || num!=num) return 0.;
  if(den<=0. || num<=0.) return 0.;
  return num/den;
}
*/
