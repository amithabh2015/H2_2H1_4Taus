#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauGetPFParticles.h"

HTauTauGetPFParticles::HTauTauGetPFParticles(edm::InputTag tag) 
{
  _pfPartColName = tag;
}
HTauTauGetPFParticles::~HTauTauGetPFParticles() {}
void HTauTauGetPFParticles::getPFParticles(const edm::Event& iEvent, patPFParticleVector & output) {
  
  using namespace edm;
  using namespace pat;

  output.clear();

  Handle<PFParticleCollection> pfPartHandle;
  iEvent.getByLabel(_pfPartColName,pfPartHandle);
  
  for (unsigned int j=0;j<pfPartHandle->size();j++) {
    PFParticleRef thisPFPartRef(pfPartHandle,j);
    output.push_back(thisPFPartRef);
  }

}
