/**
 * @class BuildMcTracks
 *
 * @brief This object is responsible for creating the various relational tables used in the
 *        Monte Carlo analysis of an event. Currently this includes the following:
 *        1) A table relating McParticles with TkrClusters
 *        2) A table relating McParticles with existing TkrCluster<->McPositionHit relations
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/McToolBox/src/BuildMcTracks.h,v 1.1.1.1 2004/02/19 22:58:18 usher Exp $
 */
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartRefVector.h"
#include "Event/MonteCarlo/McParticle.h"

#ifndef BuildMcTracks_h
#define BuildMcTracks_h

namespace Event {

class BuildMcTracks 
{
public:
    /// Standard Gaudi Tool interface constructor
    BuildMcTracks(IDataProviderSvc* dataSvc);
   ~BuildMcTracks();

private:
};

};

#endif
