/**
 * @class BuildPatCandTab
 *
 * @brief This object builds the relational tables associating the Monte Carlo with the
 *        TkrRecon pattern recognition track candidates. 
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/McToolBox/src/BuildPatCandTab.h,v 1.1.1.1 2004/02/19 22:58:18 usher Exp $
 */
#include "GaudiKernel/IDataProviderSvc.h"

#ifndef BuildPatCandTab_h
#define BuildPatCandTab_h

namespace Event {

class BuildPatCandTab 
{
public:
    /// Standard Gaudi Tool interface constructor
    BuildPatCandTab(IDataProviderSvc* dataSvc);
   ~BuildPatCandTab();

private:
};

};

#endif
