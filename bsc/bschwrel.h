#ifndef BSC_LIBHWREL_H_
#define BSC_LIBHWREL_H_

#include "libhwrel.h"

namespace libhwrel {

class BSC_HWReliabilityMonitor : public HWReliabilityMonitor {

public:

    virtual ~BSC_HWReliabilityMonitor() = default;

    virtual std::shared_ptr<Response> perform_analysis(std::shared_ptr<Request> req) override;

};



} // namespace libhwrel

#endif
