#include "bschwrel.h"

#include <iostream>

namespace libhwrel {

std::shared_ptr<Response> BSC_HWReliabilityMonitor::perform_analysis(std::shared_ptr<Request> req) {
    std::shared_ptr<Response> output;

    // INSERT YOUR CODE IN THIS METHOD
    // Fill the output variable

    // Example code:

    std::cout << "Request type = " << (int)req->get_resource_type() << std::endl;
    std::cout << "Technology type = " << (int)req->get_technology_type() << std::endl;

    auto &temperatures = req->get_temperatures();

    for (const auto &t : temperatures) {
        std::cout << "T = " << t.temperature << std::endl;
    }

    // You can also cast "req" to td::shared_ptr<RequestCPU> or td::shared_ptr<RequestMEM> or ...  

    output = std::make_shared<Response>(0.0001);
    return output;
}


} // namespace libhwrel
