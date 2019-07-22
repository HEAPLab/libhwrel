#include "libhwrel.h"

#include <iostream>

namespace libhwrel {


class MyHWReliabilityMonitor : public HWReliabilityMonitor {

public:

	virtual ~MyHWReliabilityMonitor() = default;

	virtual std::shared_ptr<Response> perform_analysis(std::shared_ptr<Request> req);

};

std::shared_ptr<Response> MyHWReliabilityMonitor::perform_analysis(std::shared_ptr<Request> req) {

	// Do something with req ...

	std::shared_ptr<Response> my_response = std::make_shared<Response>(500);	// 50%

	return my_response;
}

} // libhwrel

int main() {

	using namespace libhwrel;

	MyHWReliabilityMonitor hw_rel_monitor;

	// Build the request
	std::shared_ptr<Request> request = std::make_shared<RequestCPU>(technology_type_t::SILICON, 1200, 4, 300);
	request->push_temperature({ 40.5, std::chrono::system_clock::now() });

	// Get the response
	auto response = hw_rel_monitor.perform_analysis(request);

	std::cout << "Probability of failure: " << (response->get_fail_probability()/1000.0) << std::endl;

	return 0;
}
