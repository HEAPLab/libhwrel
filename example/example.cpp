#include "bsc/bschwrel.cpp"

#include <iostream>

using namespace libhwrel;

int main() {

    HWReliabilityMonitor *hwrel = new BSC_HWReliabilityMonitor();

    reliability_state_t state = hwrel->init(resource_type_t::CPU, technology_type_t::SILICON, 12345, 2);

    std::shared_ptr<RequestCPUCore> cpu = std::make_shared<RequestCPUCore>(technology_type_t::SILICON, 1200);
    temperature_t t1 = {59.51, std::chrono::system_clock::now()-std::chrono::seconds(60)};
    temperature_t t2 = {60.12, std::chrono::system_clock::now()};
    cpu->push_temperature(t1);
    cpu->push_temperature(t2);
    cpu->set_state(state);

    auto response = hwrel->perform_analysis(cpu);

    std::cout << "Fault probability: " << response->get_fail_probability() << std::endl;

    delete hwrel;
    return 0;
    
}
