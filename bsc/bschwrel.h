#ifndef BSC_LIBHWREL_H_
#define BSC_LIBHWREL_H_

#include "libhwrel.h"

namespace libhwrel {

class BSC_HWReliabilityMonitor : public HWReliabilityMonitor {

public:

    BSC_HWReliabilityMonitor();
    virtual ~BSC_HWReliabilityMonitor() = default;

    virtual std::shared_ptr<Response> perform_analysis(std::shared_ptr<Request> req) override;
   
  
private: 
    /*
        seconds passed in general
    */

     unsigned long long secondsCPU;
     unsigned long long secondsMEM;
     unsigned long long secondsGPU;
     unsigned long long secondsAccelerator;
     
    /*
        starting  fit value
    */

     long double fitOriginalCPU;
     long double fitOriginalMEM;
     long double fitOriginalGPU;
     long double fitOriginalAccelerator;

    /*
        size blocks
    */
    int bl_size;
     /*
        blocks
     */
     long double original_fitALU;
     long double original_fitLD_ST_AGU;
     long double original_fitROB;
     long double original_fitL1I;
     long double original_fitL1D;
     long double original_fitL2;
     long double original_fitL3;

    /*
        last fit per block calculated 
        per previous requestCPU
        * last_fitALU[0] -> alu fit per core 0 
        * last_fitALU[1] -> alu fit per core 1 
        * .. 
    */

     std::vector<long double> last_fitALU;
     std::vector<long double> last_fitLD_ST_AGU;
     std::vector<long double> last_fitROB;
     std::vector<long double> last_fitL1I;
     std::vector<long double> last_fitL1D;
     std::vector<long double> last_fitL2;
     std::vector<long double> last_fitL3;

     long double lastProbCPU;
     long double lastProbMEM;
     long double lastProbGPU;
     long double lastProbAccelerator;

     long double getAccelerationFactor(float temperature);
     long double getFIT(long double fit_previous, long double fit_original, unsigned long long accumulated, long double avf, unsigned long long time, long double af);

     long double getFitStandardCPU();
     long double getFitStandardMEM();
     long double getFitStandardGPU();
     long double getFitStandardAccelerator();

     long double getFit_processing(std::vector<core> cv ,  unsigned long long  second_chunk, float temperature_avg, long double acceleration_factor_ld);


     long double getCompressProbability(long double fit,unsigned long long second_passed);
     long double getUncompressFit(long double probability,resource_type_t res_type, technology_type_t tech_type);

    void checkPast(std::shared_ptr<Request> req);
    void checkTimeRequest(std::shared_ptr<Request> req);
    void checkTemperatureSize(std::shared_ptr<Request> req);
    /*
        get single counter value
    */
    unsigned long long getALU(std::vector<PerfCounter> v);
    unsigned long long getLD(std::vector<PerfCounter> v);
    unsigned long long getST(std::vector<PerfCounter> v);
    unsigned long long getAGU(std::vector<PerfCounter> v);
    unsigned long long getL1IH(std::vector<PerfCounter> v);
    unsigned long long getL1IM(std::vector<PerfCounter> v);
    unsigned long long getL2H(std::vector<PerfCounter> v);
    unsigned long long getL2M(std::vector<PerfCounter> v);
    unsigned long long getL3H(std::vector<PerfCounter> v);
    unsigned long long getL3M(std::vector<PerfCounter> v);
    unsigned long long getUopsIssued(std::vector<PerfCounter> v);
    unsigned long long getUopsRetired(std::vector<PerfCounter> v);

    unsigned long long getCyclesMax(std::vector<core> v);

    void init_FIT_original_block(int cores_physical_number, int socket_number);
    void init_FIT_last_block(int size);
    void incrementPastSecond(unsigned long long value,resource_type_t res_type, technology_type_t tech_type);
    
    void reset();
    unsigned long long getPastSecondCPU();
    unsigned long long getPastSecondMEM();
    unsigned long long getPastSecondGPU();
    unsigned long long getPastSecondAccelerator();

    long double getFIT_auxiliar();

    int getSKT_number(std::vector<core> cv);     
    int getPhyCoreSize(int skt, std::vector<core> &cv);
    int getSMT(int skt, std::vector<core> &cv);
    std::vector<int> getSocketVector(std::vector<core> v);
    void removeSMT(std::vector<core> &cv);
	
	
	/*debug*/
    void stamp(std::vector<core> c);
    
    



};



} // namespace libhwrel

#endif
