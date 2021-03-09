#ifndef BSC_LIBHWREL_H_
#define BSC_LIBHWREL_H_

#include "libhwrel.h"

namespace libhwrel {

class BSC_HWReliabilityMonitor : public HWReliabilityMonitor {

public:

    BSC_HWReliabilityMonitor();
    virtual ~BSC_HWReliabilityMonitor() = default;

    virtual std::shared_ptr<Response> perform_analysis(std::shared_ptr<Request> req) override;
    virtual reliability_state_t init(resource_type_t resource, technology_type_t tech, long double initial_fit, unsigned int nr_cores) override;
   
  
private: 
  

    struct corep {
        int phys;

        int bl_size = 7;
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

        long double last_fitALU;
        long double last_fitLD_ST_AGU;
        long double last_fitROB;
        long double last_fitL1I;
        long double last_fitL1D;
        long double last_fitL2;
        long double last_fitL3;

        unsigned long long second_past;
    };

    struct memory {
        /*
        * Should be aggregate FIT, whole mem skt
        */
        long double original_fitMEM;
        long double last_fitMEM;
        unsigned long long second_past;
    };

    struct gpu {
        
        int bl_size = 8;

        long double original_fit_CUDA_CORE;
        long double original_fit_slot_instr;
        long double original_fit_DP_CORE;
        long double original_fit_SFU_CORE;
        long double original_fit_LDST_CORE;
        long double original_fit_L1_unified;
        long double original_fit_shared_sm;
        long double original_fit_l2;
        
        long double last_fit_CUDA_CORE;
        long double last_fit_slot_instr;
        long double last_fit_DP_CORE;
        long double last_fit_SFU_CORE;
        long double last_fit_LDST_CORE;
        long double last_fit_L1_unified;
        long double last_fit_shared_sm;
        long double last_fit_l2;
        
        unsigned long long second_past;
    };

    struct gpu_mem {
        long double last_fitGPU_MEM;
        long double original_fitGPU_MEM;
        unsigned long long second_past;
    };



    long double getAccelerationFactor(float temperature);
    long double getFIT(long double fit_previous, long double fit_original, unsigned long long accumulated, long double avf, unsigned long long time, long double af);

    long double getCompressProbability(long double fit,unsigned long long second_passed);


};



} // namespace libhwrel

#endif
