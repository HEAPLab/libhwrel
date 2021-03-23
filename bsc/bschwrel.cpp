#include "bschwrel.h"

#include <iostream>
#include <cmath>
#include <iomanip>

namespace libhwrel
{

BSC_HWReliabilityMonitor::BSC_HWReliabilityMonitor()
{


}



reliability_state_t BSC_HWReliabilityMonitor::init(resource_type_t resource, technology_type_t tech, long double initial_fit, unsigned int nr_cores)
{

    reliability_state_t state_reliability ;
    corep state_corep;
    memory state_mem;
    gpu state_gpu;
    gpu_mem state_memgpu ;

    /*
    *
    * 
    *  resource_type_t::CPU
    * 
    * 
    */
    if( resource == resource_type_t::CPU && tech == technology_type_t::SILICON){
        

        state_corep.phys = nr_cores;
        state_corep.original_fitALU = initial_fit / ( state_corep.phys * state_corep.bl_size );
        state_corep.original_fitLD_ST_AGU = initial_fit / ( state_corep.phys * state_corep.bl_size );
        state_corep.original_fitROB = initial_fit / ( state_corep.phys * state_corep.bl_size );
        state_corep.original_fitL1I = initial_fit / ( state_corep.phys * state_corep.bl_size );
        state_corep.original_fitL1D = initial_fit / ( state_corep.phys * state_corep.bl_size );
        state_corep.original_fitL2 = initial_fit / ( state_corep.phys * state_corep.bl_size );
        state_corep.original_fitL3 = initial_fit / ( state_corep.phys * state_corep.bl_size );
        
        state_corep.last_fitALU = state_corep.original_fitALU;
        state_corep.last_fitLD_ST_AGU = state_corep.original_fitLD_ST_AGU;
        state_corep.last_fitROB = state_corep.original_fitROB;
        state_corep.last_fitL1I = state_corep.original_fitL1I;
        state_corep.last_fitL1D = state_corep.original_fitL1D;
        state_corep.last_fitL2 = state_corep.original_fitL2;
        state_corep.last_fitL3 = state_corep.original_fitL3;

        state_corep.second_past = 0;
        
        state_reliability.failure_probability = 0.0 ;
        state_reliability.state = std::make_shared<corep>(state_corep);
        state_reliability.state_size = sizeof(state_corep);


    }
    /*
    *
    * 
    *  resource_type_t::MEMORY
    * 
    * 
    */
    else if (resource == resource_type_t::MEMORY && tech == technology_type_t::SILICON) {

        state_mem.original_fitMEM = initial_fit;
        state_mem.last_fitMEM = initial_fit;

        state_mem.second_past = 0;

        state_reliability.failure_probability = 0.0 ;
        state_reliability.state = std::make_shared<memory>(state_mem);
        state_reliability.state_size = sizeof(state_mem);
    }
    /*
    *
    * 
    *  resource_type_t::GPU
    * 
    * 
    */
   else if(resource == resource_type_t::GPU && tech == technology_type_t::SILICON){

        state_gpu.original_fit_CUDA_CORE = initial_fit / state_gpu.bl_size ;
        state_gpu.original_fit_slot_instr= initial_fit / state_gpu.bl_size ;
        state_gpu.original_fit_DP_CORE = initial_fit / state_gpu.bl_size ; 
        state_gpu.original_fit_SFU_CORE= initial_fit / state_gpu.bl_size ;
        state_gpu.original_fit_LDST_CORE= initial_fit / state_gpu.bl_size ; 
        state_gpu.original_fit_L1_unified= initial_fit / state_gpu.bl_size ; 
        state_gpu.original_fit_shared_sm = initial_fit / state_gpu.bl_size ; 
        state_gpu.original_fit_l2 = initial_fit / state_gpu.bl_size ;
        
        state_gpu.last_fit_CUDA_CORE = state_gpu.original_fit_CUDA_CORE ; 
        state_gpu.last_fit_slot_instr = state_gpu.original_fit_slot_instr;
        state_gpu.last_fit_DP_CORE =  state_gpu.original_fit_DP_CORE ; 
        state_gpu.last_fit_SFU_CORE = state_gpu.original_fit_SFU_CORE;
        state_gpu.last_fit_LDST_CORE = state_gpu.original_fit_LDST_CORE; 
        state_gpu.last_fit_L1_unified = state_gpu.original_fit_L1_unified; 
        state_gpu.last_fit_shared_sm = state_gpu.original_fit_shared_sm ; 
        state_gpu.last_fit_l2 = state_gpu.original_fit_l2 ; 
        
        state_gpu.second_past = 0 ;

        state_reliability.failure_probability = 0.0 ;
        state_reliability.state = std::make_shared<gpu>(state_gpu);
        state_reliability.state_size = sizeof(state_gpu);
 
   }



    /*
    *
    * 
    *  resource_type_t::MEMORY_GPU
    * 
    * 
    */
    else if(resource == resource_type_t::MEMORY_GPU && tech == technology_type_t::SILICON){
        
        state_memgpu.original_fitGPU_MEM = initial_fit;
        state_memgpu.last_fitGPU_MEM = initial_fit;

        state_memgpu.second_past = 0;
        

        state_reliability.failure_probability = 0.0 ;
        state_reliability.state = std::make_shared<gpu_mem>(state_memgpu);
        state_reliability.state_size = sizeof(state_memgpu);

    }
    
    return state_reliability ;
}

std::shared_ptr<Response> BSC_HWReliabilityMonitor::perform_analysis(std::shared_ptr<Request> req)
{
    
    float fail_probability;

    reliability_state_t actual_state;
    
    unsigned long long  second_chunk;
    unsigned long long seconds_past_loc;
    float temperature_avg;
    long double acceleration_factor_ld;
    long double fail_probability_previous;
    long double fit_calculated;
    long double fit_previous;
    long double prob;


    if (req->get_resource_type() == resource_type_t::CPU && req->get_technology_type()==technology_type_t::SILICON)
    {   
        std::shared_ptr<RequestCPUCore> core_request = std::dynamic_pointer_cast<RequestCPUCore>(req);
        auto &temperatures_vector = core_request->get_temperatures();
        fail_probability_previous = core_request->get_state().failure_probability;

        corep * state_prev = (corep * ) core_request->get_state().state.get();
        corep corep_actual ; 
        
        corep_actual.second_past = (*state_prev).second_past;
        corep_actual.last_fitALU = (*state_prev).last_fitALU;
        corep_actual.last_fitLD_ST_AGU = (*state_prev).last_fitLD_ST_AGU;
        corep_actual.last_fitROB = (*state_prev).last_fitROB;
        corep_actual.last_fitL1I = (*state_prev).last_fitL1I;
        corep_actual.last_fitL1D = (*state_prev).last_fitL1D;
        corep_actual.last_fitL2 = (*state_prev).last_fitL2;
        corep_actual.last_fitL3 = (*state_prev).last_fitL3;
        corep_actual.original_fitALU = ((*state_prev).original_fitALU);
        corep_actual.original_fitLD_ST_AGU = ((*state_prev).original_fitLD_ST_AGU);
        corep_actual.original_fitROB = ((*state_prev).original_fitROB);
        corep_actual.original_fitL1I = ((*state_prev).original_fitL1I);
        corep_actual.original_fitL1D = ((*state_prev).original_fitL1D);
        corep_actual.original_fitL2 = ((*state_prev).original_fitL2);
        corep_actual.original_fitL3 = ((*state_prev).original_fitL3);

        seconds_past_loc  =  (*state_prev).second_past;
        
        for (int t = 0; t < temperatures_vector.size() - 1; t++)
        {
            /*Starting data*/
            second_chunk = (std::chrono::duration<double>(temperatures_vector[t + 1].epoch.time_since_epoch())).count() -
                           (std::chrono::duration<double>(temperatures_vector[t].epoch.time_since_epoch())).count();
            temperature_avg = (temperatures_vector[t].temperature + temperatures_vector[t + 1].temperature) / 2;
            acceleration_factor_ld = getAccelerationFactor(temperature_avg);

            long double fit_core;
            unsigned long long  cycles = core_request->get_PC(perf_counter_type_t::CYCLES).get_value();
            long double fit_processing = 0 ;
            
            /*
            *   Activity blocks
            */
            long double act_alu=0;
            long double act_ld_st_agu=0;
            long double ratio_issued=0;
            long double ratio_retired=0;
            long double act_rob=0;
            long double act_L1I=0;
            long double act_L1D=0;
            long double act_L2=0;
            long double act_L3=0;

            fit_calculated=0;


            act_alu = (long double)(core_request->get_PC(perf_counter_type_t::PORT_0).get_value() +
                                    core_request->get_PC(perf_counter_type_t::PORT_1).get_value() +
                                    core_request->get_PC(perf_counter_type_t::PORT_5).get_value() +
                                    core_request->get_PC(perf_counter_type_t::PORT_6).get_value()
                                    )/ (cycles * 4);

            act_ld_st_agu = (long double)(
                                            /*
                                            * LD
                                            * */
                                            core_request->get_PC(perf_counter_type_t::PORT_2).get_value() + 
                                            core_request->get_PC(perf_counter_type_t::PORT_3).get_value() + 
                                            /*
                                            * ST
                                            */
                                            core_request->get_PC(perf_counter_type_t::PORT_4).get_value() +
                                            /*
                                            *  AGU
                                            */
                                            core_request->get_PC(perf_counter_type_t::PORT_7).get_value() 
                                        ) / (cycles * 4);

            /*
            * ROB
            */
            ratio_issued = (long double)core_request->get_PC(perf_counter_type_t::UOPS_ISSUED_ANY).get_value()  / (4 * cycles);
            ratio_retired = (long double)core_request->get_PC(perf_counter_type_t::UOPS_RETIRED).get_value() / (4 * cycles);
            act_rob = ratio_retired / ratio_issued;

            /*
            *   L1 
            */

            act_L1I = (long double)core_request->get_PC(perf_counter_type_t::L1I_HIT).get_value() / cycles;
            act_L1I +=(long double)core_request->get_PC(perf_counter_type_t::L1I_MISS).get_value() / cycles;

            act_L1D = (long double)core_request->get_PC(perf_counter_type_t::PORT_2).get_value() / cycles;
            act_L1D += (long double)core_request->get_PC(perf_counter_type_t::PORT_3).get_value() / cycles;
           
            
            /*
            * L2
            */

            act_L2 = (long double)( /*Hit*/
                                     (core_request->get_PC(perf_counter_type_t::L2_RQST_REFERENCES).get_value() - 
                                        core_request->get_PC(perf_counter_type_t::L2_RQST_MISS).get_value()) 
                                    ) / (cycles * 3);
            act_L2 += (long double)core_request->get_PC(perf_counter_type_t::L2_RQST_MISS).get_value() / cycles;
            /*
            * L3
            */
            act_L3 = (long double)( /*Hit*/
                                    (core_request->get_PC(perf_counter_type_t::L3_REFERENCES).get_value() - 
                                        core_request->get_PC(perf_counter_type_t::L3_MISS).get_value()) 
                                    ) / cycles;
            act_L3 += (long double)core_request->get_PC(perf_counter_type_t::L3_MISS).get_value() / cycles;

            /*
            * functionality unity core
            */
            corep_actual.last_fitALU=getFIT( corep_actual.last_fitALU,corep_actual.original_fitALU, seconds_past_loc, act_alu , second_chunk, acceleration_factor_ld);
            corep_actual.last_fitLD_ST_AGU=getFIT(corep_actual.last_fitLD_ST_AGU,corep_actual.original_fitLD_ST_AGU, seconds_past_loc, act_ld_st_agu, second_chunk, acceleration_factor_ld);
            /*
            Rob
            */
            corep_actual.last_fitROB = getFIT(corep_actual.last_fitROB,corep_actual.original_fitROB,seconds_past_loc, act_rob, second_chunk, acceleration_factor_ld);
            /*
            Caches
            */
            corep_actual.last_fitL1I = getFIT(corep_actual.last_fitL1I,corep_actual.original_fitL1I ,seconds_past_loc,act_L1I, second_chunk, acceleration_factor_ld);
            corep_actual.last_fitL1D = getFIT(corep_actual.last_fitL1I,corep_actual.original_fitL1D ,seconds_past_loc,act_L1D, second_chunk, acceleration_factor_ld);
            corep_actual.last_fitL2 = getFIT(corep_actual.last_fitL2,corep_actual.original_fitL2, seconds_past_loc,act_L2, second_chunk, acceleration_factor_ld);
            corep_actual.last_fitL3 = getFIT(corep_actual.last_fitL3,corep_actual.original_fitL3, seconds_past_loc,act_L3, second_chunk, acceleration_factor_ld);

            fit_calculated+=corep_actual.last_fitALU;
            fit_calculated+=corep_actual.last_fitLD_ST_AGU;
            fit_calculated+=corep_actual.last_fitROB;
            fit_calculated+=corep_actual.last_fitL1I;
            fit_calculated+=corep_actual.last_fitL1D;
            fit_calculated+=corep_actual.last_fitL2;
            fit_calculated+=corep_actual.last_fitL3;   
        
            prob = getCompressProbability(fit_calculated,  seconds_past_loc + second_chunk);
            fail_probability = (prob * 1000);
            seconds_past_loc = seconds_past_loc + second_chunk;
            corep_actual.second_past = seconds_past_loc;

            
        }
    
        /*
        *   Saving to generate response
        */
       actual_state.failure_probability = fail_probability;
       actual_state.state = std::make_shared<corep>(corep_actual);
       actual_state.state_size = sizeof(corep);

    }
    
    else if(req->get_resource_type() == resource_type_t::MEMORY && req->get_technology_type() == technology_type_t::SILICON){
        
        std::shared_ptr<RequestMEM> mem_request = std::dynamic_pointer_cast<RequestMEM>(req);
        auto &temperatures_vector = mem_request->get_temperatures(); 
        fail_probability_previous = mem_request->get_state().failure_probability;

        memory * state_prev = (memory * ) mem_request->get_state().state.get();
        memory mem_actual = (*state_prev) ; 
        
        seconds_past_loc  =  (*state_prev).second_past;

        for (int t = 0; t < temperatures_vector.size() - 1; t++)
        {
            /*Starting data*/
            second_chunk = (std::chrono::duration<double>(temperatures_vector[t + 1].epoch.time_since_epoch())).count() -
                           (std::chrono::duration<double>(temperatures_vector[t].epoch.time_since_epoch())).count();
            temperature_avg = (temperatures_vector[t].temperature + temperatures_vector[t + 1].temperature) / 2;
            acceleration_factor_ld = getAccelerationFactor(temperature_avg);
            fit_calculated = 0; /*reset*/

        
            long double usage = 0;
            unsigned long long total_max = (mem_request->get_band_per_skt());
            unsigned long long bytes =  (mem_request->get_PC(perf_counter_type_t::CAS_READ).get_value() 
                                            + mem_request->get_PC(perf_counter_type_t::CAS_WRITE).get_value() ) 
                                            * 64 ; 

            
            usage  = (long double)(bytes / second_chunk) / total_max ;

            mem_actual.last_fitMEM = getFIT(   mem_actual.last_fitMEM,
                                                mem_actual.original_fitMEM, 
                                                seconds_past_loc, 
                                                usage,
                                                second_chunk, 
                                                acceleration_factor_ld);

            fit_calculated = mem_actual.last_fitMEM ;
             
            prob = getCompressProbability(fit_calculated, seconds_past_loc + second_chunk);
            fail_probability = (prob * 1000);
            seconds_past_loc = seconds_past_loc + second_chunk ;
            mem_actual.second_past = seconds_past_loc;
     
        }
        
        /*
        *   Saving to generate response
        */
       actual_state.failure_probability = fail_probability;
       actual_state.state = std::make_shared<memory>(mem_actual);
       actual_state.state_size = sizeof(memory);

        
    }else if(req->get_resource_type() == resource_type_t::GPU && req->get_technology_type()==technology_type_t::SILICON){
        
        std::shared_ptr<RequestGPU> gpu_request = std::dynamic_pointer_cast<RequestGPU>(req);
        auto &temperatures_vector = gpu_request->get_temperatures();
        fail_probability_previous = gpu_request->get_state().failure_probability;
        bool runtime = gpu_request->get_runtime_mode();

        gpu * state_prev = (gpu * ) gpu_request->get_state().state.get();
        gpu gpu_actual = (*state_prev); 
        
        seconds_past_loc  =  (*state_prev).second_past;


        /*
        * Not runtime values 
        */
        long double CC_SP_utilization = runtime ? 0 : ((long double) gpu_request->get_PC(perf_counter_type_t::SINGLE_PRECISION_FU_UTILIZATION).get_value())/10 ;
        long double CC_FP_utlization =  runtime ? 0 : ((long double) gpu_request->get_PC(perf_counter_type_t::HALF_PRECISION_FU_UTILIZATION).get_value())/10;
        long double DP_utilization = runtime ? 0 : ((long double) gpu_request->get_PC(perf_counter_type_t::DOUBLE_PRECISION_FU_UTILIZATION).get_value())/10;
        long double LDST_utilization = runtime ? 0 : ((long double) gpu_request->get_PC(perf_counter_type_t::LDST_FU_UTILIZATION).get_value())/10 ;
        long double SFU_utilization = runtime ? 0 : ((long double) gpu_request->get_PC(perf_counter_type_t::SPECIAL_FU_UTILIZATION).get_value())/10 ;
        long double SLOT_utilization = runtime ? 0 : ((long double) gpu_request->get_PC(perf_counter_type_t::ISSUE_SLOT_UTILIZATION).get_value())/10 ;
        long double L1_uni_utilization = runtime ? 0 : ((long double) gpu_request->get_PC(perf_counter_type_t::TEX_UTILIZATION).get_value())/10 ;
        long double SH_utilization = runtime ? 0 : ((long double) gpu_request->get_PC(perf_counter_type_t::SHARED_UTILIZATION).get_value())/10 ;
        long double l2_utilization = runtime ? 0 : ((long double) gpu_request->get_PC(perf_counter_type_t::L2_UTILIZATION).get_value())/10 ;

        long double CC_acf = CC_SP_utilization + CC_FP_utlization ;
        if( CC_acf > 1.0 ) CC_acf = 1.0;


        /*
        *   Runtime value
        */
        long double gpu_utilization = runtime ? ((long double) gpu_request->get_PC(perf_counter_type_t::GPU_UTILIZATION).get_value())/100 : 0;

        for (int t = 0; t < temperatures_vector.size() - 1; t++)
        {
            /*Starting data*/
            second_chunk = (std::chrono::duration<double>(temperatures_vector[t + 1].epoch.time_since_epoch())).count() -
                           (std::chrono::duration<double>(temperatures_vector[t].epoch.time_since_epoch())).count();
            temperature_avg = (temperatures_vector[t].temperature + temperatures_vector[t + 1].temperature) / 2;
            acceleration_factor_ld = getAccelerationFactor(temperature_avg);
         
            fit_calculated = 0;
            
            if(runtime){

                /*
                * global general previous GPU fit
                */
            
                fit_calculated += gpu_actual.last_fit_CUDA_CORE;
                fit_calculated += gpu_actual.last_fit_DP_CORE;
                fit_calculated += gpu_actual.last_fit_SFU_CORE;
                fit_calculated += gpu_actual.last_fit_LDST_CORE;
                fit_calculated += gpu_actual.last_fit_slot_instr;
                fit_calculated += gpu_actual.last_fit_L1_unified;
                fit_calculated += gpu_actual.last_fit_shared_sm;
                fit_calculated += gpu_actual.last_fit_l2;

                /*
                * proportions fit for every block to general actual fit
                */

                long double  CUDA_CORE_fit_percentage =  gpu_actual.last_fit_CUDA_CORE / fit_calculated;
                long double  DP_CORE_fit_percentage =  gpu_actual.last_fit_DP_CORE / fit_calculated;
                long double  SFU_CORE_fit_percentage =  gpu_actual.last_fit_SFU_CORE / fit_calculated;
                long double  LDST_CORE_fit_percentage =  gpu_actual.last_fit_LDST_CORE / fit_calculated;
                long double  slot_instr_fit_percentage =  gpu_actual.last_fit_slot_instr / fit_calculated;
                long double  L1_unified_fit_percentage =  gpu_actual.last_fit_L1_unified / fit_calculated;
                long double  shared_fit_percentage =  gpu_actual.last_fit_shared_sm / fit_calculated;
                long double  L2_fit_percentage =  gpu_actual.last_fit_l2 / fit_calculated;

                long double fit_original_gpu = 0;
                
                /*
                * getting original gpu fit above initial blocks
                */

                fit_original_gpu += gpu_actual.original_fit_CUDA_CORE;
                fit_original_gpu += gpu_actual.original_fit_DP_CORE;
                fit_original_gpu += gpu_actual.original_fit_SFU_CORE;
                fit_original_gpu += gpu_actual.original_fit_LDST_CORE;
                fit_original_gpu += gpu_actual.original_fit_slot_instr;
                fit_original_gpu += gpu_actual.original_fit_L1_unified;
                fit_original_gpu += gpu_actual.original_fit_shared_sm;
                fit_original_gpu += gpu_actual.original_fit_l2;

                fit_calculated = getFIT(fit_calculated ,fit_original_gpu ,seconds_past_loc, gpu_utilization, second_chunk, acceleration_factor_ld); 

                /*
                * saving to maintain same proportions
                */

                gpu_actual.last_fit_CUDA_CORE = fit_calculated * CUDA_CORE_fit_percentage;
                gpu_actual.last_fit_DP_CORE = fit_calculated * DP_CORE_fit_percentage;
                gpu_actual.last_fit_SFU_CORE = fit_calculated * SFU_CORE_fit_percentage;
                gpu_actual.last_fit_LDST_CORE = fit_calculated * LDST_CORE_fit_percentage;
                gpu_actual.last_fit_slot_instr = fit_calculated * slot_instr_fit_percentage;
                gpu_actual.last_fit_L1_unified = fit_calculated * L1_unified_fit_percentage;
                gpu_actual.last_fit_shared_sm = fit_calculated * CUDA_CORE_fit_percentage;
                gpu_actual.last_fit_l2 = fit_calculated * L2_fit_percentage;


            }else{
                gpu_actual.last_fit_CUDA_CORE = getFIT( gpu_actual.last_fit_CUDA_CORE,gpu_actual.original_fit_CUDA_CORE,seconds_past_loc,CC_acf, second_chunk, acceleration_factor_ld);
                gpu_actual.last_fit_DP_CORE = getFIT( gpu_actual.last_fit_DP_CORE,gpu_actual.original_fit_DP_CORE,seconds_past_loc,DP_utilization , second_chunk, acceleration_factor_ld);
                gpu_actual.last_fit_SFU_CORE = getFIT( gpu_actual.last_fit_SFU_CORE,gpu_actual.original_fit_SFU_CORE,seconds_past_loc,SFU_utilization , second_chunk, acceleration_factor_ld);
                gpu_actual.last_fit_LDST_CORE = getFIT( gpu_actual.last_fit_LDST_CORE,gpu_actual.original_fit_LDST_CORE,seconds_past_loc,LDST_utilization , second_chunk, acceleration_factor_ld);
                gpu_actual.last_fit_slot_instr = getFIT( gpu_actual.last_fit_slot_instr,gpu_actual.original_fit_slot_instr,seconds_past_loc,SLOT_utilization , second_chunk, acceleration_factor_ld);
                gpu_actual.last_fit_L1_unified = getFIT( gpu_actual.last_fit_L1_unified,gpu_actual.original_fit_L1_unified,seconds_past_loc,L1_uni_utilization , second_chunk, acceleration_factor_ld);
;                gpu_actual.last_fit_l2 = getFIT( gpu_actual.last_fit_l2,gpu_actual.original_fit_l2,seconds_past_loc,l2_utilization , second_chunk, acceleration_factor_ld);

                fit_calculated += gpu_actual.last_fit_CUDA_CORE;
                fit_calculated += gpu_actual.last_fit_DP_CORE;
                fit_calculated += gpu_actual.last_fit_SFU_CORE;
                fit_calculated += gpu_actual.last_fit_LDST_CORE;
                fit_calculated += gpu_actual.last_fit_slot_instr;
                fit_calculated += gpu_actual.last_fit_L1_unified;
                fit_calculated += gpu_actual.last_fit_shared_sm;
                fit_calculated += gpu_actual.last_fit_l2;
            }

            

            prob = getCompressProbability(fit_calculated, seconds_past_loc + second_chunk);
            fail_probability = (prob * 1000);
            seconds_past_loc = seconds_past_loc + second_chunk ;
            fit_previous = fit_calculated;
            gpu_actual.second_past = seconds_past_loc;
        }

        /*
        *   Saving to generate response
        */
        actual_state.failure_probability = fail_probability;
        actual_state.state = std::make_shared<gpu>(gpu_actual);
        actual_state.state_size = sizeof(gpu);



    }else if(req->get_resource_type() == resource_type_t::MEMORY_GPU && req->get_technology_type()==technology_type_t::SILICON){

            std::shared_ptr<RequestMEM_GPU> mem_gpu_request = std::dynamic_pointer_cast<RequestMEM_GPU>(req);
            auto &temperatures_vector = mem_gpu_request->get_temperatures();
            fail_probability_previous = mem_gpu_request->get_state().failure_probability;
            bool runtime = mem_gpu_request->get_runtime_mode();


            gpu_mem * state_prev = (gpu_mem * ) mem_gpu_request->get_state().state.get();
            gpu_mem gpumem_actual = (*state_prev); 
            
            seconds_past_loc  =  (*state_prev).second_past;
            fit_previous = (*state_prev).last_fitGPU_MEM;

            /**
             * 
             *  from 0 to 10 
             */
            
            long double memory_gpu = runtime ?
                                     ((long double)(int)mem_gpu_request->get_PC(perf_counter_type_t::MEMORY_UTILIZATION).get_value()) /100 :
                                     ((long double)(int)mem_gpu_request->get_PC(perf_counter_type_t::DRAM_UTILIZATION).get_value()) /10;

            for (int t = 0; t < temperatures_vector.size() - 1; t++)
            {
                /*Starting data*/
                second_chunk = (std::chrono::duration<double>(temperatures_vector[t + 1].epoch.time_since_epoch())).count() -
                            (std::chrono::duration<double>(temperatures_vector[t].epoch.time_since_epoch())).count();
                temperature_avg = (temperatures_vector[t].temperature + temperatures_vector[t + 1].temperature) / 2;
                acceleration_factor_ld = getAccelerationFactor(temperature_avg);

                fit_calculated = getFIT(fit_previous, (*state_prev).original_fitGPU_MEM ,seconds_past_loc, memory_gpu, second_chunk, acceleration_factor_ld); 

                prob = getCompressProbability(fit_calculated, seconds_past_loc + second_chunk);
                fail_probability = (prob * 1000);
                seconds_past_loc = seconds_past_loc + second_chunk;
                fit_previous = fit_calculated;
            }
            
            
            gpumem_actual.last_fitGPU_MEM = fit_previous;
            gpumem_actual.original_fitGPU_MEM =  (*state_prev).original_fitGPU_MEM;
            gpumem_actual.second_past =  seconds_past_loc;

            actual_state.failure_probability = fail_probability;
            actual_state.state = std::make_shared<gpu_mem>(gpumem_actual);
            actual_state.state_size = sizeof(gpu_mem);


        
    }else if(req->get_resource_type() == resource_type_t::MEMORY && req->get_technology_type()==technology_type_t::FPGA){

    //     std::shared_ptr<RequestAccelerator> accelerator_request = std::dynamic_pointer_cast<RequestAccelerator>(req);
    //     auto &temperatures_vector = accelerator_request->get_temperatures();
    //     fail_probability_previous = accelerator_request->get_previous_state().failure_probability;

    //     if(fail_probability_previous == 0.0) fit_previous = getFitStandardAccelerator();    
    //     else {

    //         fit_previous = getUncompressFit((long double) this->lastProbAccelerator/1000,
	// 										accelerator_request->get_resource_type(),
	// 										accelerator_request->get_technology_type());
    //     }
    //     for (int t = 0; t < temperatures_vector.size() - 1; t++)
    //     {
    //         /*Starting data*/
    //         second_chunk = (std::chrono::duration<double>(temperatures_vector[t + 1].epoch.time_since_epoch())).count() -
    //                        (std::chrono::duration<double>(temperatures_vector[t].epoch.time_since_epoch())).count();
    //         temperature_avg = (temperatures_vector[t].temperature + temperatures_vector[t + 1].temperature) / 2;
    //         acceleration_factor_ld = getAccelerationFactor(temperature_avg);

    //         fit_calculated = getFIT(fit_previous,this->fitOriginalAccelerator, getPastSecondAccelerator(),(double)accelerator_request->get_occupancy()/1000, second_chunk, acceleration_factor_ld);
    //         prob = getCompressProbability(fit_calculated, getPastSecondGPU() + second_chunk);
    //         this->lastProbAccelerator = prob * 1000;
    //         fail_probability = (float)(prob * 1000);
    //         incrementPastSecond(second_chunk,accelerator_request->get_resource_type(),accelerator_request->get_technology_type());
    //         fit_previous = fit_calculated;
    //     }  

    }


    return std::make_shared<Response>(actual_state);
}

long double BSC_HWReliabilityMonitor::getAccelerationFactor(float temperature)
{

    long double act_energy = 0.7;
    long double boltzmann = 8.617333262145 * (1 / powl(10, 5));
    long double t_use = 328.0;
    long double t_stress = (long double)temperature + 273.15;

    long double acceleration_factor = (powl(M_El, ((act_energy / boltzmann) * (1.0 / t_use - 1.0 / t_stress))));

    return acceleration_factor;
}
long double BSC_HWReliabilityMonitor::getFIT(long double fit_previous,long double fit_original ,unsigned long long accumulated, long double avf, unsigned long long time, long double af)
{
    /*
    *   avf 0 > 1 
    * 
    */
    
   
    long double f_rate_second = fit_original / (1000000000L * 3600L);

    unsigned long long t_standard = 1000000000L * 3600L;
    unsigned long long t_previous = (fit_previous / fit_original) * t_standard;

    unsigned long long t_add = time + (time * avf * af);
    unsigned long long t_left_std = t_standard - accumulated;

    unsigned long long t_calculated = (t_left_std + accumulated - time) + (t_previous - t_standard) + t_add;
    long double fit_new = f_rate_second * t_calculated;
    


    return fit_new;
}


long double BSC_HWReliabilityMonitor::getCompressProbability(long double fit, unsigned long long second_passed)
{
    long double f_rate_second = fit / (1000000000L * 3600L);
    long double prob = 1 - powl((long double)(M_E), -1 * (f_rate_second * second_passed));

    return prob;
}



}
