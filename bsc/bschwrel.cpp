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
    corep state_skt;
    memory state_mem;
    // gpu state_g;
    // gpu_mem gm ;

    /*
    *
    * 
    *  resource_type_t::CPU
    * 
    * 
    */
    if( resource == resource_type_t::CPU && tech == technology_type_t::SILICON){
        

        state_skt.phys = nr_cores;
        state_skt.original_fitALU = initial_fit / ( state_skt.phys * state_skt.bl_size );
        state_skt.original_fitLD_ST_AGU = initial_fit / ( state_skt.phys * state_skt.bl_size );
        state_skt.original_fitROB = initial_fit / ( state_skt.phys * state_skt.bl_size );
        state_skt.original_fitL1I = initial_fit / ( state_skt.phys * state_skt.bl_size );
        state_skt.original_fitL1D = initial_fit / ( state_skt.phys * state_skt.bl_size );
        state_skt.original_fitL2 = initial_fit / ( state_skt.phys * state_skt.bl_size );
        state_skt.original_fitL3 = initial_fit / ( state_skt.phys * state_skt.bl_size );
        
        state_skt.last_fitALU = state_skt.original_fitALU;
        state_skt.last_fitLD_ST_AGU = state_skt.original_fitLD_ST_AGU;
        state_skt.last_fitROB = state_skt.original_fitROB;
        state_skt.last_fitL1I = state_skt.original_fitL1I;
        state_skt.last_fitL1D = state_skt.original_fitL1D;
        state_skt.last_fitL2 = state_skt.original_fitL2;
        state_skt.last_fitL3 = state_skt.original_fitL3;
        state_skt.second_past = 0;
        
        state_reliability.failure_probability = 0.0 ;
        state_reliability.state = std::make_shared<corep>(state_skt);
        state_reliability.state_size = sizeof(state_skt);


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

 
   }



    /*
    *
    * 
    *  resource_type_t::MEMORY_GPU
    * 
    * 
    */
    else if(resource == resource_type_t::MEMORY_GPU && tech == technology_type_t::SILICON){
        

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
        memory mem_actual ; 
       
        mem_actual.second_past = (*state_prev).second_past;
        mem_actual.last_fitMEM =  (* state_prev).last_fitMEM;
        mem_actual.original_fitMEM =  (* state_prev).original_fitMEM;
        
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
        
       


    }else if(req->get_resource_type() == resource_type_t::MEMORY_GPU && req->get_technology_type()==technology_type_t::SILICON){




        
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
