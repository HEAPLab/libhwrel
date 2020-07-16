#include "bschwrel.h"

#include <iostream>
#include <cmath>
#include <iomanip>

namespace libhwrel
{

BSC_HWReliabilityMonitor::BSC_HWReliabilityMonitor()
{

    /*
        Second pass from the starting 
        to the moment
    */
    this->secondsCPU = 0L;    
    this->secondsMEM = 0L;
    this->secondsGPU = 0L;
    this->secondsAccelerator = 0L;

    /*
        Standard FIT value for CPU and MEMORY, tuning
        as manufacturing says. Standard FIT is on
        1.000.000.000 hours. Standard CPU concept of all
        single CPU in one
    */                 
    this->fitOriginalCPU = 19750.510; 
    this->fitOriginalMEM = 36150.601;
    this->fitOriginalGPU = 25000.32;
    this->fitOriginalGPU = 350140.87;

    /*
        supposing the fit original divide in the 
        * block size
    */
    this->bl_size = 7;
    
    /*
        Last probability request calculated
    */

    this->lastProbCPU = 0L;
    this->lastProbMEM = 0L;
    this->lastProbGPU = 0L;
    this->lastProbAccelerator = 0L;
    
}

std::shared_ptr<Response> BSC_HWReliabilityMonitor::perform_analysis(std::shared_ptr<Request> req)
{
    /*
        Checking Request
        if temperatures vector size == 0 

    */

    if(req->get_temperatures().size() == 0){
        if(req->get_resource_type() == resource_type_t::CPU && req->get_technology_type()==technology_type_t::SILICON){
            return std::make_shared<Response>((float)this->lastProbCPU);
        }else if(req->get_resource_type() == resource_type_t::MEMORY && req->get_technology_type() == technology_type_t::SILICON){
            return std::make_shared<Response>((float)this->lastProbMEM);
        }else if(req->get_resource_type() == resource_type_t::GPU && req->get_technology_type() == technology_type_t::SILICON){
            return std::make_shared<Response>((float)this->lastProbGPU);
        }else if(req->get_resource_type() == resource_type_t::MEMORY && req->get_technology_type() == technology_type_t::FPGA){
            return std::make_shared<Response>((float)this->lastProbAccelerator);
        }
    }
    
    /*
     * check if some epoch of temperature vector is < 
     * of previous state epoch 
     * */
    checkPast(req); 
    
    /*
     * check if of temperatures vector epoch are 
     * in consecutio time
     * */
    
    checkTimeRequest(req);
    
    /*
     * check if temperatures vector size is <= 1
     * */
    
    checkTemperatureSize(req);



    std::shared_ptr<Response> output;
    float fail_probability;

    // INSERT YOUR CODE IN THIS METHOD
    // Fill the output variable

    // Example code:

    unsigned long long  second_chunk;
    float temperature_avg;
    long double acceleration_factor_ld;
    float fail_probability_previous;
    long double fit_calculated;
    long double fit_previous;
    long double prob;


    if (req->get_resource_type() == resource_type_t::CPU && req->get_technology_type()==technology_type_t::SILICON)
    {   

        std::shared_ptr<RequestCPU> cpu_request = std::dynamic_pointer_cast<RequestCPU>(req);
        auto &temperatures_vector = cpu_request->get_temperatures();
        fail_probability_previous = cpu_request->get_previous_state().failure_probability;

        std::vector<core> cv = cpu_request->getCoreVector();
        /*
         * if the core are hyperthreading, we remove
         * the logical core
         * transform to physical core vector
         * */
        removeSMT(cv); 
        int socket = getSKT_number(cv);

        double fit_processing=0;

		/*
		 * if the previous state is empty
		 * init some variable
		 * */
        if(fail_probability_previous == 0.0){    
            init_FIT_original_block(cv.size(),socket);
            init_FIT_last_block(cv.size());
        }  
        
        /*
         * perform calculation for cpus 
         * */
        for (int t = 0; t < temperatures_vector.size() - 1; t++)
        {
            /*Starting data*/
            second_chunk = (std::chrono::duration<double>(temperatures_vector[t + 1].epoch.time_since_epoch())).count() -
                           (std::chrono::duration<double>(temperatures_vector[t].epoch.time_since_epoch())).count();
            temperature_avg = (temperatures_vector[t].temperature + temperatures_vector[t + 1].temperature) / 2;
            acceleration_factor_ld = getAccelerationFactor(temperature_avg);
            /*
             * fit calculated is the sum of all
             * blocks in each core physical
             * */
            fit_calculated = getFit_processing(cv,second_chunk,temperature_avg, acceleration_factor_ld);
            prob = getCompressProbability(fit_calculated, getPastSecondCPU() + second_chunk);
            this->lastProbCPU = (prob * 1000);
            fail_probability = (float)(prob * 1000);
            incrementPastSecond(second_chunk, cpu_request->get_resource_type(),cpu_request->get_technology_type());

        }
    }
    
    else if(req->get_resource_type() == resource_type_t::MEMORY && req->get_technology_type() == technology_type_t::SILICON){
        
        std::shared_ptr<RequestMEM> mem_request = std::dynamic_pointer_cast<RequestMEM>(req);
        auto &temperatures_vector = mem_request->get_temperatures(); 
        fail_probability_previous = mem_request->get_previous_state().failure_probability;

       if(fail_probability_previous == 0.0){ 
            fit_previous = getFitStandardMEM();
        } else {
            /*
            After the first time.. it is 
            needed / 1000 , for 1000 format
            */
            fit_previous = getUncompressFit((long double)this->lastProbMEM /1000 , 
                                            mem_request->get_resource_type(),
                                            mem_request->get_technology_type());
        }  

        for (int t = 0; t < temperatures_vector.size() - 1; t++)
        {
            /*Starting data*/
            second_chunk = (std::chrono::duration<double>(temperatures_vector[t + 1].epoch.time_since_epoch())).count() -
                           (std::chrono::duration<double>(temperatures_vector[t].epoch.time_since_epoch())).count();
            temperature_avg = (temperatures_vector[t].temperature + temperatures_vector[t + 1].temperature) / 2;
            acceleration_factor_ld = getAccelerationFactor(temperature_avg);
            fit_calculated = getFIT(fit_previous,this->fitOriginalMEM, getPastSecondMEM(), (double)mem_request->get_occupancy()/1000 , 
                                    second_chunk, acceleration_factor_ld);
            prob = getCompressProbability(fit_calculated, getPastSecondMEM() + second_chunk);
            this->lastProbMEM = prob * 1000;
            fail_probability = (float)(prob * 1000);
            incrementPastSecond(second_chunk,mem_request->get_resource_type(),mem_request->get_technology_type());
            fit_previous = fit_calculated;
     
        }
        
    }else if(req->get_resource_type() == resource_type_t::GPU && req->get_technology_type()==technology_type_t::SILICON){
        
        std::shared_ptr<RequestGPU> gpu_request = std::dynamic_pointer_cast<RequestGPU>(req);
        auto &temperatures_vector = gpu_request->get_temperatures();
        fail_probability_previous = gpu_request->get_previous_state().failure_probability;

        if(fail_probability_previous == 0.0) fit_previous = getFitStandardGPU();    
        else {

            fit_previous = getUncompressFit((long double)this->lastProbGPU/1000,
											gpu_request->get_resource_type(),
											gpu_request->get_technology_type());
        }  
        
        for (int t = 0; t < temperatures_vector.size() - 1; t++)
        {
            /*Starting data*/
            second_chunk = (std::chrono::duration<double>(temperatures_vector[t + 1].epoch.time_since_epoch())).count() -
                           (std::chrono::duration<double>(temperatures_vector[t].epoch.time_since_epoch())).count();
            temperature_avg = (temperatures_vector[t].temperature + temperatures_vector[t + 1].temperature) / 2;
            acceleration_factor_ld = getAccelerationFactor(temperature_avg);

            fit_calculated = getFIT(fit_previous,this->fitOriginalGPU, getPastSecondGPU(), (double)gpu_request->get_activity() / 1000, second_chunk, acceleration_factor_ld);
            prob = getCompressProbability(fit_calculated, getPastSecondGPU() + second_chunk);
            this->lastProbGPU = prob * 1000; 
            fail_probability = (float)(prob * 1000);
            incrementPastSecond(second_chunk, gpu_request->get_resource_type(),gpu_request->get_technology_type());
            fit_previous = fit_calculated;
        }
    }else if(req->get_resource_type() == resource_type_t::MEMORY && req->get_technology_type()==technology_type_t::FPGA){

        std::shared_ptr<RequestAccelerator> accelerator_request = std::dynamic_pointer_cast<RequestAccelerator>(req);
        auto &temperatures_vector = accelerator_request->get_temperatures();
        fail_probability_previous = accelerator_request->get_previous_state().failure_probability;

        if(fail_probability_previous == 0.0) fit_previous = getFitStandardAccelerator();    
        else {

            fit_previous = getUncompressFit((long double) this->lastProbAccelerator/1000,
											accelerator_request->get_resource_type(),
											accelerator_request->get_technology_type());
        }
        for (int t = 0; t < temperatures_vector.size() - 1; t++)
        {
            /*Starting data*/
            second_chunk = (std::chrono::duration<double>(temperatures_vector[t + 1].epoch.time_since_epoch())).count() -
                           (std::chrono::duration<double>(temperatures_vector[t].epoch.time_since_epoch())).count();
            temperature_avg = (temperatures_vector[t].temperature + temperatures_vector[t + 1].temperature) / 2;
            acceleration_factor_ld = getAccelerationFactor(temperature_avg);

            fit_calculated = getFIT(fit_previous,this->fitOriginalAccelerator, getPastSecondAccelerator(),(double)accelerator_request->get_occupancy()/1000, second_chunk, acceleration_factor_ld);
            prob = getCompressProbability(fit_calculated, getPastSecondGPU() + second_chunk);
            this->lastProbAccelerator = prob * 1000;
            fail_probability = (float)(prob * 1000);
            incrementPastSecond(second_chunk,accelerator_request->get_resource_type(),accelerator_request->get_technology_type());
            fit_previous = fit_calculated;
        }  

    }
    return std::make_shared<Response>(fail_probability);
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

    long double f_rate_second = fit_original / (1000000000L * 3600L);

    unsigned long long t_standard = 1000000000L * 3600L;
    unsigned long long t_previous = (fit_previous / fit_original) * t_standard;

    unsigned long long t_add = time + (time * avf * af);
    unsigned long long t_left_std = t_standard - this->secondsCPU;

    unsigned long long t_calculated = (t_left_std + this->secondsCPU - time) + (t_previous - t_standard) + t_add;
    long double fit_new = f_rate_second * t_calculated;

    return fit_new;
}
long double BSC_HWReliabilityMonitor::getFitStandardCPU()
{
    return this->fitOriginalCPU;
}
long double BSC_HWReliabilityMonitor::getFitStandardMEM()
{
    return this->fitOriginalMEM;
}
long double BSC_HWReliabilityMonitor::getFitStandardGPU()
{
    return this->fitOriginalGPU;
}
long double BSC_HWReliabilityMonitor::getFitStandardAccelerator()
{
    return this->fitOriginalAccelerator;
}

unsigned long long BSC_HWReliabilityMonitor::getPastSecondCPU()
{
    return this->secondsCPU;
}
unsigned long long BSC_HWReliabilityMonitor::getPastSecondMEM()
{
    return this->secondsMEM;
}
unsigned long long BSC_HWReliabilityMonitor::getPastSecondGPU()
{
    return this->secondsGPU;
}
unsigned long long BSC_HWReliabilityMonitor::getPastSecondAccelerator()
{
    return this->secondsAccelerator;
}
void BSC_HWReliabilityMonitor::incrementPastSecond(unsigned long long value, resource_type_t res_type, technology_type_t tech_type)
{
    if( res_type == resource_type_t::CPU && tech_type == technology_type_t::SILICON){
        this->secondsCPU = this->secondsCPU + value;
    }else if(res_type == resource_type_t::MEMORY && tech_type == technology_type_t::SILICON ){
        this->secondsMEM = this->secondsMEM + value;
    }else if(res_type == resource_type_t::GPU && tech_type == technology_type_t::SILICON){
        this->secondsGPU = this->secondsGPU + value;
    }else if(res_type == resource_type_t::MEMORY && tech_type == technology_type_t::FPGA){
        this->secondsGPU = this->secondsAccelerator + value;
    }   

}
long double BSC_HWReliabilityMonitor::getCompressProbability(long double fit, unsigned long long second_passed)
{
    long double f_rate_second = fit / (1000000000L * 3600L);
    long double prob = 1 - powl((long double)(M_E), -1 * (f_rate_second * second_passed));

    return prob;
}
long double BSC_HWReliabilityMonitor::getUncompressFit(long double probability, resource_type_t res_type, technology_type_t tech_type)
{
    

    long double base = (long double)1.0 - probability;
    long double espo;
    if(res_type == resource_type_t::CPU && tech_type == technology_type_t::SILICON) 
        espo = ((long double)1.0 / (-1.0 * (long double )getPastSecondCPU()));
    else if(res_type == resource_type_t::MEMORY && tech_type == technology_type_t::SILICON)    
        espo = ((long double)1.0 / (-1.0 * (long double )getPastSecondMEM()));
    else if(res_type == resource_type_t::GPU && tech_type == technology_type_t::SILICON)       
        espo = ((long double)1.0 / (-1.0 * (long double )getPastSecondGPU()));
    else if(res_type == resource_type_t::MEMORY && tech_type == technology_type_t::FPGA)
        espo = ((long double)1.0 / (-1.0 * (long double )getPastSecondAccelerator()));

    long double base_log = powl(base, espo);
    long double fit_rate_second = logl(base_log);
    long double fit_std = fit_rate_second * 1000000000L * 3600L;

    return fit_std;
}

int BSC_HWReliabilityMonitor::getSKT_number(std::vector<core> cv)
{
    std::vector<int> tmp;
    bool flag = false;
    
    for(int core = 0; core < cv.size(); core++){
        flag = false;
        for ( int j=0; j<tmp.size(); j++){
            if(cv[core].id_socket == tmp[j]) flag  = true;
        }
        if(!flag)tmp.push_back(cv[core].id_socket);
   } 
    return tmp.size();
}

int BSC_HWReliabilityMonitor::getPhyCoreSize(int skt, std::vector<core> &cv)
{
    std::vector<int> tmp;

    for(int c=0; c<cv.size();c++){
        if(cv[c].id_socket == skt){
            bool flag = false;
            for(int j=0;j<tmp.size();j++){
                if(tmp[j] == cv[c].id_physical)flag=true;
            }
            if(!flag)tmp.push_back(cv[c].id_physical);
        }
    }
    return tmp.size();
}

int BSC_HWReliabilityMonitor::getSMT(int skt, std::vector<core> &cv)
{
    int l = 0;
    for(int c=0; c<cv.size();c++){
        if(cv[c].id_socket == skt){
            l++;
        }
    }
    return l/getPhyCoreSize(skt,cv);
}

void BSC_HWReliabilityMonitor::removeSMT(std::vector<core> &cv)
{

    int socket = getSKT_number(cv);

    std::vector<std::vector<core>> mat;
    std::vector<bool> available(cv.size(),true);
    std::vector<int> socketVector = getSocketVector(cv);
    
    /**
     *  Building a matrix where each row is composed
     *  by logical core. Two core logical in the same
     *  physical core reside in diffrent rows 
     */
     
    for(int skt = 0 ; skt < socketVector.size() ; skt++){
        int level = getSMT(socketVector[skt],cv);
        for(int lev=0;lev<level;lev++){
            std::vector<core> tmp;
            for(int log=0;log<cv.size();log++){
                if(socketVector[skt]==cv[log].id_socket && available[log]==true){
                    bool present=false;
                    for(int i=0;i<tmp.size() && !present;i++){
                        if(tmp[i].id_physical == cv[log].id_physical)present = true;
                    }
                    if(!present){
                        tmp.push_back(cv[log]);
                        available[log]=false;
                    }
                }
            }
            mat.push_back(tmp);
        }
    }

    std::vector<core> tempo;


    for(int r=0;r<mat.size();){
        int smt = getSMT(mat[r][0].id_socket, cv);
        int phy_size = getPhyCoreSize(mat[r][0].id_socket,cv);
        /*
            Should be fixed counter 
            for all logical core
            All counter must have
            the same order
        */

        /*
            matrix 
            1 dimension = size counter used of a general core, enum size
            2 dimension = core physical size of a socket
        */
        std::vector<std::vector<unsigned long long>> tmp(mat[0][0].counter.size(),std::vector<unsigned long long>(phy_size,0));

        for(int c=0;c<mat[r].size();c++){ /*all core in the current row*/
            for(int counter = 0; counter<tmp.size();counter++){ /*checking all counter of actual core*/
         
                tmp[counter][c]+=mat[r][c].counter[counter].get_value();

                for(int l=r+1; l<(r+1+smt-1);l++){
                    for(int k=0;k<mat[l].size();k++){
                        if(mat[l][k].id_physical == mat[r][c].id_physical){
                            /*
                                if we check cycles counter we get 
                                the maxinum value of two logical core 
                                in the same physical core
                                for the other 
                            */
                            if(mat[l][k].counter[counter].get_type() == perf_counter_type_t::CYCLES){
                                if(mat[l][k].counter[counter].get_value() >  tmp[counter][c]){
                                    tmp[counter][c]=mat[l][k].counter[counter].get_value();
                            }
                            }else{
                                tmp[counter][c]+=mat[l][k].counter[counter].get_value();
                            }

                        }
                    }
                }
               

            }
            std::vector<PerfCounter> counter_vector;
            for(int perf = 0; perf < tmp.size(); perf++){
                counter_vector.push_back(PerfCounter(mat[r][c].counter[perf].get_type(),tmp[perf][c]));
            }
            unsigned int id_logical = (unsigned int)mat[r][c].id_physical;
            unsigned int id_physical = (unsigned int)mat[r][c].id_physical;
            unsigned int id_socket = (unsigned int)mat[r][c].id_socket;
            tempo.push_back({id_logical, id_physical, id_socket, counter_vector});

        }
    r+=smt;
    }

    cv.clear();
    std::copy(tempo.begin(), tempo.end(), back_inserter(cv));    
}

unsigned long long BSC_HWReliabilityMonitor::getALU(std::vector<PerfCounter> v){
    unsigned long long alu = 0;
    for(int i=0;i<v.size();i++){
        if(v[i].get_type() == perf_counter_type_t::PORT_0 ||
            v[i].get_type() == perf_counter_type_t::PORT_1 ||
            v[i].get_type() == perf_counter_type_t::PORT_5 ||
            v[i].get_type() == perf_counter_type_t::PORT_6){
  
            alu = alu + v[i].get_value();
        }
    }
return alu;
}

unsigned long long BSC_HWReliabilityMonitor::getLD(std::vector<PerfCounter> v){
    unsigned long long ld = 0;
    for(int i=0;i<v.size();i++){
        if(v[i].get_type() == perf_counter_type_t::PORT_2 ||
            v[i].get_type() == perf_counter_type_t::PORT_3 ){
  
            ld = ld + v[i].get_value();
        }
    }
return ld;
}

unsigned long long BSC_HWReliabilityMonitor::getST(std::vector<PerfCounter> v){
    unsigned long long st = 0;
    for(int i=0;i<v.size();i++){
        if(v[i].get_type() == perf_counter_type_t::PORT_4){
            st = st + v[i].get_value();
        }
    }
return st;
}

unsigned long long BSC_HWReliabilityMonitor::getAGU(std::vector<PerfCounter> v){
    unsigned long long agu = 0;
    for(int i=0;i<v.size();i++){
        if(v[i].get_type() == perf_counter_type_t::PORT_7){
            agu = agu + v[i].get_value();
        }
    }
return agu;
}

unsigned long long BSC_HWReliabilityMonitor::getL1IH(std::vector<PerfCounter> v){
    unsigned long long value = 0;
    for(int i=0;i<v.size();i++){
        if(v[i].get_type() == perf_counter_type_t::L1I_HIT){
            value = value + v[i].get_value();
        }
    }
return value;
}

unsigned long long BSC_HWReliabilityMonitor::getL1IM(std::vector<PerfCounter> v){
    unsigned long long value = 0;
    for(int i=0;i<v.size();i++){
        if(v[i].get_type() == perf_counter_type_t::L1I_MISS){
            value = value + v[i].get_value();
        }
    }
return value;
}

unsigned long long BSC_HWReliabilityMonitor::getL2M(std::vector<PerfCounter> v){
    unsigned long long value = 0;
    for(int i=0;i<v.size();i++){
        if(v[i].get_type() == perf_counter_type_t::L2_RQST_MISS){
            value = value + v[i].get_value();
        }
    }
return value;
}

unsigned long long BSC_HWReliabilityMonitor::getL2H(std::vector<PerfCounter> v){
    unsigned long long value = 0;
    unsigned long long value2 = 0 ;

    for(int i=0;i<v.size();i++){
        if(v[i].get_type() == perf_counter_type_t::L2_RQST_MISS){
            value = value + v[i].get_value();
        }
        if(v[i].get_type() == perf_counter_type_t::L2_RQST_REFERENCES){
            value2 = value2 + v[i].get_value();
        }

    }
return (value2 - value);
}
    
unsigned long long BSC_HWReliabilityMonitor::getL3M(std::vector<PerfCounter> v){
    unsigned long long value = 0;

    for(int i=0;i<v.size();i++){
        if(v[i].get_type() == perf_counter_type_t::L3_MISS){
            value = value + v[i].get_value();
        }
    }
return value;
}   

unsigned long long BSC_HWReliabilityMonitor::getL3H(std::vector<PerfCounter> v){
    unsigned long long value = 0;
    unsigned long long value2 = 0;
    
    for(int i=0;i<v.size();i++){
        if(v[i].get_type() == perf_counter_type_t::L3_MISS){
            value = value + v[i].get_value();
        }
        if(v[i].get_type() == perf_counter_type_t::L3_REFERENCES){
            value2 = value2 + v[i].get_value();
        }
    }
return (value2 - value);
}
unsigned long long BSC_HWReliabilityMonitor::getCyclesMax(std::vector<core> v){

    unsigned long long cycles = 0;
    for(int i=0;i<v[0].counter.size();i++){
        if(v[0].counter[i].get_type() == perf_counter_type_t::CYCLES){
            cycles = (unsigned long long)v[0].counter[i].get_value();
        }
    }
    for(int i=1;i<v.size();i++){
        for(int c=0;c<v[i].counter.size();c++){
            if(v[i].counter[i].get_type() == perf_counter_type_t::CYCLES && v[i].counter[i].get_value() > cycles ){
                cycles = (unsigned long long)v[i].counter[i].get_value();
            }   
        }        
    }

    return cycles;
}

unsigned long long BSC_HWReliabilityMonitor::getUopsIssued(std::vector<PerfCounter> v){
    unsigned long long value = 0;
    
    for(int i=0;i<v.size();i++){
        if(v[i].get_type() == perf_counter_type_t::UOPS_ISSUED_ANY){
            value = value + v[i].get_value();
        }
    }
return value;
}

unsigned long long BSC_HWReliabilityMonitor::getUopsRetired(std::vector<PerfCounter> v){
    unsigned long long value = 0;
    
    for(int i=0;i<v.size();i++){
        if(v[i].get_type() == perf_counter_type_t::UOPS_RETIRED){
            value = value + v[i].get_value();
        } 
    }
return value;
}

void BSC_HWReliabilityMonitor::init_FIT_original_block(int cores_physical_number, int socket_number)
{

    this->original_fitALU = this->fitOriginalCPU / ( this->bl_size * (cores_physical_number)/socket_number) ;
    this->original_fitLD_ST_AGU = this->fitOriginalCPU / ( this->bl_size * (cores_physical_number)/socket_number) ;
    this->original_fitROB = this->fitOriginalCPU / ( this->bl_size * (cores_physical_number)/socket_number) ;
    this->original_fitL1I = this->fitOriginalCPU / ( this->bl_size * (cores_physical_number)/socket_number) ;
    this->original_fitL1D = this->fitOriginalCPU / ( this->bl_size * (cores_physical_number)/socket_number)  ;
    this->original_fitL2 = this->fitOriginalCPU / ( this->bl_size * (cores_physical_number)/socket_number) ;
    this->original_fitL3 = this->fitOriginalCPU / ( this->bl_size * (cores_physical_number)/socket_number) ;

}

void BSC_HWReliabilityMonitor::init_FIT_last_block( int size)
{
    this->last_fitALU.clear();
    this->last_fitLD_ST_AGU.clear();
    this->last_fitROB.clear();
    this->last_fitL1I.clear();
    this->last_fitL1D.clear();
    this->last_fitL2.clear();
    this->last_fitL3.clear();


    for(int i=0;i<size;i++){
    /*
        Set in the first call of perfom calculation 
        last fit block for every block. We use the
        default value of original FIT per block.
    */  
        /*                    Core Calculation Units                    */
        this->last_fitALU.push_back(this->original_fitALU);
        this->last_fitLD_ST_AGU.push_back(this->original_fitLD_ST_AGU);
        /*                    Rob                    */
        this->last_fitROB.push_back(this->original_fitROB);
        /*                    Caches                    */
        this->last_fitL1I.push_back(this->original_fitL1I);
        this->last_fitL1D.push_back(this->original_fitL1D);
        this->last_fitL2.push_back(this->original_fitL2);
        this->last_fitL3.push_back(this->original_fitL3);
}
}

void BSC_HWReliabilityMonitor::stamp(std::vector<core> c){

    for(int core=0;core<c.size();core++){
        std::cout << "Core " << c[core].id_logical << " " << c[core].id_physical << " " << c[core].id_socket << " " ;
        for(int counter = 0; counter < c[core].counter.size(); counter++){
            std::cout << "*" << (int)c[core].counter[counter].get_type() << " " << c[core].counter[counter].get_value() << "" ;
        } 
        std::cout << std::endl;
    }
}

void BSC_HWReliabilityMonitor::reset(){
    this->secondsCPU = 0L;    
    this->secondsMEM = 0L;
    this->secondsGPU = 0L;
    this->secondsAccelerator = 0L;
    
    this->lastProbCPU = 0;
    this->lastProbMEM = 0;
    this->lastProbGPU = 0;
    this->lastProbAccelerator = 0;
}



void BSC_HWReliabilityMonitor::checkPast(std::shared_ptr<Request> req){

    for(int i=0;i<req->get_temperatures().size();i++){
        if(req->get_previous_state().epoch >  req->get_temperatures()[i].epoch){
           throw std::invalid_argument("Request not subsequent in time.");
        }
    }
}

void BSC_HWReliabilityMonitor::checkTimeRequest(std::shared_ptr<Request> req){
    for(int i=0;i<req->get_temperatures().size()-1;i++){
        if(req->get_temperatures()[i].epoch > req->get_temperatures()[i+1].epoch){
           throw std::invalid_argument("Request vector not subsequent in time.");
        }
    }
}
void BSC_HWReliabilityMonitor::checkTemperatureSize(std::shared_ptr<Request> req){

    if(req->get_temperatures().size()<=1){
        throw std::invalid_argument("Request vector size error <= 1");
    }
}

std::vector<int> BSC_HWReliabilityMonitor::getSocketVector(std::vector<core> v){
    std::vector<int> socketVector ; 
    bool flag = false;

    for(int log=0;log<v.size();log++){
        flag = false;
        for(int i=0;i<socketVector.size() && !flag;i++){
            if(socketVector[i] == v[log].id_socket)flag = true;
        }
        if(!flag)socketVector.push_back(v[log].id_socket);
    }

    return socketVector;
}

long double BSC_HWReliabilityMonitor::getFIT_auxiliar()
{
    return this->lastProbCPU;
}
// You can also cast "req" to td::shared_ptr<RequestCPU> or td::shared_ptr<RequestMEM> or ...

long double BSC_HWReliabilityMonitor::getFit_processing(std::vector<core> cv, unsigned long long  second_chunk, float temperature_avg, long double acceleration_factor_ld )
{
    unsigned long long  cycles = getCyclesMax(cv);
    long double fit_processing = 0 ;
    /*Activiti block variable*/
    long double act_alu;
    long double act_ld_st_agu;
    long double ratio_issued;
    long double ratio_retired;
    long double act_rob;
    long double act_L1I;
    long double act_L1D;
    long double act_L2;
    long double act_L3;

     for(int p_core=0;p_core<cv.size();p_core++){

                    /*
                        activity calculation
                    */
                    /*alu, ld_st_agu*/
                    act_alu = (long double)getALU(cv[p_core].counter) / (cycles * 4);
                    act_ld_st_agu = (long double)(getLD(cv[p_core].counter) + getST(cv[p_core].counter) + getAGU(cv[p_core].counter) ) / (cycles * 4);
                    
                    /*ROB*/
                    ratio_issued = (long double)getUopsIssued(cv[p_core].counter) / (4 * cycles);
			        ratio_retired =  (long double)getUopsRetired(cv[p_core].counter) / (4 * cycles);
                    act_rob = ratio_retired / ratio_issued;

                    /*Caches*/
                    act_L1I = (long double)getL1IH(cv[p_core].counter) / cycles;
                    act_L1I +=(long double)getL1IM(cv[p_core].counter) / cycles;
                    act_L1D = (long double)getL1IH(cv[p_core].counter) / cycles;
                    act_L1D +=(long double)getL1IM(cv[p_core].counter) / cycles;
                    act_L2 = (long double)(getLD(cv[p_core].counter) + getST(cv[p_core].counter)) / (cycles * 3);
                    act_L3 = (long double)getL3H(cv[p_core].counter) / cycles;
                    act_L3 += (long double)getL3M(cv[p_core].counter) / cycles;

                /*
                    Set in the first call of perfom calculation 
                    last fit block for every block. We use the
                    default value of original FIT per block.
                */  
                    /*
                    Core Calculation Units
                    */
                    this->last_fitALU[p_core]=getFIT(this->last_fitALU[p_core],this->original_fitALU, getPastSecondCPU(), act_alu , second_chunk, acceleration_factor_ld);
                    this->last_fitLD_ST_AGU[p_core]=getFIT(this->last_fitLD_ST_AGU[p_core],this->original_fitLD_ST_AGU, getPastSecondCPU(), act_ld_st_agu, second_chunk, acceleration_factor_ld);
                    /*
                    Rob
                    */
                    this->last_fitROB[p_core] = getFIT(this->last_fitROB[p_core],this->original_fitROB,getPastSecondCPU(), act_rob, second_chunk, acceleration_factor_ld);
                    /*
                    Caches
                    */
                    this->last_fitL1I[p_core] = getFIT(this->last_fitL1I[p_core],this->original_fitL1I ,getPastSecondCPU(),act_L1I, second_chunk, acceleration_factor_ld);
                    this->last_fitL1D[p_core] = getFIT(this->last_fitL1I[p_core],this->original_fitL1D ,getPastSecondCPU(),act_L1D, second_chunk, acceleration_factor_ld);
                    this->last_fitL2[p_core] = getFIT(this->last_fitL2[p_core],this->original_fitL2, getPastSecondCPU(),act_L2, second_chunk, acceleration_factor_ld);
                    this->last_fitL3[p_core] = getFIT(this->last_fitL3[p_core],this->original_fitL3, getPastSecondCPU(),act_L3, second_chunk, acceleration_factor_ld);
                    
                    
                    fit_processing+=this->last_fitALU[p_core];
                    fit_processing+=this->last_fitLD_ST_AGU[p_core];
                    fit_processing+=this->last_fitROB[p_core];
                    fit_processing+=this->last_fitL1I[p_core];
                    fit_processing+=this->last_fitL1D[p_core];
                    fit_processing+=this->last_fitL2[p_core];
                    fit_processing+=this->last_fitL3[p_core];
                    
            }

            return fit_processing;
    
}

}
