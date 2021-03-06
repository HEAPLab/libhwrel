/** @file libhwrel.h
 * The main header file for libhwrel.
 *
 * This header contains the main classes to be used for interfacing the resource manager with the
 * underlying software that computes the hardware reliability information.
 */

#ifndef LIBHWREL_H_
#define LIBHWREL_H_

#include "perf_counter.h"

#include <chrono>
#include <map>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>
#include <iostream>

namespace libhwrel {

//
// ---------------------------------- ENUMERATIONS ---------------------------------- 
//

/**
 * @brief Enumeration for resource types.
 *
 * This enum represents the possible values for the type of resource.
 */
enum class resource_type_t {
    CPU,        ///< General Purpose CPU
    GPU,        ///< GPGPU
    MEMORY_GPU,    ///< GPGPU MEMORY
    MEMORY      ///< Memory (usually DRAM)
};

/** 
 * @brief Enumeration for technology types.
 *
 * This enum represents the possible values for the technology used to build the specific resource.
 */
enum class technology_type_t {
    SILICON,    ///< Custom ASIC chip
    FPGA        ///< FPGA-implemented resource
};

//
// ------------------------------------ STRUCTS ------------------------------------ 
//

/** 
 * @brief The struct containing the failure probability previously provided by the the HW monitor.
 * 
 * This struct contains the failure probability previously returned by the HW reliability monitor,
 * together with the time point when it was generated.
 */
typedef struct reliability_state_t {
    std::chrono::time_point<std::chrono::system_clock> epoch;
    long double failure_probability;          /**< The failure probability CDF on 1000% */
    
    std::shared_ptr<void> state;
    size_t state_size;
} reliability_state_t;

/**
 * @brief The struct containing a single temperature element. 
 *
 * This struct contains a single temperature with the associated time point. This is used to
 * provide to the hw monitor the list of previous and current values of temperatures. 
 */
struct temperature_t {
    float temperature;                              /** The temperature in Celsius */
    std::chrono::time_point<std::chrono::system_clock> epoch;
};

//
// ------------------------------------ CLASSES ------------------------------------ 
//


/** 
 * @brief The class representing the object returned by the HW reliability monitor.
 *
 * This class is constructed and returned by the HW reliability monitor as a response of a subclass
 * of Request class.
 * @note This class is **immutable**.
 */
class Response {

public:
    /**
     * @brief The class constructor
     * @param state The state (containing the failure probability)
     * @throw std::invalid_argument If state.failure_probability < 0 or state.failure_probability > 1000.
      */


    Response(const reliability_state_t &state) : state(state) {
        if(state.failure_probability > 1000.0 || state.failure_probability < 0.0) {
            throw std::invalid_argument("Failure probability invalid value (>1000 or <0).");
        }
    }


    /** @brief Getter for the failure probability */
    float get_fail_probability() const noexcept {
        return this->state.failure_probability;
    }

    /** @brief Getter for the state */
    reliability_state_t get_state() const noexcept {
        return this->state;
    }
    

private:
    reliability_state_t state;

};

/** 
 * @brief The (parent) class representing the object sended to the HW reliability monitor by the
 *        resource manager.
 *
 * This class is constructed and sended to the HW reliability monitor by the resource manager to ask
 * for a computation of the failure probability. The objects beloning to this class are usually
 * initialized with a specific subclass depending on resource types.
 */
class Request {

public:
    /**
     * @brief The Request class constructor
     * @param res_type  The type of resource
     * @param tech_type The type of technology
      */
    Request(resource_type_t res_type, technology_type_t tech_type) noexcept
        : res_type(res_type), tech_type(tech_type)  {

    }

    /**
     * @brief The default virtual destructor (no dynamic memory used)
      */
    virtual ~Request() = default;

    /** @brief Getter for resource type */
    inline resource_type_t get_resource_type() const noexcept {
        return this->res_type;
    }

    /** @brief Getter for technology type */
    inline technology_type_t get_technology_type() const noexcept {
        return this->tech_type;
    }

    /** 
          * @brief Getter for the previous temperatures array. 
      * @note The caller must not try to edit the array returned by this function directly
          */
    const std::vector<temperature_t> & get_temperatures() const noexcept
    {
        return this->temperatures;
    }

    /** @brief Getter for the previous state returned by the HW reliability monitor. If the `epoch`
      *        of the previous state is 0, it should be interpret as no previous state exists. */
    reliability_state_t get_state() const noexcept {
        return this->prev_state;
    }

    /** @brief Setter for the previous state returned by the HW reliability monitor */
    void set_state(reliability_state_t prev_state) noexcept {
        this->prev_state = prev_state;
    }

    /** 
     * @brief Add a new temperature to the array. This function should be called only on RM-side 
     * @note This should be used by the Resource Manager only!
     */
    void push_temperature(const temperature_t& temp) {
        this->temperatures.push_back(temp);
    }

    /** @brief Clear the whole temperature array. This function should be called only on RM-side */
    void clear_temperatures() noexcept {
        this->temperatures.clear();
    }

    /** @brief Return the value of a performance counter by type
     *  @throws std::out_of_range if not found
     */
    const PerfCounter &get_PC(perf_counter_type_t type) const {
        return perf_counters.at(type);
    }

    /** @brief Add a new value for a given performance number
      * @note This should be used by the Resource Manager only!
     */
    void add_PC(perf_counter_type_t type, unsigned long value) noexcept {
        perf_counters.emplace(type, PerfCounter(type, value));
    }


private:
    const resource_type_t   res_type;
    const technology_type_t tech_type;

    std::map<perf_counter_type_t, PerfCounter> perf_counters;

    std::vector<temperature_t> temperatures;

    reliability_state_t prev_state;

};

/** 
 * @brief The specialied Request class for memories.
 *
 */
class RequestMEM : public Request {
public:

    /**
     * @brief The RequestMEM class constructor
     * @param tech_type The type of technology
     * @param band_skt_max  maximum value of bandwith supported by the CPU SKT, supposing
     *                      the DIMM provided max perfomance. If property maxinum bandwithc
     *                      of DIMM are available set this. But CPU SKT vision. [B/s]
     * 
     * 
     * @throw std::invalid_argument If occupancy > 1000.
      */
    RequestMEM(technology_type_t tech_type, unsigned long long  band_skt_max)
    : Request(resource_type_t::MEMORY, tech_type), band_skt_max(band_skt_max)
    {

    }

    /**
     * @brief The default virtual destructor (no dynamic memory used)
      */
    virtual ~RequestMEM() = default;


    inline unsigned long long get_band_per_skt() const noexcept {
        return this->band_skt_max;
    }



private:
    unsigned long long band_skt_max;        /* Max bandwith per socket*/

};


/** 
 * @brief The specialied Request class for memories gpu.
 *
 */
class RequestMEM_GPU : public Request {
public:

    /**
     * @brief The RequestMEM class constructor
     * @param tech_type The type of technology
     * @param band_max
     * @param runtime_mode If true, it uses the approximate mode (by considering only the MEMORY_UTILIZATION
     *                     counter). Otherwise, it uses the whole set of counters.
     * 
     * 
     * @throw std::invalid_argument If occupancy > 1000.
      */
    RequestMEM_GPU(technology_type_t tech_type, unsigned long long  band_max, bool runtime_mode)
    : Request(resource_type_t::MEMORY_GPU, tech_type), band_max(band_max), runtime_mode(runtime_mode)
    {

    }

    /**
     * @brief The default virtual destructor (no dynamic memory used)
      */
    virtual ~RequestMEM_GPU() = default;


    inline unsigned long long get_band_per_gpu() const noexcept {
        return this->band_max;
    }

    inline bool get_runtime_mode() const {
        return this->runtime_mode;
    }

private:
    unsigned long long band_max;        /* Max bandwith per gpu*/
    bool runtime_mode;
};

/** 
 * @brief The specialied Request class for *one single physical core* of the CPU.
 *
 */
class RequestCPUCore : public Request {
public:

    /**
     * @brief The RequestCPU class constructor
     * @param tech_type The type of technology
     * @param clock_frequency The clock frequency in MHz
      */
    RequestCPUCore(technology_type_t tech_type , unsigned int clock_frequency) noexcept
    : Request(resource_type_t::CPU, tech_type)
    {
        this->clock_frequency = clock_frequency ; 
    }
    
    void set_clock_frequency( unsigned int clock_frequency){
        this->clock_frequency = clock_frequency; 
    }

    inline unsigned int get_clock_frequency() const noexcept {
        return this->clock_frequency ; 
    }
    /**
     * @brief The default virtual destructor (no dynamic memory used)
      */
    virtual ~RequestCPUCore() = default;


private:
    unsigned int clock_frequency;

};

/** 
 * @brief The specialied Request class for CPUs.
 *
 */
class RequestGPU : public Request {
public:

    /**
     * @brief The RequestCPU class constructor
     * @param tech_type The type of technology
     * @param clock_frequency The clock frequency in MHz
     * @param runtime_mode If true, it uses the approximate mode (by considering only the GPU_UTILIZATION
     *                     counter). Otherwise, it uses the whole set of counters  
     * @throw std::invalid_argument If activity > 1000.
      */
    RequestGPU(technology_type_t tech_type, unsigned int clock_frequency, bool runtime_mode)
    : Request(resource_type_t::GPU, tech_type), clock_frequency(clock_frequency), runtime_mode(runtime_mode)
    {
        // if(activity > 1000) {
        //     throw std::invalid_argument("Activity invalid value (>1).");
        // }
    }

    /**
     * @brief The default virtual destructor (no dynamic memory used)
      */
    virtual ~RequestGPU() = default;

    /** @brief Getter for the clock frequency in MHz. */
    inline unsigned int get_clock_frequency() const {
        return this->clock_frequency;
    }

    // inline unsigned short get_activity(){
    //     return this->activity;
    // }
    /** 
     * @brief Setter for the clock frequency in MHz. 
     * @note This method should be used by the Resource Manager only!
     */
    inline void set_clock_frequency(unsigned int clock_frequency) {
        this->clock_frequency = clock_frequency;
    }
    
    inline bool get_runtime_mode() const {
        return this->runtime_mode;
    }


private:
    unsigned int clock_frequency;   ///< The clock frequency in MHz
    bool runtime_mode;
};


/** 
 * @brief The specialied Request class for memories.
 *
 */
class RequestAccelerator : public Request {
public:

    /**
     * @brief The RequestFPGA class constructor
     * @param tech_type The type of technology
     * @param size      The size in MB
     * @param occupancy The level of occupancy for the memory in per-mille format (0-1000). The
     *                  value 1000 means 100%.
     * @throw std::invalid_argument If occupancy > 1000.
      */
    RequestAccelerator(technology_type_t tech_type, unsigned int occupancy)
    : Request(resource_type_t::MEMORY, tech_type), occupancy(occupancy) 
    {
        if(occupancy > 1000) {
            throw std::invalid_argument("Occupancy invalid value (>1).");
        }
    }

    /**
     * @brief The default virtual destructor (no dynamic memory used)
      */
    virtual ~RequestAccelerator() = default;

    /** @brief Getter for the occupancy. The level of occupancy for the memory in per-mille format (0-1000).  */
    inline unsigned short get_occupancy() const noexcept {
        return this->occupancy;
    }

    /**
     * @brief Setter for the occupancy. The level of occupancy for the memory in per-mille format (0-1000).
     * @note This should be used by the Resource Manager only!
     */
    inline void set_occupancy(unsigned int occupancy) noexcept {
        this->occupancy = occupancy;
    }

private:
    unsigned short occupancy;    // 0-1000
};


/**
 * @brief The main class to be inherited and implemented by the HW reliability monitor.
 *
 * @note This class is abstract and must be inherited to be instanced.
 */
class HWReliabilityMonitor {

public:

    /** @brief The default virtual destructor. */
    virtual ~HWReliabilityMonitor() = default;

    /**
     * @brief The main method to be implemented that will perform the analysis.
     */
    virtual std::shared_ptr<Response> perform_analysis(std::shared_ptr<Request> req) = 0;

    /**
     * @brief Initialize the reliability monitor and return the initial state
     */
    virtual reliability_state_t init(resource_type_t resource, technology_type_t tech, long double initial_fit, unsigned int nr_cores=0) = 0;

};

}    // namespace libhwrel

#endif // LIBHWREL_H_

