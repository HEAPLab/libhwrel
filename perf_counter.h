#ifndef PERF_COUNTER_H_
#define PERF_COUNTER_H_

namespace libhwrel {

enum class perf_counter_type_t {
    
    /**
     *  Need to be aggregate by core physical
     * 
     */

    L1I_HIT, //0x8301
    L1I_MISS, //0x8302
    L2_RQST_MISS , //x0243f
    L2_RQST_REFERENCES, //0x24ff
    L3_MISS, //0x2e41
    L3_REFERENCES, //0X2E4F

    PORT_0, //0xA101
    PORT_1, //0xA102
    PORT_2, //0xA104
    PORT_3, //0xA108
    PORT_4, //0xA110
    PORT_5, //0xA120
    PORT_6, //0xA140
    PORT_7,  //0xA180

    UOPS_RETIRED, //0XC202
    UOPS_ISSUED_ANY, //0X0E01

    /*
    * Need to be max of cpu logical 0 
    * and cpu logical 1. Elapsed since time 0
    * to time 1. 
    */

    CYCLES, 


    /*
    *   Should be the sum of all values CAS_READ and
    *   CAS_WRITE for a whole socket.
    */
    CAS_READ,   //event = 0x4 umask=0x3  
    CAS_WRITE,   //event = 0x4 umask=0xC

    /*
        For temperature counter
        value = counter_bits[22:16]
    */
    TEMPERATURE_SKT, //0x01B1  
    TEMPERATURE_CORE, //0x019C  


    /*
    * 
    * General Purpose Graphics Processing Unit 
    * GPU
    *  The values are from 0 to 10 
    *  issue_slot_utilization need to be converted to same format
    */
    single_precision_fu_utilization,
    half_precision_fu_utlization,
    double_precision_fu_utilization,
    ldst_fu_utilization,
    special_fu_utilization,
    issue_slot_utilization , 
    tex_utilization,
    shared_utilization,
    l2_utilization,
    dram_utilization
    
};


/** 
 * @brief The class representing a performance counter.
 * This class is built by the resource manager and provided to the library
 * @note This class is **immutable**.
 */
class PerfCounter {

public:

    /**
     * @brief The class constructor
      */
    PerfCounter(perf_counter_type_t type, unsigned long value) noexcept : type(type), value(value)
    {
    }
    
    /** @brief Getter for the PC value */
    unsigned long get_value() const noexcept {
        return this->value;
    }

    /** @brief Getter for the PC type */
    perf_counter_type_t get_type() const noexcept {
        return this->type;
    }

private:
    const perf_counter_type_t type;
    const unsigned long value;
};

} // namespace libhwrel

#endif
