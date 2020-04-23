#ifndef PERF_COUNTER_H_
#define PERF_COUNTER_H_

namespace libhwrel {

enum class perf_counter_type_t {
    CPU_L1_MISS,
    CPU_L1_HIT,
    CPU_L2_MISS,
    CPU_L2_HIT,
    
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
