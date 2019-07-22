/** @file libhwrel.h
 * The main header file for libhwrel.
 *
 * This header contains the main classes to be used for interfacing the resource manager with the
 * underlying software that computes the hardware reliability information.
 */

#ifndef LIBHWREL_H_
#define LIBHWREL_H_

#include <chrono>
#include <stdexcept>
#include <utility>
#include <vector>

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
	CPU,		///< General Purpose CPU
	GPU,		///< GPGPU
	MEMORY		///< Memory (usually DRAM)
};

/** 
 * @brief Enumeration for technology types.
 *
 * This enum represents the possible values for the technology used to build the specific resource.
 */
enum class technology_type_t {
	SILICON,	/**< Custom ASIC chip */
	FPGA		/**< FPGA-implemented resource */
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
struct previous_reliability_t {
	unsigned short failure_probability;		  /**< The failure probability in
							                       per-mille format (0-1000) */
	std::chrono::time_point<std::chrono::system_clock> epoch;
};

/**
 * @brief The struct containing a single temperature element. 
 *
 * This struct contains a single temperature with the associated time point. This is used to
 * provide to the hw monitor the list of previous and current values of temperatures. 
 */
struct temperature_t {
	float temperature;	                          /** The temperature in Celsius */
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
	 * @param failure_probability The failure probability in 0-1000 per-mille format (e.g.
     *                            10 means 0.001).
	 * @throw std::invalid_argument If failure_probability < 0 or failure_probability > 1000.
 	 */
	explicit Response(unsigned short failure_probability) : fail_probability(fail_probability) {
		if(failure_probability > 1000) {
			throw std::invalid_argument("Failure probability invalid value (>1).");
		}
	}

	/** @brief Getter for the failure probability */
	inline unsigned short get_fail_probability() const noexcept {
		return this->fail_probability;
	}

private:
	const unsigned short fail_probability;	// 0-1000 in per mille

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
	explicit Request(resource_type_t res_type, technology_type_t tech_type) noexcept
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

	/** @brief Getter for the previous state returned by the HW reliability monitor */
	previous_reliability_t get_previous_state() const noexcept {
		return this->prev_state;
	}

	/** @brief Add a new temperature to the array. This function should be called only on RM-side */
	void push_temperature(const temperature_t& temp) {
		this->temperatures.push_back(temp);
	}

	/** @brief Clear the whole temperature array. This function should be called only on RM-side */
	void clear_temperatures() {
		this->temperatures.clear();
	}

private:
	const resource_type_t   res_type;
	const technology_type_t tech_type;

	std::vector<temperature_t> temperatures;

	previous_reliability_t prev_state;

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
	 * @param size      The size in MB
	 * @param occupancy The level of occupancy for the memory in per-mille format (0-1000). The
     *                  value 1000 means 100%.
	 * @throw std::invalid_argument If occupancy > 1000.
 	 */
	RequestMEM(technology_type_t tech_type, unsigned int size, unsigned int occupancy)
	: Request(resource_type_t::MEMORY, tech_type), size(size), occupancy(occupancy) 
	{
		if(occupancy > 1000) {
			throw std::invalid_argument("Occupancy invalid value (>1).");
		}
	}

	/**
	 * @brief The default virtual destructor (no dynamic memory used)
 	 */
	virtual ~RequestMEM() = default;

	/** @brief Getter for the size in MB. This value is always constant. */
	inline unsigned int get_size() const noexcept {
		return this->size;
	}

	/** @brief Getter for the occupancy. The level of occupancy for the memory in per-mille format (0-1000).  */
	inline unsigned short get_occupancy() const noexcept {
		return this->occupancy;
	}

	/** @brief Setter for the occupancy. The level of occupancy for the memory in per-mille format (0-1000).  */
	inline void set_occupancy(unsigned int occupancy) noexcept {
		this->occupancy = occupancy;
	}

private:
	const unsigned int size; 	// in MiB
	unsigned short occupancy;	// 0-1000
};

/** 
 * @brief The specialied Request class for CPUs.
 *
 */
class RequestCPU : public Request {
public:

	/**
	 * @brief The RequestCPU class constructor
	 * @param tech_type The type of technology
	 * @param clock_frequency The clock frequency in MHz
	 * @param nr_cores The number of physical processing elements.
	 * @param activity The activity level of the whole CPU in per-mille format (0-1000). The
     *                  value 1000 means 100% (all cores).
	 * @throw std::invalid_argument If activity > 1000.
 	 */
	RequestCPU(technology_type_t tech_type, unsigned int clock_frequency, unsigned int nr_cores,
		   unsigned int activity)
	: Request(resource_type_t::CPU, tech_type), clock_frequency(clock_frequency), 
	  nr_cores(nr_cores), activity(activity)
	{
		if(activity > 1000) {
			throw std::invalid_argument("Activity invalid value (>1).");
		}
	}

	/**
	 * @brief The default virtual destructor (no dynamic memory used)
 	 */
	virtual ~RequestCPU() = default;

	/** @brief Getter for the clock frequency in MHz. */
	inline unsigned int get_clock_frequency() const {
		return this->clock_frequency;
	}

	/** @brief Setter for the clock frequency in MHz. */
	inline void set_clock_frequency(unsigned int clock_frequency) {
		this->clock_frequency = clock_frequency;
	}

	/** @brief Getter for the number of cores. */
	inline unsigned int get_nr_cores() const {
		return this->nr_cores;
	}

	/** @brief Setter for the number of cores. */
	inline void set_nr_cores(unsigned int nr_cores) {
		this->nr_cores = nr_cores;
	}

private:
	unsigned int clock_frequency;   ///< The clock frequency in MHz
	unsigned short nr_cores;	///< The number of cores
	unsigned short activity;	///< The per-mille value of current activity
};

using RequestGPU = RequestCPU; 		///< GPU currently has the same attributes of CPU

};	// libhwrel

#endif // LIBHWREL_H_
