/**
 * @file
 * Header for the AirborneTransmissionProfile class.
 */

#pragma once

#include "pop/Person.h"
#include "util/RnHandler.h"

#include <boost/property_tree/ptree.hpp>
#include <vector>
#include <numeric>


namespace stride {

/**
 * AirborneTransmission probability from disease data.
 */
class AirborneTransmissionProfile
{
public:
	TransmissionProfile(): m_per_person_viral_shedding(1),
						   m_linking_hazard_virus(0.226),
						 {}

	/// Initialize.
	void Initialize(const boost::property_tree::ptree& configPT, const boost::property_tree::ptree& diseasePt);

	/// Return per-person viral shedding.
	double GetPerPersonViralShedding() const;

	/// Return coefficient linking hazard rate to virus load.
	double GetLinkingHazardVirus() const;

private:

    double            			m_per_person_viral_shedding; ///< Per-person viral shedding
    double             			m_linking_hazard_virus; ///< Coefficient linking harard rate to virus load

};

} // namespace stride

