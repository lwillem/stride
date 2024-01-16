
/**
 * @file
 * Implementation for the AirborneTransmissionProfile class.
 */

#include "AirborneTransmissionProfile.h"

#include "util/StringUtils.h"


namespace stride {

using namespace std;
using namespace boost::property_tree;
using namespace stride::util;

void AirborneTransmissionProfile::Initialize(const ptree& configPt, const ptree& diseasePt)
{
    // 1. setup general transmission aspects
    m_per_person_viral_shedding  = diseasePt.get<double>("disease.k", 1);
    m_linking_hazard_virus     = diseasePt.get<double>("disease.delta", 0.226);

}

double TransmissionProfile::GetPerPersonViralShedding() const {
	return m_per_person_viral_shedding;
}

double TransmissionProfile::GetLinkingHazardVirus() const {
	return m_linking_hazard_virus;
}

} // namespace stride
