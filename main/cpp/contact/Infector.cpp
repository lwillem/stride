/*
 *  This is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *  The software is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License
 *  along with the software. If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright 2020, Willem L, Kuylen E, Broeckhove J
 */

/**
 * @file
 * Implementation of Infector algorithms.
 */

#include "Infector.h"

#include "ContactPool.h"
#include "calendar/Calendar.h"
#include "pop/Person.h"

using namespace std;

namespace {

/// Primary LOG_POLICY policy, implements LogMode::None.
/// \tparam LL
template <EventLogMode::Id LL>
class LOG_POLICY
{
public:
        static void Contact(const std::shared_ptr<spdlog::logger>&, const Person*, const Person*, ContactType::Id,
                            unsigned short int, const double, const double)
        {
        }

        static void Trans(const std::shared_ptr<spdlog::logger>&, const Person*, const Person*, ContactType::Id,
                          unsigned short int, unsigned int)
        {
        }
};

/// Specialized LOG_POLICY policy LogMode::Incidence.
template <>
class LOG_POLICY<EventLogMode::Id::Incidence>
{
public:
        static void Contact(const std::shared_ptr<spdlog::logger>&, const Person*, const Person*, ContactType::Id,
                            unsigned short int, const double, const double)
        {
        }

        // p1: infector & p2:infectee
        static void Trans(const std::shared_ptr<spdlog::logger>& logger, const Person* p1, const Person* p2,
                          ContactType::Id type, unsigned short int sim_day, unsigned int id_index_case)
        {
                logger->info("[TRAN_M] {} {} {} {} {}",
							 p2->GetAge(),
							 sim_day,
							 p2->GetHealth().GetStartInfectiousness(),
							 p2->GetHealth().GetStartSymptomatic(),
							 p2->GetHealth().GetEndSymptomatic()
							 );
        }
};

/// Specialized LOG_POLICY policy LogMode::Transmissions.
template <>
class LOG_POLICY<EventLogMode::Id::Transmissions>
{
public:
        static void Contact(const std::shared_ptr<spdlog::logger>&, const Person*, const Person*, ContactType::Id,
                            unsigned short int, const double, const double)
        {
        }

        // p1: infector & p2:infectee
        static void Trans(const std::shared_ptr<spdlog::logger>& logger, const Person* p1, const Person* p2,
                          ContactType::Id type, unsigned short int sim_day, unsigned int id_index_case)
        {
                logger->info("[TRAN] {} {} {} {} {} {} {} {} {} {} {} {} {} {}", p2->GetId(), p1->GetId(), p2->GetAge(), p1->GetAge(),
                             ToString(type), sim_day, id_index_case,
							 p2->GetHealth().GetStartInfectiousness(),p2->GetHealth().GetEndInfectiousness(),
							 p2->GetHealth().GetStartSymptomatic(),p2->GetHealth().GetEndSymptomatic(),
							 p1->GetHealth().IsSymptomatic(),
							 p2->GetHealth().GetRelativeInfectiousness(),
							 p2->GetHealth().GetRelativeSusceptibility());
        }
};

/// Specialized LOG_POLICY policy LogMode::All.
template <>
class LOG_POLICY<EventLogMode::Id::All>
{
public:
        static void Contact(const std::shared_ptr<spdlog::logger>& logger, const Person* p1, const Person* p2,
                            ContactType::Id type, unsigned short int sim_day, const double cProb, const double tProb)
        {
                if (p1->IsSurveyParticipant()) {
                        logger->info("[CONT] {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}", p1->GetId(), p1->GetAge(),
                                     p2->GetAge(), static_cast<unsigned int>(type == ContactType::Id::Household),
                                     static_cast<unsigned int>(type == ContactType::Id::K12School),
                                     static_cast<unsigned int>(type == ContactType::Id::College),
                                     static_cast<unsigned int>(type == ContactType::Id::Workplace),
                                     static_cast<unsigned int>(type == ContactType::Id::PrimaryCommunity),
                                     static_cast<unsigned int>(type == ContactType::Id::SecondaryCommunity),
									 static_cast<unsigned int>(type == ContactType::Id::HouseholdCluster),
									 static_cast<unsigned int>(type == ContactType::Id::Collectivity),
									 sim_day,
									 cProb, tProb,p2->GetHealth().IsSymptomatic(),p1->GetHealth().IsSymptomatic());
                }
        }

        static void Trans(const std::shared_ptr<spdlog::logger>& logger, const Person* p1, const Person* p2,
                          ContactType::Id type, unsigned short int sim_day, unsigned int id_index_case)
        {
                logger->info("[TRAN] {} {} {} {} {} {} {} {} {} {} {} {} {} {}", p2->GetId(), p1->GetId(), p2->GetAge(), p1->GetAge(),
                             ToString(type), sim_day, id_index_case,
							 p2->GetHealth().GetStartInfectiousness(),p2->GetHealth().GetEndInfectiousness(),
							 p2->GetHealth().GetStartSymptomatic(),p2->GetHealth().GetEndSymptomatic(),
							 p1->GetHealth().IsSymptomatic(),  p1->GetHealth().GetRelativeInfectiousness(),
							 p2->GetHealth().GetRelativeSusceptibility());
        }
};

} // namespace

namespace {

using namespace stride;
using namespace stride::ContactType;
using namespace stride::util;

inline double GetContactProbability(const AgeContactProfile& profile, const Person* p1, const Person* p2,
		size_t pool_size, const ContactType::Id pType,
		double cnt_reduction_school, double cnt_reduction_collectivity, double cnt_reduction_intergeneration, unsigned int cnt_reduction_intergeneration_cutoff,
		std::shared_ptr<Population>& population, double cnt_intensity_householdCluster,std::shared_ptr<Calendar> calendar)
{

		// initiate a contact adjustment factor, to account for physical distancing and/or contact intensity
    	double cnt_adjustment_factor = 1;

    		// Check if one of the persons is a non-complier to social distancing measures in this particular pooltype
    		if ((not p1->IsNonComplier(pType)) and (not p2->IsNonComplier(pType))) {

    			// account for physical distancing at work
    			if(pType == Id::Workplace){

    				double workplace_distancing_factor = calendar->GetWorkplaceDistancingFactor();
    				cnt_adjustment_factor = (1-workplace_distancing_factor);
    			}

    			// account for physical distancing in the community
    			if((pType == Id::PrimaryCommunity || pType == Id::SecondaryCommunity)){

    				double community_distancing_factor = calendar->GetCommunityDistancingFactor();

    				// apply inter-generation distancing factor if age cutoff is > 0 and at least one age is > cutoff
    				if((cnt_reduction_intergeneration > 0) &&
    					((p1->GetAge() > cnt_reduction_intergeneration_cutoff) || (p2->GetAge() > cnt_reduction_intergeneration_cutoff))){
    					cnt_adjustment_factor = (1-cnt_reduction_intergeneration);
    				} else {
    					// apply uniform community distancing
    					cnt_adjustment_factor = (1-community_distancing_factor);
    				}
    			}
    			// account for physical distancing at school
    			if(pType == Id::K12School){

    				// get distancing at school
    				double school_distancing_factor = cnt_reduction_school ;

//    				// work in progress //TODO
//    				double school_distancing_factor = calendar->GetSchoolDistancingFactor(p1->GetAge()) ;
//    				double school_distancing_factor2 = calendar->GetSchoolDistancingFactor(p2->GetAge()) ;
//					if(school_distancing_factor < school_distancing_factor2){
//						school_distancing_factor = school_distancing_factor2;
//					}
//
//					if(cnt_reduction_school != school_distancing_factor){
//						cout << cnt_reduction_school << " - "<< school_distancing_factor  << " - "<< school_distancing_factor2<< endl;
//
//					}

    				cnt_adjustment_factor = (1-school_distancing_factor);
    			}
    			// account for physical distancing in the collectivity
				if(pType == Id::Collectivity){
					cnt_adjustment_factor = (1-cnt_reduction_collectivity);
				}

    		}


		// account for contact intensity in household clusters
		if(pType == Id::HouseholdCluster){
			cnt_adjustment_factor = cnt_intensity_householdCluster;
		}


		// get the reference number of contacts, given age and age-contact profile
		double reference_num_contacts_p1{profile[EffectiveAge(static_cast<unsigned int>(p1->GetAge()))]};
        double reference_num_contacts_p2{profile[EffectiveAge(static_cast<unsigned int>(p2->GetAge()))]};
        const double potential_num_contacts{static_cast<double>(pool_size - 1)};

        // adjust contact for distancing
        reference_num_contacts_p1 *= cnt_adjustment_factor;
        reference_num_contacts_p2 *= cnt_adjustment_factor;


        // special case: reduce the number of community contacts if part of a HouseholdCluster
        if(cnt_intensity_householdCluster > 0 && (pType == Id::PrimaryCommunity || pType == Id::SecondaryCommunity)){

        	// get the number of non-household members in the HouseholdCluster
        	double householdCluster_non_household_members_p1 = static_cast<double>(population->GetPoolSize(Id::HouseholdCluster,p1)) -
        											                  population->GetPoolSize(Id::Household,p1);

        	double householdCluster_non_household_members_p2 = static_cast<double>(population->GetPoolSize(Id::HouseholdCluster,p2)) -
        	        											      population->GetPoolSize(Id::Household,p2);

        	// account for negative values ==>> not part of a HouseholdCluster
        	householdCluster_non_household_members_p1 = householdCluster_non_household_members_p1 < 0 ? 0 : householdCluster_non_household_members_p1;
        	householdCluster_non_household_members_p2 = householdCluster_non_household_members_p2 < 0 ? 0 : householdCluster_non_household_members_p2;

        	// account for the HouseholdCluster size and contact intensity
        	reference_num_contacts_p1 -= householdCluster_non_household_members_p1 * cnt_intensity_householdCluster;
        	reference_num_contacts_p2 -= householdCluster_non_household_members_p2 * cnt_intensity_householdCluster;

        	// account for negative values
        	reference_num_contacts_p1 = reference_num_contacts_p1 < 0 ? 0 : reference_num_contacts_p1;
        	reference_num_contacts_p2 = reference_num_contacts_p2 < 0 ? 0 : reference_num_contacts_p2;
        }

        // calculate the contact probability based on the (possible) number of contacts
        double individual_contact_probability_p1 = reference_num_contacts_p1 / potential_num_contacts;
        double individual_contact_probability_p2 = reference_num_contacts_p2 / potential_num_contacts;

        // use the minimum of both age-specific probabilities
        double contact_probability = individual_contact_probability_p1 ;
		if(individual_contact_probability_p2 < individual_contact_probability_p1){
			contact_probability = individual_contact_probability_p2;
		}

	    // limit probability to 0.999
        if (contact_probability >= 1) {
        	contact_probability = 0.999;
        }


    	// assume fully connected households
		if(pType == Id::Household){
			contact_probability = 0.999;
		}

		// exclude contacts with household members within household cluster
		if(pType == Id::HouseholdCluster){
			contact_probability = (p1->GetPoolId(Id::Household) == p2->GetPoolId(Id::Household)) ? 0.0 : cnt_intensity_householdCluster;
		}


        return contact_probability;
}

} // namespace

namespace stride {

using namespace stride::util;

//-------------------------------------------------------------------------------------------------
// Definition for ContactLogMode::Contacts,
// both with track_index_case false and true.
//-------------------------------------------------------------------------------------------------
template <EventLogMode::Id LL, bool TIC, bool TO>
void Infector<LL, TIC, TO>::Exec(ContactPool& pool, const AgeContactProfile& profile,
                                 const TransmissionProfile& transProfile, util::RnHandler& rnHandler,
                                 unsigned short int simDay, shared_ptr<spdlog::logger> eventLogger,
								 double cnt_reduction_school,
								 double cnt_reduction_collectivity,
								 double cnt_reduction_intergeneration, unsigned int cnt_reduction_intergeneration_cutoff,
								 std::shared_ptr<Population> population, double m_cnt_intensity_householdCluster,
								 std::shared_ptr<Calendar> calendar)
{
        using LP = LOG_POLICY<LL>;

        // set up some stuff
        const auto  pType    = pool.m_pool_type;
        const auto& pMembers = pool.m_members;
        const auto  pSize    = pMembers.size();

        // check all contacts
        for (size_t i_person1 = 0; i_person1 < pSize; i_person1++) {
                // check if member is present today
                const auto p1 = pMembers[i_person1];
                if (!p1->IsInPool(pType)) {
                        continue;
                }
                // loop over possible contacts (contacts can be initiated by each member)
                for (size_t i_person2 = i_person1 + 1; i_person2 < pSize; i_person2++) {
                        // check if not the same person
                        if (i_person1 == i_person2) {
                                continue;
                        }
                        // check if member is present today
                        const auto p2 = pMembers[i_person2];
                        if (!p2->IsInPool(pType)) {
                                continue;
                        }
                        // check for contact
                        const double cProb = GetContactProbability(profile, p1, p2, pSize, pType,
                        		cnt_reduction_school, cnt_reduction_collectivity, cnt_reduction_intergeneration,
								cnt_reduction_intergeneration_cutoff,population,m_cnt_intensity_householdCluster,calendar);
                        if (rnHandler.Binomial(cProb)) {
								const auto  tProb_p1_p2    = transProfile.GetProbability(p1,p2);
								const auto  tProb_p2_p1    = transProfile.GetProbability(p2,p1);

                                // log contact if person 1 is participating in survey
                                LP::Contact(eventLogger, p1, p2, pType, simDay, cProb, tProb_p1_p2);
                                // log contact if person 2 is participating in survey
                                LP::Contact(eventLogger, p2, p1, pType, simDay, cProb, tProb_p2_p1);


                                // if track&trace is in place, option to register (both) contact(s)
                                p1->RegisterContact(p2);
                                p2->RegisterContact(p1);

                                // transmission & infection.
                                // note: no tertiary infections with TIC; so mark new case as 'recovered'
                                auto& h1 = p1->GetHealth();
                                auto& h2 = p2->GetHealth();

								// if h1 infectious, account for susceptibility of p2
								if (h1.IsInfectious() && h2.IsSusceptible() &&
									rnHandler.Binomial(tProb_p1_p2)) {

										double rel_inf = transProfile.GetIndividualInfectiousness(rnHandler);
										h2.StartInfection(h1.GetIdIndexCase(),p1->GetId(),rel_inf);

										if (TIC)
												h2.StopInfection();
										LP::Trans(eventLogger, p1, p2, pType, simDay, h1.GetIdIndexCase());
								}

								// if h2 infectious, account for susceptibility of p1
								if (h2.IsInfectious() && h1.IsSusceptible() &&
									rnHandler.Binomial(tProb_p2_p1)) {

									double rel_inf = transProfile.GetIndividualInfectiousness(rnHandler);
									h1.StartInfection(h2.GetIdIndexCase(),p2->GetId(), rel_inf);

										if (TIC)
												h1.StopInfection();
										LP::Trans(eventLogger, p2, p1, pType, simDay, h2.GetIdIndexCase());
								}
                        }
                }
        }
}

//-------------------------------------------------------------------------------------------
// Definition for ContactLogMode::None and ContactLogMode::Transmission
// both with track_index_case false and true.
//-------------------------------------------------------------------------------------------
template <EventLogMode::Id LL, bool TIC>
void Infector<LL, TIC, true>::Exec(ContactPool& pool, const AgeContactProfile& profile,
                                   const TransmissionProfile& transProfile, util::RnHandler& rnHandler,
                                   unsigned short int simDay, shared_ptr<spdlog::logger> eventLogger,
								   double cnt_reduction_school,
								   double cnt_reduction_collectivity,
								   double cnt_reduction_intergeneration, unsigned int cnt_reduction_intergeneration_cutoff,
								   std::shared_ptr<Population> population, double m_cnt_intensity_householdCluster,
								   std::shared_ptr<Calendar> calendar)
{
        using LP = LOG_POLICY<LL>;

        // check for infected members and sort
        bool   infectious_cases;
        size_t num_cases;
        tie(infectious_cases, num_cases) = pool.SortMembers();

        if (!infectious_cases) {
                return;
        }

        // set up some stuff
        const auto  pType    = pool.m_pool_type;
        const auto  pImmune  = pool.m_index_immune;
        const auto& pMembers = pool.m_members;
        const auto  pSize    = pMembers.size();

        // match infectious and susceptible members, skip last part (immune members)
        for (size_t i_infected = 0; i_infected < num_cases; i_infected++) {
                // check if member is present today
                const auto p1 = pMembers[i_infected];
                if (!p1->IsInPool(pType)) {
                        continue;
                }
                auto& h1 = p1->GetHealth();
                if (h1.IsInfectious()) {
                        // loop over possible susceptible contacts
                        for (size_t i_contact = num_cases; i_contact < pImmune; i_contact++) {
                                // check if member is present today
                                const auto p2 = pMembers[i_contact];
                                if (!p2->IsInPool(pType)) {
                                        continue;
                                }
                                const double cProb_p1 = GetContactProbability(profile, p1, p2, pSize, pType,
                                								cnt_reduction_school,
																cnt_reduction_collectivity, cnt_reduction_intergeneration, cnt_reduction_intergeneration_cutoff,
															population, m_cnt_intensity_householdCluster,calendar);
                                const auto  tProb_p1_p2   = transProfile.GetProbability(p1,p2);
                                if (rnHandler.Binomial(cProb_p1, tProb_p1_p2)) {

                                        auto& h2 = p2->GetHealth();
                                        if (h1.IsInfectious() && h2.IsSusceptible()) {
                                                double rel_inf = transProfile.GetIndividualInfectiousness(rnHandler);
                                                h2.StartInfection(h1.GetIdIndexCase(),p1->GetId(), rel_inf);

                                                // if track&trace is in place, option to register (both) contact(s)
                                                p1->RegisterContact(p2); //TODO: make use of "log policy" template

                                                // No secondary infections with TIC; just mark p2 'recovered'
                                                if (TIC)
                                                        h2.StopInfection();
                                                LP::Trans(eventLogger, p1, p2, pType, simDay, h1.GetIdIndexCase());
                                        }
                                }
                        }
                }
        }
}

//--------------------------------------------------------------------------
// All explicit instantiations.
//--------------------------------------------------------------------------
template class Infector<EventLogMode::Id::None, false>;
template class Infector<EventLogMode::Id::None, true>;
template class Infector<EventLogMode::Id::Incidence, false>;
template class Infector<EventLogMode::Id::Incidence, true>;
template class Infector<EventLogMode::Id::Transmissions, false>;
template class Infector<EventLogMode::Id::Transmissions, true>;
template class Infector<EventLogMode::Id::All, false>;
template class Infector<EventLogMode::Id::All, true>;


} // namespace stride
