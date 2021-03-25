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
 *  Copyright 2020, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Header for the Health class.
 */

#pragma once

namespace stride {

/// Enumerate the various health states with respect to the infection.
enum class HealthStatus : unsigned short int
{
        Susceptible              = 0U,
        Exposed                  = 1U,
        Infectious               = 2U,
        Symptomatic              = 3U,
        InfectiousAndSymptomatic = 4U,
        Recovered                = 5U,
};

/// Holds a person's health data.
class Health
{
public:
        ///
        explicit Health(unsigned short int start_infectiousness = 0U, unsigned int short start_symptomatic = 0U,
                        unsigned short int time_infectious = 0U, unsigned short int time_symptomatic = 0U,
						double sympt_cnt_reduction_work_school = 0U, double sympt_cnt_reduction_community=0U,
						double relative_susceptibility = 0U);

        ///
        unsigned short int GetEndInfectiousness() const { return m_end_infectiousness; }

        ///
        unsigned short int GetEndSymptomatic() const { return m_end_symptomatic; }

        ///
        unsigned short int GetStartInfectiousness() const { return m_start_infectiousness; }

        ///
        unsigned short int GetStartSymptomatic() const { return m_start_symptomatic; }

        ///
        unsigned int GetIdIndexCase() const { return m_id_index_case; }

        ///
        unsigned int GetIdInfector() const { return m_id_infector; }

        /// Is this person infected?
        bool IsInfected() const
        {
                return m_status == HealthStatus::Exposed || m_status == HealthStatus::Infectious ||
                       m_status == HealthStatus::InfectiousAndSymptomatic || m_status == HealthStatus::Symptomatic;
        }

        /// Is this person infectious.
        bool IsInfectious() const
        {
                return m_status == HealthStatus::Infectious || m_status == HealthStatus::InfectiousAndSymptomatic;
        }

        /// Is this person recovered?
        bool IsRecovered() const { return m_status == HealthStatus::Recovered; }

        /// Is this person susceptible?
        bool IsSusceptible() const { return m_status == HealthStatus::Susceptible; }

        /// Is this person symptomatic?
        bool IsSymptomatic() const
        {
                return m_status == HealthStatus::Symptomatic || m_status == HealthStatus::InfectiousAndSymptomatic;
        }

        /// Is this person exposed?
		bool IsExposed() const
		{
				return m_status == HealthStatus::Exposed;
		}

        /// Have the symptoms started today?
        bool SymptomsStartedToday() const { return GetDiseaseCounter() == m_start_symptomatic; }

        /// Have the symptoms started X days before?
        bool NumberDaysSymptomatic(unsigned int days_before) const;

        /// Is infected X days before?
        bool NumberDaysInfected(unsigned int days_before) const;

        /// Set health state to susceptible
        void SetSusceptible() { m_status = HealthStatus::Susceptible; }

        /// Start the infection.
        void StartInfection(unsigned int id_index_case, unsigned int id_infector, 
                double relative_infectiousness);

        /// Stop the infection.
        void StopInfection();

        /// Update progress of the disease.
        void Update();

        /// Get contact reduction in school/work pools when symptomatic infected
        double GetSymptomaticCntReductionWorkSchool() const { return m_sympt_cnt_reduction_work_school; };

        /// Get contact reduction in community pools when symptomatic infected
        double GetSymptomaticCntReductionCommunity() const { return m_sympt_cnt_reduction_community; };

		/// Is this individual PCR detectable?
		bool IsPcrDetectable(unsigned int detectable_delay) const { return GetDiseaseCounter() >= detectable_delay; }

        /// Get relative infectiousness
        double GetRelativeInfectiousness() const {
        		return m_relative_infectiousness;
        }

        /// Get relative probability to acquire infection
		double GetRelativeSusceptibility() const {
				return m_relative_susceptibility;
		}


private:
        /// Get the disease counter.
        unsigned short int GetDiseaseCounter() const { return m_disease_counter; }

        /// Increment disease counter.
        void IncrementDiseaseCounter() { m_disease_counter++; }

        /// Reset the disease counter.
        void ResetDiseaseCounter() { m_disease_counter = 0U; }

private:
        unsigned short int m_disease_counter; ///< The disease counter.
        HealthStatus       m_status;          ///< The current status of the person w.r.t. the disease.

        unsigned short int m_start_infectiousness; ///< Days after infection to become infectious.
        unsigned short int m_start_symptomatic;    ///< Days after infection to become symptomatic.
        unsigned short int m_end_infectiousness;   ///< Days after infection to end infectious state.
        unsigned short int m_end_symptomatic;      ///< Days after infection to end symptomatic state.

        unsigned int       m_id_index_case;        ///< ID of the index case, given infection
        unsigned int       m_id_infector;          ///< ID of the infector, given infection

        double             m_sympt_cnt_reduction_work_school;  ///< Proportional reduction of presence in work/school pool when symptomatic
        double             m_sympt_cnt_reduction_community;    ///< Proportional reduction of presence in the community pools when symptomatic

        double			   m_relative_infectiousness;   ///< Relative probability of transmission when infected [0-1]
        double			   m_relative_susceptibility;   ///< Relative probability of acquiring infection upon exposure [0-1]

};

} // namespace stride
