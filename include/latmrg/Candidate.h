#ifndef LATMRG_CANDIDATE_H
#define LATMRG_CANDIDATE_H

#include <list>

#include "latticetester/FlexTypes.h"

namespace LatMRG {

template<typename Int>
struct Candidate
{
    Int modulus;
    IntVec aa;
    double merit = 0.0; // Also used when no FOM is calculated.

    const Int& getModulus() const { return modulus; }
    const IntVec& getAa() const { return aa; }
    double getMerit() const { return merit; }
};

template<typename Int>
class CandidateList
{
public:

    CandidateList(unsigned long maxLength, bool best)
        : m_max(maxLength), m_best(best)
    {
    }

    void add(const Candidate<Int>& candidate)
    {
        // List is ordered from best to worst.
        auto it = m_candidates.begin();

        for (; it != m_candidates.end(); ++it) {
            if (m_best ? candidate.merit > it->merit
                       : candidate.merit < it->merit)
                break;
        }

        if (m_candidates.size() >= m_max && it == m_candidates.end())
            return;

        m_candidates.insert(it, candidate);

        if (m_candidates.size() > m_max)
            m_candidates.pop_back();
    }

    void addWithoutMerit(const Candidate<Int>& candidate)
    {
        if (!full())
            m_candidates.push_back(candidate);
    }

    bool full() const
    {
        return m_candidates.size() >= m_max;
    }

    const std::list<Candidate<Int>>& getList() const
    {
        return m_candidates;
    }

    // Merit of the worst candidate currently retained.
    double getCutoffMerit() const
    {
        return m_candidates.back().merit;
    }

private:
    std::list<Candidate<Int>> m_candidates;
    unsigned long m_max;
    bool m_best;
};

} // namespace LatMRG

#endif