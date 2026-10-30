
template<typename Int, typename Real>  void SearchMRG<Int, Real>::resetSearchState() {
        currentRandomRegion = 0;
        currentSubIndex = 0;
        hasLastAk = false;
        lastMaxPeriodResult = false;
        mrg.reset();
        currentAkTrial = 0;
        currentAkSample = 0;
        hasCurrentAk = false;
        WorkingLattice.reset();
        aa.SetLength(this->order + 1);
        alpha.SetLength(this->order + 1);
        if (lowBoundaries.length() == this->order + 1 && highBoundaries.length() == this->order + 1)  prepareEqualityGroups();
        exhaustiveIndex = Int(0);
    }
 /**
  * nextGeneratorExhaustive() performs an exhaustive enumeration of the admissible MRG 
  * generators defined by the configuration. The coefficients of the recurrence 
  * are obtained by decoding the current generator index into the corresponding 
  * multiplier values within the prescribed ranges. Each admissible set of 
  * coefficients is used to construct an MRG lattice that is subsequently analyzed 
  * by the search procedure defined in the base class.
  * If the search is restricted to generators with maximal period, then nextGeneratorExhaustive()
  * gives the next generator with this property (otherwise they are skipped).
  */
template<typename Int, typename Real>
bool SearchMRG<Int, Real>::nextGeneratorExhaustive()
{
    Int akStepSize(1), offsetK(0), tmp, range, offset;
    for (int i = 1; i < this->order; ++i)
        akStepSize *= highBoundaries[i] - lowBoundaries[i] + 1;
    // Enumerate vectors until the candidate limit or search region is exhausted.
    // The independent exhaustiveIndex always advances through the region. For each 
    // new a_k, test the period condition once. If it fails, skip the remaining vectors
    // with that a_k; otherwise return the next vector in aa.
    while (!this->reachedNumGeneratorLimit() && exhaustiveIndex < getNumMultipliers())
    {
        tmp = exhaustiveIndex;
        // Decode exhaustiveIndex as a mixed-radix number. For coefficient a_i,
        // the radix is the number of values in [lowBoundaries[i], highBoundaries[i]].
        // Taking the remainder selects its value; dividing by the radix prepares
        // tmp for the next coefficient. Thus each index identifies one vector.
        for (int i = 1; i <= this->order; ++i) {
            range = highBoundaries[i] - lowBoundaries[i] + 1;
            offset = tmp % range;
            tmp /= range;
            aa[i] = lowBoundaries[i] + offset;
            if (i == this->order)
                offsetK = offset;
        }
        ++exhaustiveIndex;
        ++this->numCandidatesExamined;
        // Test once when enumeration reaches a new value of a_k.
        if (this->maxPeriod) {
            if (!hasLastAk || aa[this->order] != lastAk) {
                lastMaxPeriodResult = getMRGComponent().maxPeriod(aa);
                lastAk = aa[this->order];
                hasLastAk = true;
            }
            if (!lastMaxPeriodResult) {
                // Start at the first vector with the next value of a_k.
                exhaustiveIndex = (offsetK + 1) * akStepSize;
                continue;
            }
            ++this->numMaxPeriodCandidates;
        }
        return true;
    }
    return false;
}
  //===========================================================================

  /**
  * nextGeneratorRandom() generates admissible MRG generators by randomly selecting 
  * the recurrence coefficients within the ranges specified in the configuration. 
  * It provides the random search analagoue of nextGeneratorExhaustive() and can be used by 
  * performSearch() when a stochastic search. 
  * If the search is restricted to generators with maximal period, then 
  * nextGeneratorRandom() only yields generators having this property.
  */
  template<typename Int, typename Real> bool SearchMRG<Int, Real>::nextGeneratorRandom()
  {    
    // Keep drawing coefficient vectors until one passes the requested period
    // test or the generator limit is reached. A failed vector is skipped;
    // a successful vector remains in aa when the function returns true.
    while (!this->reachedNumGeneratorLimit())
    {
        drawRandomCoefficients();
        ++this->numCandidatesExamined;
        if (this->maxPeriod) {
            if (!getMRGComponent().maxPeriod(aa)) continue;
            ++this->numMaxPeriodCandidates;
        }
        return true;
    }
    return false;
  }
  //===========================================================================

  /**
   * nextGeneratorRandomBox() performs a random-box search over the admissible
   * MRG coefficient region. It first selects random subregions, or boxes, within
   * the prescribed coefficient bounds. Then, each selected box is checked 
   * exhaustively. The side lengths of the boxes are specified by the subregion 
   * sizes h_i in the configuration.
   *
   * If the search is restricted to generators with maximal period, the maximal
   * period condition is checked before constructing the corresponding MRG
   * lattice. When possible, candidates sharing the same last coefficient a_k are
   * skipped together if this value cannot lead to maximal period, avoiding
   * unnecessary lattice constructions.
  */
  template<typename Int, typename Real> bool SearchMRG<Int, Real>::nextGeneratorRandomBox()
  {
    Int subregionTotal(1);
    Int akStepSize(1);    
    // Calculate subregion size
    for (int i = 1; i <= this->order; i++) {
       Int hi = subregionSize[i];
        subregionTotal *= hi;
        if (i < this->order) akStepSize *= hi;
    }
    // Visit the selected random boxes one at a time. Within each box,
    // currentSubIndex enumerates the coefficient vectors. When a box is
    // exhausted, draw the next one. Test the period condition once per
    // distinct a_k within the box; if it fails, skip the remaining vectors
    // with that a_k. Return true with the next accepted vector in aa, or
    // false when the search limits are reached.
    while (!this->reachedNumGeneratorLimit() && currentRandomRegion < numRandomRegions)
    {
        // If the current box is exhausted, move to the next random box.
        if (currentSubIndex >= subregionTotal) {
            ++currentRandomRegion;
            currentSubIndex = 0;
            hasLastAk = false;
            continue;
        }
        // Draw a new random box.
        if (currentSubIndex == 0) {
            for (int i = 1; i <= this->order; i++) {
                alpha[i] = randInt(lowBoundaries[i], highBoundaries[i] - subregionSize[i] + 1);
            }
            hasLastAk = false;
        }
        // Decode currentSubIndex into a vector inside the current box.
        Int tmp = currentSubIndex;
        Int offset;
        Int offsetK(0);
        for (int i = 1; i <= this->order; i++) {
            Int hi = subregionSize[i];
            offset = tmp % hi;
            tmp /= hi;
            aa[i] = alpha[i] + offset;
            // Save what is the current offset in the last variable
            if (i == this->order) offsetK = offset;
        }        
        ++this->numCandidatesExamined;
        // Check maximal period only when a_k has changed.
        if (this->maxPeriod) {
            bool akChanged = !hasLastAk || aa[this->order] != lastAk;
            if (akChanged) {
                setModulusIntP<Int>(this->modulus);
                lastMaxPeriodResult = getMRGComponent().maxPeriod(aa);
                lastAk = aa[this->order];
                hasLastAk = true;
            }
            // If this a_k fails, skip all remaining vectors with this a_k.
            if (!lastMaxPeriodResult) {
                // Skip the unexamined vectors with this a_k
                currentSubIndex = (offsetK + 1) * akStepSize;
                continue;
            }
            ++this->numMaxPeriodCandidates;
        }
        ++currentSubIndex;
        return true;
    }    
    return false;
}
//===========================================================================

/**
 * nextGeneratorRandomByAk() performs a random search organized around the
 * last multiplier a_k. It first draws values of a_k at random and
 * retains only those that satisfy the maximal-period requirement, if this
 * requirement is chosen by the user.
 *
 * For each accepted value of a_k, the function generates a prescribed number of
 * random coefficient vectors by drawing the remaining coefficients
 * a_1, ..., a_{k-1} independently within their configured bounds. 
 */
template<typename Int, typename Real> bool SearchMRG<Int, Real>::nextGeneratorRandomByAk()
{  
    // Continue until the generator limit or the limit on sampled a_k values
    // is reached. For each accepted a_k, generate up to
    // numRandomPointsPerAk vectors by drawing the other coefficients.
    // Return true with the next vector in aa, or false when a limit is reached.
    while (!this->reachedNumGeneratorLimit() && currentAkTrial < numRandomAk)
    {
        // Draw new values of a_k until one passes the period test or a search
        // limit is reached. Test each drawn a_k once. A failed value is skipped;
        // an accepted value is kept for the next numRandomPointsPerAk vectors.
        if (!hasCurrentAk || currentAkSample >= numRandomPointsPerAk) {
            hasCurrentAk = false;
            currentAkSample = 0;            
            while (currentAkTrial < numRandomAk && !this->reachedNumGeneratorLimit())
            {
                currentAk = randInt(lowBoundaries[this->order], highBoundaries[this->order]);
                aa[this->order] = currentAk;
                ++currentAkTrial;
                /*
                 * Check whether this a_k has maximal period.
                 */
                if (this->maxPeriod) {
                    if (!getMRGComponent().maxPeriod(aa)) {                        
                        ++this->numCandidatesExamined;
                        continue;
                    }
                }
                hasCurrentAk = true;
                break;
            }
            if (!hasCurrentAk)
                return false;
        }
        // Keeps the class of coefficients which equal a_k fixed and draws only the other coefficients
        drawRandomCoefficients(this->order);
        ++currentAkSample;
        ++this->numCandidatesExamined;   
        if (this->maxPeriod) ++this->numMaxPeriodCandidates;    
        return true;
    }
    return false;
}
//===========================================================================

/**
 * nextGeneratorRandomByAkBox() performs a hybrid random-box search organized
 * around the las mulitpiliers a_k. It first draws values of a_k at
 * random and retains only those that satisfy the maximal-period requirement, if
 * wanted by the user.
 *
 * For each accepted value of a_k, the function draws a random box for the
 * remaining coefficients a_1, ..., a_{k-1}. The box is defined by the configured
 * subregion sizes h_i for i = 1, ..., k-1, while a_k remains fixed, corresponding
 * to h_k = 1. The points in this box are then checked exhaustively.
 */
template<typename Int, typename Real> bool SearchMRG<Int, Real>::nextGeneratorRandomByAkBox()
{
    /*
     * h_k is assumed to be 1. Thus the box is only over the coordinates
     * a_1, ..., a_{k-1}.
     */
    Int boxSize(1);
    for (int i = 1; i < this->order; i++) {
        Int hi = subregionSize[i];
        boxSize *= hi;
    }
    // Continue while the generator limit has not been reached and either
    // a box remains to be examined or another a_k may be drawn. For each
    // accepted a_k, enumerate its box before moving to the next value.
    // Return true with the next vector in aa, or false when the search ends.
    while (!this->reachedNumGeneratorLimit() && (hasCurrentAk || currentAkTrial < numRandomAk))
    {
        /*
         * If there is no current good a_k, or if the current box has been
         * exhausted, draw a new a_k and a new random box for the other
         * coordinates.
         */
        if (!hasCurrentAk || currentSubIndex >= boxSize) {
            hasCurrentAk = false;
            currentSubIndex = 0;
            while (currentAkTrial < numRandomAk && !this->reachedNumGeneratorLimit())
            {
                /*
                 * Draw the lower corner of the box for coordinates
                 * a_1, ..., a_{k-1}.
                 */
                for (int i = 1; i < this->order; i++) {
                    alpha[i] = randInt(lowBoundaries[i], highBoundaries[i] - subregionSize[i] + 1);
                    aa[i] = alpha[i];
                }
                /*
                 * Draw a_k independently.
                 */
                currentAk = randInt(lowBoundaries[this->order], highBoundaries[this->order]);
                aa[this->order] = currentAk;
                ++currentAkTrial;
                /*
                 * Check whether this a_k is good. The check is performed
                 * only once, using the first point of the box.
                 */
                if (this->maxPeriod) {
                    if (!getMRGComponent().maxPeriod(aa)) {                        
                        ++this->numCandidatesExamined;
                        continue;
                    }
                }
                hasCurrentAk = true;
                break;
            }
            if (!hasCurrentAk)
                return false;
        }
        /*
         * Enumerate the current box exhaustively for a_1, ..., a_{k-1},
         * while keeping a_k fixed. The same logic as for the exhaustive search applies.
         */
        Int tmp = currentSubIndex;
        for (int i = 1; i < this->order; i++) {
            Int hi = subregionSize[i];
            Int offset = tmp % hi;
            tmp /= hi;
            aa[i] = alpha[i] + offset;
        }
        aa[this->order] = currentAk;
        ++currentSubIndex;
        ++this->numCandidatesExamined;
        if (this->maxPeriod) ++this->numMaxPeriodCandidates;
        return true;
    }
    return false;
}
//===========================================================================

/**
 * Generates MRG candidates whose coefficients satisfy a sufficient condition
 * for approximate factoring.
 *
 * For each coefficient a_i, the function draws an integer j uniformly from
 * 0 < j < sqrt(m). It then randomly chooses |a_i| to be either j or
 * floor(m / j), and independently chooses the sign. The usual coefficient
 * boundaries and equality constraints are not used.
 */
template<typename Int, typename Real> bool SearchMRG<Int, Real>::nextGeneratorRandomApproxFactoring()
{
    Int j, maxJ, draw;
    long bits;
    NTL::SqrRoot(maxJ, this->modulus);
    if (maxJ * maxJ == this->modulus) --maxJ;
    // Repeatedly draw a complete coefficient vector in approximate-factoring
    // form until one passes the requested period test or the generator limit
    // is reached. Return true with the accepted vector in aa, or false when
    // no further candidates may be drawn.
    while (!this->reachedNumGeneratorLimit()) {
        for (int i = 1; i <= this->order; ++i) {
            // A single uniform draw encodes both j and two choice bits.
            // Dividing by 4 gives j in [1, maxJ]; the remainder selects
            // either j or floor(modulus / j) and determines its sign.
            draw = randInt(Int(0), 4 * maxJ - 1);
            j = draw / 4 + 1;
            bits = NTL::conv<long>(draw % 4);
            aa[i] = (bits & 1) ? this->modulus / j : j;
            if (bits & 2) aa[i] = -aa[i];
        }
        ++this->numCandidatesExamined;
        if (this->maxPeriod) {
            if (!getMRGComponent().maxPeriod(aa)) continue;
            ++this->numMaxPeriodCandidates;
        }
        return true;
    }
    return false;
}
//===========================================================================

/**
 * Like nextGeneratorRandomApproxFactoring() but keeps an
 * accepted a_k fixed for several random points.
 */
template<typename Int, typename Real> bool SearchMRG<Int, Real>::nextGeneratorRandomApproxFactoringByAk()
{
    Int j, maxJ, draw;
    long bits;
    NTL::SqrRoot(maxJ, this->modulus);
    if (maxJ * maxJ == this->modulus) --maxJ;
    // Continue while the generator limit permits more candidates and either
    // an accepted a_k has points remaining or another a_k may be drawn.
    // For each accepted a_k, draw fresh approximate-factoring values for
    // the other coefficients and return the resulting vector in aa.
    while (!this->reachedNumGeneratorLimit() && (hasCurrentAk || currentAkTrial < numRandomAk))
    {
        if (!hasCurrentAk ||
            currentAkSample >= numRandomPointsPerAk) {
            hasCurrentAk = false;
            currentAkSample = 0;
            // Draw a_k in approximate-factoring form until one passes the period
            // check or a search limit is reached. 
            while (currentAkTrial < numRandomAk && !this->reachedNumGeneratorLimit())
            {
                draw = randInt(Int(0), 4 * maxJ - 1);
                j = draw / 4 + 1;
                bits = NTL::conv<long>(draw % 4);
                currentAk = (bits & 1) ? this->modulus / j : j;
                if (bits & 2) currentAk = -currentAk;
                aa[this->order] = currentAk;
                ++currentAkTrial;
                if (this->maxPeriod &&
                    !getMRGComponent().maxPeriod(aa)) {
                    ++this->numCandidatesExamined;
                    continue;
                }
                hasCurrentAk = true;
                break;
            }
            if (!hasCurrentAk)
                return false;
        }
        // Draw fresh coefficients while keeping the max period a_k fixed.
        for (int i = 1; i < this->order; ++i) {
            draw = randInt(Int(0), 4 * maxJ - 1);
            j = draw / 4 + 1;
            bits = NTL::conv<long>(draw % 4);
            aa[i] = (bits & 1) ? this->modulus / j : j;
            if (bits & 2) aa[i] = -aa[i];
        }
        aa[this->order] = currentAk;
        ++currentAkSample;
        ++this->numCandidatesExamined;
        if (this->maxPeriod) ++this->numMaxPeriodCandidates;
        return true;
    }
    return false;
}
//===========================================================================

/**
 * Generates random MRG coefficients as signed sums of powers of two.
 *
 * For each a_i, chooses between 1 and n_i distinct exponents from
 * 0, ..., p_i, then independently chooses the sign of each term.
 * Equality constraints and coefficient boundaries are not used.
 */
template<typename Int, typename Real> bool SearchMRG<Int, Real>::nextGeneratorRandomPowersOfTwo()
{
    Int power, value;
    long ni, pi, choice, maxTerms, terms;
    // Repeatedly draw a complete coefficient vector whose entries are signed
    // sums of powers of two. Skip vectors that fail the requested period test
    // and return true with the next accepted vector in aa. Return false when
    // the generator limit is reached. 
    while (!this->reachedNumGeneratorLimit()) {
        for (int i = 1; i <= this->order; ++i) {
            // These std::vectors are indexed from 0: entry i-1 describes a_i.
            ni = numPowerTwo[i - 1];
            pi = maxPowerTwo[i - 1];
            maxTerms = std::min(ni, pi + 1);
            terms = randInt(1,maxTerms);
            // Draw distinct exponents without replacement.
            std::vector<long> exponents(pi + 1);
            for (long e = 0; e <= pi; ++e)
                exponents[e] = e;
            value = 0;
            for (long j = 0; j < terms; ++j) {
                choice = randInt(j, pi);
                std::swap(exponents[j], exponents[choice]);
                power = NTL::power(NTL::ZZ(2), exponents[j]);
                if (randInt(0,1) == 0) value += power;
                else value -= power;
            }
            aa[i] = value;
        }
        ++this->numCandidatesExamined;
        if (this->maxPeriod) {
            if (!getMRGComponent().maxPeriod(aa)) continue; // Tested and counted; draw another vector.
            ++this->numMaxPeriodCandidates;
        }
        return true;
    }
    return false;
}
//===========================================================================

/**
 * Like nextGeneratorRandomPowerofTwo but for each accepted a_k 
 * it draws several independent vectors for the remaining
 * coefficients while keeping a_k fixed.
 */
template<typename Int, typename Real> bool SearchMRG<Int, Real>::nextGeneratorRandomPowersOfTwoByAk()
{
    // Continue while the generator limit allows more candidates and either
    // an accepted a_k has points remaining or another a_k may be drawn.
    // For each accepted a_k, draw fresh signed-power-of-two values for the
    // other coefficients and return the resulting vector in aa.
    while (!this->reachedNumGeneratorLimit() && (hasCurrentAk || currentAkTrial < numRandomAk))
    {
        // Find a new acceptable a_k when needed.
        if (!hasCurrentAk ||
            currentAkSample >= numRandomPointsPerAk) {
            hasCurrentAk = false;
            currentAkSample = 0;
            // Draw a_k as a signed sum of powers of two until one passes the period
            // check or a search limit is reached.
            while (currentAkTrial < numRandomAk && !this->reachedNumGeneratorLimit())
            {
                currentAk = drawPowerOfTwoCoefficient(this->order);
                aa[this->order] = currentAk;
                ++currentAkTrial;
                if (this->maxPeriod && !getMRGComponent().maxPeriod(aa)) {
                    ++this->numCandidatesExamined;
                    continue;
                }
                hasCurrentAk = true;
                break;
            }
            if (!hasCurrentAk)
                return false;
        }
        // Draw a fresh vector for the accepted a_k.
        for (int i = 1; i < this->order; ++i)
            aa[i] = drawPowerOfTwoCoefficient(i);
        aa[this->order] = currentAk;
        ++currentAkSample;
        ++this->numCandidatesExamined;
        if (this->maxPeriod)
            ++this->numMaxPeriodCandidates;
        return true;
    }

    return false;
}
//===========================================================================
