/* ============================================================================
 * I B E X - ibex_utils.h
 * ============================================================================
 * Copyright   : IMT Atlantique (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Antoine Marendet, Gilles Chabert
 * Created     : May 4, 2018
 * ---------------------------------------------------------------------------- */
 
#ifndef __SIP_IBEX_UTILS_H__
#define __SIP_IBEX_UTILS_H__

#include "ibex_Function.h"
#include "ibex_IntervalVector.h"
#include "ibex_Vector.h"

#include <vector>
#include <set>
#include <map>

namespace ibex {
Interval centeredFormEval(const Function& function, const IntervalVector& arg);
std::vector<IntervalVector> bisectAllDim(const IntervalVector& iv);

bool isfinite(const Vector& v);

IntervalVector sip_to_ext_box(const IntervalVector& box, const Interval& g);
IntervalVector sip_from_ext_box(const IntervalVector& ext_box);
Interval sip_goal_from_ext_box(const IntervalVector& ext_box);

Vector sip_to_ext_box(const Vector& box, double g);
Vector sip_from_ext_box(const Vector& ext_box);
double sip_goal_from_ext_box(const Vector& ext_box);


std::string print_mma(const Vector& iv);
std::string print_mma(const IntervalVector& iv);
std::string print_mma(const std::vector<IntervalVector>& path);
std::string print_mma_path(const std::vector<IntervalVector>& path);

int symbol_array_dim(const Array<const ExprSymbol>& array);

template<class T, class C, class A>
bool set_contains(const std::set<T, C, A>& s, const T& value) {
	// return find(s.begin(), s.end(), value) != s.end();
    return s.find(value) != s.end();
}

template<class T, class C, class A>
bool multiset_contains(const std::multiset<T, C, A>& s, const T& value) {
	return find(s.begin(), s.end(), value) != s.end();
}

template<class K, class V>
bool map_contains_key(const std::map<K, V>& m, const K& key) {
    return m.find(key) != m.end();
}

template<class K, class V>
V map_get_with_default(const std::map<K, V>& m, const K& key, V def) {
    auto it = m.find(key);
	if(it == m.cend()) {
		return def;
	} else {
		return it->second;
	}
}

} // end namespace ibex

#endif // __SIP_IBEX_UTILS_H__
