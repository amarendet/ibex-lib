//============================================================================
//                                  I B E X
// File        : ibex_NodeType.h
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : 18 Aug 2014
//============================================================================

#ifndef __IBEX_NODE_TYPE_H__
#define __IBEX_NODE_TYPE_H__

namespace ibex {

/** \ingroup iset */
/*@{*/

/**
 * \brief Status of a node in the representation of an i-set.
 *
 * __IBEX_IN__       : the node box is inside the set
 * __IBEX_OUT__      : the node box is outside the set
 * __IBEX_UNK__      : the node is of unknown status (typically, corresponds to a boundary box)
 * __IBEX_UNK_IN     : the subtree contains IN and UNK leaves
 * __IBEX_UNK_OUT    : the subtree contains OUT and UNK leaves
 * __IBEX_UNK_IN_OUT : the subtree contains IN, OUT and UNK leaves. Used also for a temporary leaf (too large to be considered as boundary).
 */
typedef enum { __IBEX_IN__,
	           __IBEX_OUT__,
	           __IBEX_UNK__,
	           __IBEX_UNK_IN__,
	           __IBEX_UNK_OUT__,
	           __IBEX_UNK_IN_OUT__ } NodeType;

/**
 * \brief Status of a union of two nodes **in the sync sense**.
 *
 * \warning The sync sense means that __IBEX_IN__ | __IBEX_OUT__ gives __IBEX_UNK_IN_OUT__ and
 * not __IBEX_IN__ !
 *
 * Ex: __IBEX_IN__ | __IBEX_UNK__ gives __IBEX_IN_UNK__.
 *
 * \see SetBisect constructor.
 */
NodeType operator|(NodeType x, NodeType y);

/**
 * \brief Status of an intersection of two nodes **in the sync sense**
 *
 * \warning The sync sense means that __IBEX_IN__ & __IBEX_OUT__ is impossible (throw NoSet) and
 * not __IBEX_OUT__ !
 *
 * Ex: __IBEX_IN__ & __IBEX_UNK__ gives __IBEX_IN__.
 *
 * \see SetBisect constructor.
 */
NodeType operator&(NodeType x, NodeType y);

/**
 * \brief Status of an intersection of two nodes
 *
 *
 * Ex: __IBEX_IN__ & __IBEX_UNK__ gives __IBEX_UNK__.
 *     __IBEX_UNK_OUT__ & __IBEX_UNK__ gives __IBEX_UNK_OUT__.
 *
 * \see SetBisect constructor.
 */
NodeType inter(NodeType x, NodeType y);

/**
 * \brief Status of an union of two nodes
 *
 *
 * Ex: __IBEX_IN__ & __IBEX_OUT__ gives __IBEX_UNK_IN_OUT__.

 *
 * \see SetBisect constructor.
 */
NodeType _union(NodeType x, NodeType y);

/**
 * \brief Status of the subset of a node
 *
 * Ex: subset(__IBEX_UNK_IN__) = __IBEX_UNK__
 *
 * \see SetBisect constructor.
 */
NodeType subset(NodeType x);


/**
 * \brief False only if the subtree contains no inner point
 */
bool possibly_contains_in(NodeType x);

/**
 * \brief False only if the subtree contains no outer point
 */
bool possibly_contains_out(NodeType x);

/**
 * \brief True only if the subtree contains inner points
 */
bool certainly_contains_in(NodeType x);

/**
 * \brief True only if the subtree contains outer points
 */
bool certainly_contains_out(NodeType x);

/**
 * \brief Convert the status to char
 */
char to_string(const NodeType& status);

/*@}*/

} // namespace ibex

#endif // __IBEX_NODE_TYPE_H__
