/* ============================================================================
 * I B E X - Symbolic diff tests
 * ============================================================================
 * Copyright   : Ecole des Mines de Nantes (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Gilles Chabert
 * Created     : Mar 14, 2013
 * ---------------------------------------------------------------------------- */

#include "TestExprDiff.h"
#include "ibex_ExprDiff.h"

using namespace std;

namespace ibex {

void TestExprDiff::linear01() {

	Variable x,y;
	Function f(x,y,2*x+3*y);
	Function df(f,Function::DIFF);
	const ExprConstant* c=dynamic_cast<const ExprConstant*>(&df.expr());
	TEST_ASSERT(c);
	TEST_ASSERT(c->dim.type()==Dim::ROW_VECTOR);
	double _grad[][2] = {{2,2},{3,3}};
	IntervalVector grad(2,_grad);
	TEST_ASSERT(c->get_vector_value()==grad);
}

void TestExprDiff::poly01() {

	Variable x("x"),y("y");
	Function f(x,y,((sqr(x)+2*x*y)+pow(y,3))+1);
	Function df(f,Function::DIFF);
	const ExprVector* v=dynamic_cast<const ExprVector*>(&df.expr());
	TEST_ASSERT(v);
	TEST_ASSERT(v->dim.type()==Dim::ROW_VECTOR);
	TEST_ASSERT(sameExpr(v->arg(0),"((2*x)+(2*y))")
			|| sameExpr(v->arg(0),"((2*y)+(2*x))"));
	TEST_ASSERT(sameExpr(v->arg(1),"((2*x)+(3*y^2))"));
}

void TestExprDiff::one_var_one_func() {
	Variable x("x");
	Function f(x,sqr(x)+1);
	Function df(f,Function::DIFF);
	TEST_ASSERT(df.expr().dim.type()==Dim::SCALAR);
	TEST_ASSERT(sameExpr(df.expr(),"(2*x)"));
}

void TestExprDiff::vec01() {

	Variable x(4,"x");
	Function f(x,x[0]+x[3]);
	Function df(f,Function::DIFF);
	const ExprConstant* c=dynamic_cast<const ExprConstant*>(&df.expr());
	TEST_ASSERT(c);
	TEST_ASSERT(c->dim.type()==Dim::ROW_VECTOR);
	double _grad[][2] = {{1,1},{0,0},{0,0},{1,1}};
	IntervalVector grad(4,_grad);

	TEST_ASSERT(c->get_vector_value()==grad);
}

void TestExprDiff::vec02() {
	Variable x(4,"x");
	Function f(x,x[0]*x[3]);
	Function df(f,Function::DIFF);
	const ExprVector* v=dynamic_cast<const ExprVector*>(&df.expr());
	TEST_ASSERT(v);
	TEST_ASSERT(v->dim.type()==Dim::ROW_VECTOR);
	TEST_ASSERT(sameExpr(v->arg(0),x[3]));
	TEST_ASSERT(sameExpr(v->arg(3),x[0]));
}

void TestExprDiff::vec03() {
	Variable x(2,"x");
	Array<const ExprNode> _vec1(x[0],x[1]);
	const ExprNode& vec1=ExprVector::new_(_vec1,false);
	Function f(x,vec1[1]);
	Function df(f,Function::DIFF);
	const ExprConstant* c=dynamic_cast<const ExprConstant*>(&df.expr());
	TEST_ASSERT(c);
	TEST_ASSERT(c->dim.type()==Dim::ROW_VECTOR);
	double _v[][2]= {{0,0},{1,1}};
	TEST_ASSERT(c->get_vector_value()==IntervalVector(2,_v));
}

void TestExprDiff::mat01() {
	Variable x(2,2,"x");
	Function f(x,x[1][0]);
	Function df(f,Function::DIFF);
	const ExprConstant* c=dynamic_cast<const ExprConstant*>(&df.expr());
	TEST_ASSERT(c!=NULL);
	//TEST_ASSERT(c->dim.type()==Dim::MATRIX);
	double _M[][2]= {{0,0},{0,0},{1,1},{0,0}};
	//TEST_ASSERT(c->get_matrix_value()==IntervalMatrix(2,2,_M));
	TEST_ASSERT(c->get_vector_value()==IntervalVector(4,_M));
}

void TestExprDiff::mat02() {
	Variable x(2,2,"x");
	Array<const ExprNode> vec1(x[0],x[1]);
	Function f(x,vec1[1][0]);
	Function df(f,Function::DIFF);
	const ExprConstant* c=dynamic_cast<const ExprConstant*>(&df.expr());
	TEST_ASSERT(c!=NULL);
	//TEST_ASSERT(c->dim.type()==Dim::MATRIX);
	double _M[][2]= {{0,0},{0,0},{1,1},{0,0}};
	//TEST_ASSERT(c->get_matrix_value()==IntervalMatrix(2,2,_M));
	TEST_ASSERT(c->get_vector_value()==IntervalVector(4,_M));
}

void TestExprDiff::apply01() {
	Variable x("x");
	Function f(x,sqr(x),"f");
	Function g(x,f(3*x));
	Function dg(g,Function::DIFF);
	TEST_ASSERT(sameExpr(f.diff().expr(),"(2*x)"));
	TEST_ASSERT(sameExpr(dg.expr(),"(3*df((3*x)))"));
}

void TestExprDiff::apply02() {
	Variable x("x");
	Function f(x,sqr(x),"f");
	Function g(x,3*f(x));
	Function dg(g,Function::DIFF);
	TEST_ASSERT(sameExpr(dg.expr(),"(df(x)*3)"));
}

void TestExprDiff::apply03() {
	Variable x("x"),y("y");
	Function f(x,y,x*y,"f");
	Function g(x,y,f(2*x,3*y));
	Function dg(g,Function::DIFF);
	TEST_ASSERT(sameExpr(dg.expr(),"((2*df((2*x),(3*y))[0]),(3*df((2*x),(3*y))[1]))"));
	double _box[][2]={{1,1},{2,2}};
	double _dg_box[][2]={{12,12},{6,6}};
	IntervalVector dg_box(dg.eval_vector(IntervalVector(2,_box)));
	TEST_ASSERT(dg_box==IntervalVector(2,_dg_box));
}

} // end namespace
