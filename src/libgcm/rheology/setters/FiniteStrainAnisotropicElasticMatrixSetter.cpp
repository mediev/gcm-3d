#include "libgcm/rheology/setters/FiniteStrainAnisotropicElasticMatrixSetter.hpp"

using namespace gcm;

unsigned int FiniteStrainAnisotropicElasticMatrixSetter::getNumberOfStates() const
{
    return 0;
}

unsigned int FiniteStrainAnisotropicElasticMatrixSetter::getStateForNode(const ICalcNode& node) const
{
    return 0;
}

void FiniteStrainAnisotropicElasticMatrixSetter::setX(gcm_matrix& a, const MaterialPtr& material, const ICalcNode& node)
{
    a.clear();

    auto p = material->getRheologyProperties();
    auto rho = node.getRho();

    a(0,3) = a(1,4) = a(2,5) = -1.0/rho;

	a(3,0) = -p.c11;	a(3,1) = -p.c16-node.sxy;					a(3,2) = -p.c15-node.sxz;
	a(4,0) = -p.c16;	a(4,1) = -p.c66-(node.syy-node.sxx)/2.0;	a(4,2) = -p.c56-node.syz/2.0;
	a(5,0) = -p.c15;	a(5,1) = -p.c56-node.syz/2.0;				a(5,2) = -p.c55-(node.szz-node.sxx)/2.0;
	a(6,0) = -p.c12;	a(6,1) = -p.c26+node.sxy;					a(6,2) = -p.c25;
	a(7,0) = -p.c14;	a(7,1) = -p.c46+node.sxz/2.0;				a(7,2) = -p.c45+node.sxy/2.0;
	a(8,0) = -p.c13;	a(8,1) = -p.c36;							a(8,2) = -p.c35+node.sxz;
	
	a(9,0) = rho;
}

void FiniteStrainAnisotropicElasticMatrixSetter::setY(gcm_matrix& a, const MaterialPtr& material, const ICalcNode& node)
{
    a.clear();

    auto p = material->getRheologyProperties();
    auto rho = node.getRho();

    a(0,4) = a(1,6) = a(2,7) = -1.0/rho;

	a(3,0) = -p.c16+node.sxy;					a(3,1) = -p.c12;	a(3,2) = -p.c14;
	a(4,0) = -p.c66-(node.sxx-node.syy)/2.0;	a(4,1) = -p.c26;	a(4,2) = -p.c46-node.sxz/2.0;
	a(5,0) = -p.c56+node.syz/2.0;				a(5,1) = -p.c25;	a(5,2) = -p.c45+node.sxy/2.0;
	a(6,0) = -p.c26-node.sxy;					a(6,1) = -p.c22;	a(6,2) = -p.c24-node.syz;
	a(7,0) = -p.c46-node.sxz/2.0;				a(7,1) = -p.c24;	a(7,2) = -p.c44-(node.szz-node.syy)/2.0;
	a(8,0) = -p.c36;							a(8,1) = -p.c23;	a(8,2) = -p.c34+node.syz;
	
	a(9,1) = rho;
}

void FiniteStrainAnisotropicElasticMatrixSetter::setZ(gcm_matrix& a, const MaterialPtr& material, const ICalcNode& node)
{
    a.clear();

    auto p = material->getRheologyProperties();
    auto rho = node.getRho();

    a(0,5) = a(1,7) = a(2,8) = -1.0/rho;

	a(3,0) = -p.c15+node.sxz;					a(3,1) = -p.c14;							a(3,2) = -p.c13;
	a(4,0) = -p.c56+node.syz/2.0;				a(4,1) = -p.c46+node.sxz/2.0;				a(4,2) = -p.c36;
	a(5,0) = -p.c55-(node.sxx-node.szz)/2.0;	a(5,1) = -p.c45-node.sxy/2.0;				a(5,2) = -p.c35;
	a(6,0) = -p.c25;							a(6,1) = -p.c24+node.syz;					a(6,2) = -p.c23;
	a(7,0) = -p.c45-node.sxy/2.0;				a(7,1) = -p.c44-(node.syy-node.szz)/2.0;	a(7,2) = -p.c34;
	a(8,0) = -p.c35-node.sxz;					a(8,1) = -p.c34-node.syz;					a(8,2) = -p.c33;
	
	a(9,2) = rho;
}
