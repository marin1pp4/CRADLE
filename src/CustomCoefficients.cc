#include "CRADLE/CustomCoefficients.hh"

#include <stdexcept>
#include <string>
#include <cmath>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "exprtk.hpp"
#include "CRADLE/Utilities.hh"
#include "CRADLE/Polarisation.hh"


namespace CRADLE {
namespace customcoefficients {

std::string GetCustomCoefficientExpression(
    const std::string& filename,
    const std::string& transition,
    const std::string& coefficient)
{
    boost::property_tree::ptree tree;
    boost::property_tree::read_json(filename, tree);
    const std::string path = transition + "." + coefficient;

    return tree.get<std::string>(path, "");
}


double EvaluateCustomCoefficient(
    const std::string& expressionString,
    const CustomCoefficientContext& ctx)
{
    double Ee = ctx.Ee;
    double rho = ctx.rho;
    double rM = ctx.rM;
    double mf = ctx.mf;
    double mgt = ctx.mgt;
    double ji = ctx.ji;
    double jf = ctx.jf;
    double Z = ctx.Z;
    double betaType = ctx.betaType;
    double xi = ctx.xi;
    double gamma = ctx.gamma;
    double gammaRatio = ctx.gammaRatio;
    double coulombCorr = ctx.coulombCorr;
    double bigLambda = ctx.bigLambda;
    double smallLambda = ctx.smallLambda;
    double CS_re = ctx.CS.real();
    double CSP_re = ctx.CSP.real();
    double CV_re = ctx.CV.real();
    double CVP_re = ctx.CVP.real();
    double CT_re = ctx.CT.real();
    double CTP_re = ctx.CTP.real();
    double CA_re = ctx.CA.real();
    double CAP_re = ctx.CAP.real();
    double CS_im = ctx.CS.imag();
    double CSP_im = ctx.CSP.imag();
    double CV_im = ctx.CV.imag();
    double CVP_im = ctx.CVP.imag();
    double CT_im = ctx.CT.imag();
    double CTP_im = ctx.CTP.imag();
    double CA_im = ctx.CA.imag();
    double CAP_im = ctx.CAP.imag();

    // --------------------------------------------------------
    // Create the ExprTk symbol table
    // --------------------------------------------------------
    exprtk::symbol_table<double> symbolTable;

    // Energy and nuclear parameters
    symbolTable.add_variable("Ee", Ee);
    symbolTable.add_variable("rho", rho);
    symbolTable.add_variable("rM", rM);
    symbolTable.add_variable("mf", mf);
    symbolTable.add_variable("mgt", mgt);
    symbolTable.add_variable("ji", ji);
    symbolTable.add_variable("jf", jf);
    symbolTable.add_variable("Z", Z);
    symbolTable.add_variable("betaType", betaType);
    symbolTable.add_variable("xi", xi);
    symbolTable.add_variable("gamma", gamma);
    symbolTable.add_variable("gammaRatio", gammaRatio);
    symbolTable.add_variable("coulombCorr", coulombCorr);
    symbolTable.add_variable("bigLambda", bigLambda);
    symbolTable.add_variable("smallLambda", smallLambda);

    // --------------------------------------------------------
    // Coupling constants: real components
    // --------------------------------------------------------

    // Short names are aliases for the real components.
    // This preserves the behavior of the original function.
    symbolTable.add_variable("CS", CS_re);
    symbolTable.add_variable("CSP", CSP_re);
    symbolTable.add_variable("CV", CV_re);
    symbolTable.add_variable("CVP", CVP_re);
    symbolTable.add_variable("CT", CT_re);
    symbolTable.add_variable("CTP", CTP_re);
    symbolTable.add_variable("CA", CA_re);
    symbolTable.add_variable("CAP", CAP_re);

    // Explicit real components
    symbolTable.add_variable("CS_re", CS_re);
    symbolTable.add_variable("CSP_re", CSP_re);
    symbolTable.add_variable("CV_re", CV_re);
    symbolTable.add_variable("CVP_re", CVP_re);
    symbolTable.add_variable("CT_re", CT_re);
    symbolTable.add_variable("CTP_re", CTP_re);
    symbolTable.add_variable("CA_re", CA_re);
    symbolTable.add_variable("CAP_re", CAP_re);

    // --------------------------------------------------------
    // Coupling constants: imaginary components
    // --------------------------------------------------------
    symbolTable.add_variable("CS_im", CS_im);
    symbolTable.add_variable("CSP_im", CSP_im);
    symbolTable.add_variable("CV_im", CV_im);
    symbolTable.add_variable("CVP_im", CVP_im);
    symbolTable.add_variable("CT_im", CT_im);
    symbolTable.add_variable("CTP_im", CTP_im);
    symbolTable.add_variable("CA_im", CA_im);
    symbolTable.add_variable("CAP_im", CAP_im);

    // --------------------------------------------------------
    // Compile the mathematical expression
    // --------------------------------------------------------
    exprtk::expression<double> expression;
    expression.register_symbol_table(symbolTable);
    exprtk::parser<double> parser;

    if (!parser.compile(expressionString, expression)) {

        throw std::runtime_error(
            "Could not compile custom coefficient expression: "
            + expressionString
        );
    }

    return expression.value();
}


// ============================================================
// Evaluate beta-decay correlation coefficients
// ============================================================

CorrelationCoefficients EvaluateCorrelationCoefficients(
    double Ee,
    const CustomCoefficientContext& ctx,
    const BetaTransitionParameters& parameters,
    const CustomCoefficientExpressions& expressions)
{
    // --------------------------------------------------------
    // 1. Prepare the energy-dependent context
    // --------------------------------------------------------
    CustomCoefficientContext localCtx = ctx;

    localCtx.Ee = Ee;
    const double m_e = utilities::EMASSC2;
    localCtx.gammaRatio = localCtx.gamma * m_e / Ee;
    localCtx.coulombCorr = utilities::FINESTRUCTURE * localCtx.Z / std::sqrt(Ee * Ee / (m_e * m_e) - 1.0);

    // --------------------------------------------------------
    // 2. Extract nuclear parameters and coupling constants
    // --------------------------------------------------------
    const double mf = localCtx.mf;
    const double mgt = localCtx.mgt;
    const double ji = localCtx.ji;
    const double jf = localCtx.jf;
    const double Z = localCtx.Z;
    const int betaType = static_cast<int>(localCtx.betaType);
    const auto& CS  = localCtx.CS;
    const auto& CSP = localCtx.CSP;
    const auto& CV  = localCtx.CV;
    const auto& CVP = localCtx.CVP;
    const auto& CT  = localCtx.CT;
    const auto& CTP = localCtx.CTP;
    const auto& CA  = localCtx.CA;
    const auto& CAP = localCtx.CAP;

    // --------------------------------------------------------
    // 3. Calculate standard correlation coefficients
    // --------------------------------------------------------
    CorrelationCoefficients coeff;

    if (parameters.calculationMode == CoefficientCalculationMode::LegacyReal) {

        coeff.a = utilities::CalculateBetaNeutrinoAsymmetry(
            CS.real(), CSP.real(), CT.real(),
            CTP.real(), CV.real(), CVP.real(),
            CA.real(), CAP.real(),
            mf, mgt,
            parameters.a_conf,
            parameters.b_conf
        );

        coeff.b = utilities::CalculateFierz(
            CS.real(), CSP.real(), CT.real(),
            CTP.real(), CV.real(), CVP.real(),
            CA.real(), CAP.real(),
            mf, mgt,
            parameters.a_conf,
            parameters.b_conf
        );

        } else {

            coeff.a = utilities::CalculateBetaNeutrinoAsymmetry(
                CS, CSP, CT, CTP,
                CV, CVP, CA, CAP,
                mf, mgt,
                parameters.a_conf,
                parameters.b_conf,
                Ee, Z, betaType
            );

            coeff.b = utilities::CalculateFierz(
                CS, CSP, CT, CTP,
                CV, CVP, CA, CAP,
                mf, mgt,
                parameters.a_conf,
                parameters.b_conf,
                Z, betaType
            );
        }

    // Spin-dependent coefficients
    if (ji > 0.0) {

        coeff.A = polarisation::CalculateBetaAssymetry(
            CS, CSP, CT, CTP,
            CV, CVP, CA, CAP,
            mf, mgt,
            ji, jf,
            betaType, Z, Ee
        );

        coeff.B = polarisation::CalculateNeutrinoAssymetry(
            CS, CSP, CT, CTP,
            CV, CVP, CA, CAP,
            mf, mgt,
            ji, jf,
            betaType, Z, Ee
        );

        coeff.D = polarisation::CalculateDTripleCorrelation(
            CS, CSP, CT, CTP,
            CV, CVP, CA, CAP,
            mf, mgt,
            ji, jf,
            betaType, Z, Ee
        );

        if (ji > 0.5) {

            coeff.c = polarisation::CalculateAlignmentCorrelation(
                CS, CSP, CT, CTP,
                CV, CVP, CA, CAP,
                mf, mgt,
                ji, jf,
                betaType, Z, Ee
            );
        }
    }

    // --------------------------------------------------------
    // 4. Apply custom coefficient expressions
    // --------------------------------------------------------
    if (!expressions.a.empty()) {
        coeff.a = EvaluateCustomCoefficient(expressions.a, localCtx);
    }

    if (!expressions.b.empty()) {
        coeff.b = EvaluateCustomCoefficient(expressions.b, localCtx);
    }

    if (ji > 0.0) {

        if (ji > 0.5 && !expressions.c.empty()) {
            coeff.c = EvaluateCustomCoefficient(expressions.c, localCtx);
        }

        if (!expressions.A.empty()) {
            coeff.A = EvaluateCustomCoefficient(expressions.A, localCtx);
        }

        if (!expressions.B.empty()) {
            coeff.B = EvaluateCustomCoefficient(expressions.B, localCtx);
        }

        // Preserve the existing CRADLE++ restriction
        if (ji == jf && !expressions.D.empty()) {
            coeff.D = EvaluateCustomCoefficient(expressions.D, localCtx);
        }
    }

    // --------------------------------------------------------
    // 5. Apply polarization and alignment factors
    // --------------------------------------------------------
    coeff.c *= -parameters.alignment;
    coeff.A *= parameters.polMag;
    coeff.B *= parameters.polMag;
    coeff.D *= parameters.polMag;

    return coeff;
}


} // namespace customcoefficients
} // namespace CRADLE
