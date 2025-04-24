#pragma once

#ifndef BOLTZMANNSOLVER__SPEC_LIB_NAME__
  #define BOLTZMANNSOLVER__SPEC_LIB_NAME__

  #ifndef __cplusplus
    #define __cplusplus
  #endif

  #include "avgsvcalculator.hpp"
  #include "config.hpp"
  #include "relicparam.hpp"

namespace __SPEC_LIB_NAME__
{

  class BoltzmannSolver : public AvgSvCalculator, public Relicparam_t
  {
public:
    short unsigned int solver{1}; //!< Switch for linear or logarithmic differential equation solver

private:
    real_t Tfo{0.}; //!< Freeze-out temperature
    real_t T_prev_boltzphi;
    real_t rhophi_prev_boltzphi;
    real_t T_prev_boltzphi_log;
    real_t rhophi_prev_boltzphi_log;

public:
    /**
     * @brief This method resets the private attributes uniques in the instance of BoltzmannSolver
     */
    void reset_prevs();

    /**
      @brief Sets the model for the QCD equation of state.
      @param x Integer representing the QCD equation of state model.
      @see Documentation of Relicparam_t.
    */
    void setQCDeosModel(const int x)
    {
      Relicparam_t::init(getMassLBSM(), x);
      reset_prevs();
    };

    /**
     * @brief Constructor with given Param_t and QCD model.
     *
     * @param param Param_t structure with input values.
     * @param x Integer corresponding to a QCD EOS model. Default x=2.
     */
    BoltzmannSolver(const Param_t& param, int x = 2) : AvgSvCalculator(param), Relicparam_t(getMassLBSM(), x)
    {
      reset_prevs();
    };

    /**
     * @brief Constructor with given AvgSvCalculator and QCD model.
     *
     * @param sop AvgSvCalculator class with input values.
     * @param x Integer corresponding to a QCD EOS model. Default x=2.
     */
    BoltzmannSolver(const AvgSvCalculator& sop, int x = 2) : AvgSvCalculator(sop), Relicparam_t(getMassLBSM(), x)
    {
      reset_prevs();
    };

    /**
     * @brief Constructor with given input, process list and QCD model
     *
     * @param param Param_t structure with input values.
     * @param proclist Process list to consider for the calculation of \f$\langle \sigma v \rangle \f$
     * @param x Integer corresponding to a QCD EOS model. Default x=2.
     */
    BoltzmannSolver(const Param_t& param, std::shared_ptr<std::vector<Process2to2>> proclist, int x = 2)
            : AvgSvCalculator(param, proclist), Relicparam_t(getMassLBSM(), x)
    {
      reset_prevs();
    }

    /**
     * @brief Method to change the input parameters of a BoltzmannSovler instance
     *
     * @param other The new structure of the input parameters
     * @return unsigned int 0 for success
     */
    unsigned int changeInput(const Param_t& other)
    {
      AvgSvCalculator::changeInput(other);
      Relicparam_t::init(other.getLightestBSMmass(), model_eff);
      reset_prevs();
      return 0;
    };

    /**
     * @brief Method to print a BoltzmannSolver instance. Overrides previously defined virtual methods.
     *
     * Prints:
     *  1. The process list
     *  2. The table of \f$\tilde W(\sqrt s)\f$, defined as \f$\tilde W(\sqrt s) = g^2_{DM} W_{eff}(\sqrt s)\f$
     *  3. The Relicparam_t elements
     * @param out The output stream. Default out=std::cout.
     */
    void print(std::ostream& out = std::cout) const override
    {
      print_procs(out);
      print_g2_wefftable(out);
      print_relicparam(out);
      out << "solver = " << solver << '\n';
      out << "Tfo = " << Tfo << '\n';
    };

    /**
     * @brief Calculates the \f$Y_{eq}(T)\f$ at a given temperature \f$T\f$.
     *
     * @param T The input temperature.
     * @return real_t The value of \f$Y_{eq}(T)\f$.
     */
    real_t Yeq(const real_t& T);

    real_t Sigmatildestar_init(const real_t& T, const real_t& delta);

    /**
     * @brief  Calculates the \f$\frac{dY_{eq}(T)}{dT}\f$ at a given temperature \f$T\f$.
     *
     * @param T The input temperature.
     * @return real_t The value of  \f$\frac{dY_{eq}(T)}{dT}\f$.
     */
    real_t dYeq_dT(const real_t& T);

    /**
     * @brief Computes the right-hand side of the Boltzmann equation with Y as unknown.
     *
     * If the computation is made without quintessence or decaying scalar field, Yphi=0.
     *
     * @param T The temperature
     * @param sv The value of \f$\langle \sigma v \rangle \f$
     * @param Y The value of Y
     * @param Yphi The value of \f$Y_\phi\f$
     * @return real_t the right-hand side of the Boltzmann equation with Y as unknown
     */
    real_t boltzright(const real_t& T, const real_t& sv, const real_t& Y, const real_t& Yphi);

    /**
     * @brief Computes the freeze-out temperature, and updates the corresponding data member.
     *
     * @param delta The \f$\delta\f$ parameter to solve the implicit equation. Default is 1.5.
     */
    void setTfo(real_t delta = 1.5);

    /**
     * @brief Returns the freeze-out temperature.
     *
     * @return real_t
     */
    real_t getTfo() const { return Tfo; };

    real_t boltzright_phi(const real_t& T, const real_t& Y, const real_t& Yphi);

    /**
     * @brief Computes the right-hand side of the logarithmic Boltzmann equation with Y as unknown.
     *
     * If the computation is made without quintessence or decaying scalar field, Yphi=0.
     *
     * @param T The temperature
     * @param sv The value of \f$\langle \sigma v \rangle \f$
     * @param Y The value of Y
     * @param Yphi The value of \f$Y_\phi\f$
     * @return real_t the right-hand side of the logarithmic Boltzmann equation with Y as unknown
     */
    real_t boltzright_log(const real_t& T, const real_t& sv, const real_t& Y, const real_t& Yphi);

    /**
     * @brief Computes the right-hand side of the logarithmic Boltzmann equation with Y as unknown,
     *        with model with \f$φ\f$
     *
     * If the computation is made without quintessence or decaying scalar field, Yphi=0.
     *
     * @param T The temperature
     * @param Y The value of Y
     * @param Yphi The value of \f$Y_\phi\f$
     * @return real_t the right-hand side of the logarithmic Boltzmann equation with Y as unknown
     */
    real_t boltzright_phi_log(const real_t& T, const real_t& Y, const real_t& Yphi);

    /**
     * @brief Solves the linear Boltzmann equation, and returns the value for the relic density \f$Ω_{DM}h^2\f$
     *
     * @return real_t \f$Ω_{DM}h^2\f$
     */
    real_t relic_density_linsolver();

    /**
     * @brief Solves the logarithmic Boltzmann equation, and returns the value for the relic density \f$Ω_{DM}h^2\f$
     *
     * @return real_t \f$Ω_{DM}h^2\f$
     */
    real_t relic_density_logsolver();

    /**
     * @brief Calls an algorithm to solve the Boltzmann equation, and returns the value for the relic density
     * \f$Ω_{DM}h^2\f$
     *
     * The algoritm's choice depends on the data memeber solver:
     *  - If solver=1 (default) uses the logarithmic solver
     *  - If solver=2 uses the linear solver
     *
     * @return real_t \f$Ω_{DM}h^2\f$
     */
    real_t relic_density();
  };

} // End of namespace __SPEC_LIB_NAME__

inline std::ostream& operator<<(std::ostream& out, const __SPEC_LIB_NAME__::BoltzmannSolver& var)
{
  var.print(out);
  return out;
}

#endif // BOLTZMANNSOLVER__SPEC_LIB_NAME__