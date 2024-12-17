#pragma once
#ifndef __cplusplus
  #define __cplusplus
#endif

#include "config.hpp"
#include <array>
#include <vector>

namespace __SPEC_LIB_NAME__
{

  constexpr const static int NTABMAX =
      1000; //!< This constexpr variable contains the maximum size for the Relicparam_t::table_rhoPD data structure
  constexpr const static unsigned int size_data_temperature =
      276; //!< This constexpr variable contains the exact size for the Relicparam_t::dataT data structure

  /** This type contains a temperature T and the corresponding values for \f$\sqrt{g_{eff}^\ast}\f$ and \f$h_{eff}\f$
   * at that temperature.
   */
  struct DataWithTemperature_t
  {
    real_t T;         //!< Temperature
    real_t sqrtgstar; //!< Square root of the effective number of relativistic degrees of freedom
    real_t heff;      //!< Effective entropy degrees of freedom

    /**
      @brief Construct an empty DataWithTemperature_t object
    */
    DataWithTemperature_t(){};

    /**
      @brief Construct a DataWithTemperature_t object with given values
      @param T_in Temperature
      @param sqrtgstar_in Square root of the effective number of relativistic degrees of freedom
      @param heff_in Effective entropy degrees of freedom
    */
    DataWithTemperature_t(const real_t& T_in, const real_t& sqrtgstar_in, const real_t& heff_in)
            : T(T_in), sqrtgstar(sqrtgstar_in), heff(heff_in){};

    /**
      @brief Print the contents of the DataWithTemperature_t object
      @param out Output stream to print to (default is std::cout)
    */
    void print(std::ostream& out = std::cout) const
    {
      out << "T = " << T << ' ';
      out << "sqrtgstar = " << sqrtgstar << ' ';
      out << "heff = " << heff << '\n';
    };
  };

  struct Relicparam_t //< structure containing the cosmological model parameters
  {

public:
    /** @brief Enumeration for representing the QCD equation of state model
     */
    enum QCDeosModel : int
    {
      A = 1,   /**<   1: Model A (ignores hadrons)*/
      B,       /**<   2: Model B (Tc = 154 MeV, free meson and hadron gas)*/
      B2,      /**<   3: Model B2 (variation of B, scaled by 0.9)*/
      B3,      /**<   4: Model B3 (variation of B, scaled by 1.1)*/
      C,       /**<   5: Model C (Tc = 185.5 MeV, free meson and hadron gas)*/
      Bonn,    /**<*/
      IdealGas /**<   Other: Old Model (ideal gas)*/
    };

protected:
    csl::InitSanitizer<int> model_eff{"model_eff"}; //!< Variable storing which QCD EOS model is being used

private:
    std::array<DataWithTemperature_t, size_data_temperature> dataT;
    /**< @brief This variable contains \f$h_{eff}\f$ and \f$\sqrt{g_{eff}^\ast}\f$
     *   as functions of the temperature.
     *   Data are stored with decreasing temperature.
     */

public:
    //  Original members from SuperIso Relic V. 4

    csl::InitSanitizer<int> entropy_model{"entropy_model"}; //!< Entropy model switch
    csl::InitSanitizer<int> energy_model{"energy_model"};   //!< Energy model switch

    //  Ref. is Eq. (33) of [Manual SuperIso Relic 3.1]
    csl::InitSanitizer<real_t> dd0{"dd0"};     //!< Dark energy density (rho_D) at BBN
    csl::InitSanitizer<real_t> ndd{"ndd"};     //!< Dark energy density decreasing exponent
    csl::InitSanitizer<real_t> Tdend{"Tdend"}; //!< Dark energy density cutoff temperature
    csl::InitSanitizer<real_t> Tddeq{"Tddeq"}; //!< Temperature at which rho_D = rho_rad (radiation energy density)

    csl::InitSanitizer<real_t> sd0{"sd0"};     //!< Dark entropy density at BBN
    csl::InitSanitizer<real_t> nsd{"nsd"};     //!< Dark entropy density decreasing exponent
    csl::InitSanitizer<real_t> Tsend{"Tsend"}; //!< Dark entropy density cutoff temperature

    csl::InitSanitizer<real_t> Sigmad0{"Sigmad0"};       //!< Dark entropy injection at BBN
    csl::InitSanitizer<real_t> nSigmad{"nSigmad"};       //!< Dark entropy injection slope
    csl::InitSanitizer<real_t> TSigmadend{"TSigmadend"}; //!< Dark entropy injection cutoff temperature

    csl::InitSanitizer<real_t> Sigmarad0{"Sigmarad0"};       //!< Standard entropy injection amplitude
    csl::InitSanitizer<real_t> nSigmarad{"nSigmarad"};       //!< Standard entropy injection slope
    csl::InitSanitizer<real_t> TSigmaradend{"TSigmaradend"}; //!< Standard entropy injection cutoff temperature

    csl::InitSanitizer<real_t> nt0{"nt0"};     //!< Non-thermal production of relics amplitude
    csl::InitSanitizer<real_t> nnt{"nnt"};     //!< Non-thermal production of relics slope
    csl::InitSanitizer<real_t> Tnend{"Tnend"}; //!< Non-thermal production of relics cutoff temperature

    csl::InitSanitizer<int> coupd{"coupd"}; //!< Dark fluid coupling to plasma switch

    csl::InitSanitizer<real_t> quintn2{"quintn2"};   //!< Effective quintessence model n2 parameter
    csl::InitSanitizer<real_t> quintn3{"quintn3"};   //!< Effective quintessence model n3 parameter
    csl::InitSanitizer<real_t> quintn4{"quintn4"};   //!< Effective quintessence model n4 parameter
    csl::InitSanitizer<real_t> quintT12{"quintT12"}; //!< Effective quintessence model T12 parameter
    csl::InitSanitizer<real_t> quintT23{"quintT23"}; //!< Effective quintessence model T23 parameter
    csl::InitSanitizer<real_t> quintT34{"quintT34"}; //!< Effective quintessence model T34 parameter

    csl::InitSanitizer<int> phi_model{"phi_model"};            //!< Decaying scalar field model switch
    csl::InitSanitizer<real_t> eta_phi{"eta_phi"};             //!< Decaying scalar field model eta parameter
    csl::InitSanitizer<real_t> Gamma_phi{"Gamma_phi"};         //!< Decaying scalar field model decay rate
    csl::InitSanitizer<real_t> rhot_phi_Tmax{"rhot_phi_Tmax"}; //!< Decaying scalar field model energy density at Tmax
    csl::InitSanitizer<real_t> n_phi{"n_phi"};                 //!< Decaying scalar field model power law index
    csl::InitSanitizer<real_t> rhot_phi0{"rhot_phi0"};         //!< Decaying scalar field model initial energy density
    csl::InitSanitizer<real_t> Tphi0{"Tphi0"};                 //!< Decaying scalar field model initial temperature


    csl::InitSanitizer<real_t> T_RH{"T_RH"};                             //!< Temperature at the reheating epoch
    csl::InitSanitizer<real_t> Sigmatildestar{"Sigmatildestar"};         //!<
    csl::InitSanitizer<real_t> Sigmatildestar_max{"Sigmatildestar_max"}; //!<
    csl::InitSanitizer<real_t> Tstdstar_max{"Tstdstar_max"};             //!<

    csl::InitSanitizer<real_t> mgravitino{"mgravitino"}; //!< Gravitino mass

    csl::InitSanitizer<real_t> relicmass{"relicmass"}; //!< Mass of the relic particles
    csl::InitSanitizer<int> scalar{"scalar"};          //!< Switch for the MSSM annihilation processes

    csl::InitSanitizer<int> solver{"solver"}; //!< Switch for linear or logarithmic differential equation solver
    csl::InitSanitizer<int> beta_samples{"beta_samples"}; //!< Number of samples for the Coulomb log

    csl::InitSanitizer<real_t> Tfo{"Tfo"};   //!< Freeze-out temperature
    csl::InitSanitizer<real_t> Tmax{"Tmax"}; //!< Maximum temperature

    csl::InitSanitizer<int> full_comput{
        "full_comput"}; //!< Switch to deactivate the fast freeze-out temperature determination

    csl::InitSanitizer<int> use_table_rhoPD{"use_table_rhoPD"}; //!< Switch for using the table of rho * PD
    csl::InitSanitizer<real_t> table_rhoPD[2][NTABMAX];
    /**< @brief Table of (Temperature T, dark energy density rho_PD)
         If this table is defined, it overrides the standard parametrisations for
         dark energy density that can be defined with the prevoius methods
    */

    csl::InitSanitizer<size_t> size_table_rhoPD{"size_table_rhoPD"}; //!< Size of the table of rho * PD

    /*---------------------*/
    /* AlterBBN parameters */
    /*---------------------*/

    csl::InitSanitizer<int> err{"err"};                                  //!< Error code
    csl::InitSanitizer<int> failsafe{"failsafe"};                        //!< Switch for the integration method
    csl::InitSanitizer<real_t> eta0{"eta0"};                             //!< Initial baryon to photon ratio
    csl::InitSanitizer<real_t> Nnu{"Nnu"};                               //!< Number of neutrinos (e+- included)
    csl::InitSanitizer<real_t> dNnu{"dNnu"};                             //!< Number of extra neutrinos (delta N_nu)
    csl::InitSanitizer<real_t> life_neutron{"life_neutron"};             //!< Neutron lifetime
    csl::InitSanitizer<real_t> life_neutron_error{"life_neutron_error"}; //!< Error on the neutron lifetime
    csl::InitSanitizer<real_t> xinu1{"xinu1"};                           //!< Degeneracy parameter for e- and neutrino
    csl::InitSanitizer<real_t> xinu2{"xinu2"};                           //!< Degeneracy parameter for muon and neutrino
    csl::InitSanitizer<real_t> xinu3{"xinu3"};                           //!< Degeneracy parameter for tau and neutrino
    csl::InitSanitizer<real_t> m_chi{"m_chi"};                           //!< Mass of WIMP
    csl::InitSanitizer<real_t> g_chi{"g_chi"};                           //!< Coupling of WIMP to SM particles
    csl::InitSanitizer<real_t> Tinit{"Tinit"};                           //!< Initial temperature
    csl::InitSanitizer<real_t> Tnudec{"Tnudec"};                         //!< Neutrino decoupling temperature
    csl::InitSanitizer<int> wimp{"wimp"};                                //!< Switch to enable (1) / disable (0) WIMPs
    csl::InitSanitizer<int> SMC_wimp{"SMC_wimp"};
    /**< @brief WIMP coupling to SM particles.
     *    - 1 for EM,
     *    - 2 for neutrino,
     *    - 3 for neutrino and equivalent neutrino
     */

    csl::InitSanitizer<int> selfConjugate{"selfConjugate"};
    csl::InitSanitizer<int> fermion{"fermion"};
    csl::InitSanitizer<int> EM_coupled{"EM_coupled"};
    csl::InitSanitizer<int> neut_coupled{"neut_coupled"};
    csl::InitSanitizer<int> neuteq_coupled{"neuteq_coupled"};
    csl::InitSanitizer<real_t> chi2{"chi2"};
    csl::InitSanitizer<int> nobs{"nobs"};
    csl::InitSanitizer<real_t> fierz{"fierz"}; //!< Fierz interference term from LQ sector
    csl::InitSanitizer<real_t> B_chi{
        "B_chi"}; //!< branching ratio of WIMP DM of mass m_p < m_chi < m_n to explain tau_n anomaly
    csl::InitSanitizer<real_t> rhob0{"rhob0"}; //!< current baryon density
    csl::InitSanitizer<real_t> b_cdm_ratio{
        "b_cdm_ratio"}; //!< current ratio of baryon density to cold dark matter density

    //  Methods for initialisation or costruction

    /**
      @brief Sets the model for the QCD equation of state.
      @param x Integer representing the QCD equation of state model.
    */
    void setQCDeosModel(const int x);

    /**
     * @brief Initialises the dataT member with the QCD EOS model corresponding to the value in model_eff.
     *        If no value is stored, the model B (model_eff=2) will be adopted
     */
    void init_heffgeff();

    /**
     * @brief Initialises the dataT member with the QCD EOS model corresponding to the value in x.
     *        It also updates the value stored in the model_eff member.
     * @param x
     */
    void init_heffgeff(const int x);

    /**
     * @brief This method initialises an instance of this object with the mass of the lightest DM candidate equal to
     * mrelic and the QCD EOS model corresponding to x, according to the enumeration QCDeosModel.
     *
     * @param mrelic Mass of the lightest DM candidate.
     * @param x Value corresponging to the QCD EOS model in the enumeration QCDeosModel.
     */
    void init(const real_t& mrelic, const int x);

    /**
     * @brief Construct a new Relicparam_t object with the specified parameters.
     *
     * @param mrelic Mass of the lightest DM candidate. If it is not provided, no value will be assigned
     *               to the corresponding member of the instance.
     * @param x Value corresponging to the QCD EOS model in the enumeration QCDeosModel. If it is not provided,
     *          the default value is 2.
     */
    Relicparam_t(const real_t& mrelic = -1., const int x = 2);

    /**
     * @brief Virtual destructor for the Relicparam_t object
     */
    virtual ~Relicparam_t(){};

    //  Metod to print the elements
    void print_relicparam(std::ostream& out = std::cout) const;
    virtual void print(std::ostream& out = std::cout) const { print_relicparam(out); };

    //  The following methods need to be called after the constructor
    //  Their purpose is to initialise this class with modified cosmological model
    //  according to the same convention of SuperIso Relic v4

    /**
     * @brief Initializes the cosmological model parameters for relic density computation.
     *
     * This procedure defines the cosmological model based on which the relic density is
     * computed. It must be called FIRST when instantiating a Relicparam_t structure.
     *
     * @param eta The baryon-to-photon ratio.
     * @param Nnu_local The number of SM neutrinos.
     * @param dNnu_local The number of extra neutrino species.
     * @param life_neutron_local The neutron lifetime.
     * @param life_neutron_error_local The error on the neutron lifetime.
     * @param xinu1_local The degeneracy parameter for electron and electron neutrino.
     * @param xinu2_local The degeneracy parameter for muon and muon neutrino.
     * @param xinu3_local The degeneracy parameter for tau and tau neutrino.
     */
    void Init_cosmomodel_param(const real_t& eta, const real_t& Nnu_local, const real_t& dNnu_local,
                               const real_t& life_neutron_local, const real_t& life_neutron_error_local,
                               const real_t& xinu1_local, const real_t& xinu2_local, const real_t& xinu3_local);

    /**
     * @brief Initializes the parameters for an included light WIMP.
     *
     * This method modifies the parameters of an included light WIMP in the cosmological model.
     * It sets various properties of the WIMP, including its mass, coupling strength, and interaction types.
     *
     * @param mass_wimp The mass of the WIMP in GeV.
     * @param EM_coupled_local Flag indicating if the WIMP couples to the electromagnetic sector (1 for yes, 0 for
     * no).
     * @param neut_coupled_local Flag indicating if the WIMP couples to neutrinos (1 for yes, 0 for no).
     * @param neuteq_coupled_local Flag indicating if the WIMP couples to neutrino-equivalent particles (1 for yes,
     * 0 for no).
     * @param fermion_local Flag indicating if the WIMP is a fermion (1 for yes, 0 for no).
     * @param selfConjugate_local Flag indicating if the WIMP is self-conjugate (1 for yes, 0 for no).
     * @param g_chi_local The coupling strength of the WIMP to Standard Model particles.
     */
    void Init_wimp(const real_t& mass_wimp, const int EM_coupled_local, const int neut_coupled_local,
                   const int neuteq_coupled_local, const int fermion_local, const int selfConjugate_local,
                   const real_t& g_chi_local);

    /**
      @brief   This procedure defines the cosmological model based on which the relic density is
      computed.

      If it is not called, no additional density will be added, and the calculation will be performed in
      the standard cosmological model.

      It adds a dark energy density according to the equation:
      \f[ \rho_D = \kappa_\rho \cdot \rho_{rad}(T_{BBN}) \cdot \left( \frac{T}{T_{BBN}} \right)^{n_\rho} \f]
      @param dd0_local \f$\kappa_\rho\f$: ratio of dark energy density to radiation energy density at \f$T_{BBN}\f$
      @param ndd_local \f$n_\rho\f$: the decreasing exponent of dark energy density
      @param T_end Reheating temperature
      @note This equation is referenced as Eq. (33) in the SuperIso Relic 3.1 Manual.
      @see SuperIso Relic Manual: https://arxiv.org/abs/1811.12736
    */
    void Init_dark_density(const real_t& dd0_local, const real_t& ndd_local, const real_t& T_end);

    /**
      @brief This procedure defines the cosmological model based on which the relic density is
      computed.

      If it is not called, no additional density will be added, and the calculation will be performed in
      the standard cosmological model.

      The dark energy density is computed as:
      \f[ \rho_D(T) = \rho_{rad}(T) \cdot \frac{g_{eff}(T_{ddeq})}{g_{eff}(T)} \cdot
      \left(\frac{h_{eff}(T)}{h_{eff}(T_{ddeq})}\right)^{n_{dd}/3} \cdot \left(\frac{T}{T_{ddeq}}\right)^{n_{dd}}
      \f]
      @param ndd_local \f$n_{dd}\f$: the decreasing exponent of dark energy density
      @param Tddeq_local \f$T_{ddeq}\f$: temperature at which \f$\rho_D = \rho_{rad}\f$ (radiation energy density)
      @param T_end Reheating temperature
      @note This parameterization provides an alternative to the one used in Init_dark_density().
      @see SuperIso Relic v4 Manual: https://arxiv.org/abs/1806.11489
    */
    void Init_dark_density2(const real_t& ndd_local, const real_t& Tddeq_local, const real_t& T_end);

    void Init_dark_coupling(const int coupD);

    /**
    @brief Initializes the quintessence model parameters for the cosmological calculation.

      This method defines a cosmological model for relic density computation that includes
      quintessence in addition to the standard cosmological scenario.

      If not called, the
      calculation uses the standard cosmological model without additional components.
      @param T12 Temperature parameter T12 for the quintessence model
      @param n2 Power law index n2 for the quintessence model
      @param T23 Temperature parameter T23 for the quintessence model
      @param n3 Power law index n3 for the quintessence model
      @param T34 Temperature parameter T34 for the quintessence model
      @param n4 Power law index n4 for the quintessence model
      @note The explanation of these parameters can be found in paragraph 5.3.2 of G. Robbin's PhD thesis.
      @see G. Robbin's PhD thesis: https://theses.hal.science/tel-02002073
    */
    void Init_quintessence(const real_t& T12, const real_t& n2, const real_t& T23, const real_t& n3, const real_t& T34,
                           const real_t& n4);

    /**
      @brief Initializes the dark entropy parameters for the cosmological model.

      This method defines the cosmological model used for relic density computation.

      If not called, the calculation uses the standard cosmological model without additional entropy.

      It adds a dark entropy component according to the equation:
      \f[ s_D = \kappa_s \cdot s_{rad}(T_{BBN}) \cdot \left( \frac{T}{T_{BBN}} \right)^{n_s} \f]
      corresponding to the "dark entropy production" scenario.

      @param sd0_local \f$\kappa_s\f$: ratio of dark entropy density to radiation entropy density at \f$T_{BBN}\f$
      @param nsd_local \f$n_s\f$: the decreasing exponent of dark entropy density
      @param T_end Reheating temperature
      @note This equation is referenced as Eq. (34) in the SuperIso Relic 3.1 Manual.
      @see SuperIso Relic Manual: https://arxiv.org/abs/1811.12736
    */
    void Init_dark_entropy(const real_t& sd0, const real_t& nsd, const real_t& T_end);

    /**
      @brief Initializes the dark entropy injection parameters for the cosmological model.

      This method defines the cosmological model used for relic density computation.

      If not called, the calculation uses the standard cosmological model without additional entropy injection.

      It adds a dark entropy injection component according to the equation:
      \f[ \Sigma_D = \kappa_\Sigma \cdot \Sigma_{rad}(T_{BBN}) \cdot \left( \frac{T}{T_{BBN}} \right)^{n_\Sigma} \f]

      @param Sigmad0_local \f$\kappa_\Sigma\f$: ratio of dark entropy injection to radiation entropy injection at
      \f$T_{BBN}\f$
      @param nSigmad_local \f$n_\Sigma\f$: the decreasing exponent of dark entropy injection
      @param T_end Reheating temperature
      @note This equation is referenced as Eq. (35) in the SuperIso Relic 3.1 Manual and Eq. (A8) in the SuperIso
      Relic 4 Manual.
      @see SuperIso Relic 3.1 Manual: https://arxiv.org/abs/1811.12736
      @see SuperIso Relic 4 Manual: https://arxiv.org/abs/1806.11489
    */
    void Init_dark_entropySigmaD(const real_t& Sigmad0, const real_t& nSigmad, const real_t& T_end);

    /**
      @brief Initializes the standard entropy injection parameters for the cosmological model.

      This method defines the cosmological model used for relic density computation.

      If not called, the calculation uses the standard cosmological model without additional entropy injection.

      It adds a standard entropy injection component according to the equation:
      \f[ \Sigma_{rad} = \kappa_{\Sigma_{rad}} \cdot \Sigma_{rad}^{eff}(T_{BBN}) \cdot \left( \frac{T}{T_{BBN}}
      \right)^{n_{\Sigma_{rad}}} \f] This modifies the radiation entropy density \f$s_{rad}\f$ following: \f[
      \dot{s}_{rad} = -3 H s_{rad} + \Sigma_{rad}(T) \f] or, equivalently: \f[ \dot{s}_{rad} = -3 H ( 1 -
      \tilde\Sigma_{rad}^\ast) s_{rad}\f]
      @param Sigmarad0_local \f$\kappa_{\Sigma_{rad}}\f$: ratio of standard entropy injection to effective standard
      entropy injection at \f$T_{BBN}\f$
      @param nSigmarad_local \f$n_{\Sigma_{rad}}\f$: the decreasing exponent of standard entropy injection
      @param T_end Reheating temperature
      @note This parameterization is based on equations (A10) and (A11) in the SuperIso Relic 4 Manual.
      @see SuperIso Relic 4 Manual: https://arxiv.org/abs/1806.11489
    */
    void Init_entropySigmarad(const real_t& Sigmarad0, const real_t& nSigmarad, const real_t& T_end);

    void Init_nonthermal(const real_t& nt0, const real_t& nnt, const real_t& T_end);
    void Init_gravitino(const real_t& mgravitino);
    void Init_scalarfield(const real_t& rhotilde_phi_Tmax, const real_t& Tmax, const real_t& T_RH,
                          const real_t& eta_phi, const real_t& n_phi);

    /**
      @brief Initializes a dark energy density table for the cosmological model.

      This method resets all parameters to those of the standard cosmological model,
      and then copies the input table of temperature and dark energy density pairs
      into the table_rhoPD data member.
      @param table[2][NTABMAX] A 2D array containing pairs of (Temperature, rho_D) values.
      @param nlines The number of rows in the input table.
      @note This method overrides any previously set dark energy density parameterizations, to allow the user to use
            a custom definition.
    */
    void Init_dark_density_table(real_t table[2][NTABMAX], const size_t nlines);

    void Init_neutron_decay(const real_t& tau, const real_t& tau_err, const real_t& fierz, const real_t& m_chi,
                            const real_t& B_chi);


    //  Methods to get thermodynamical quantities

    /**
     * @brief Computes \f$h_{eff}(T)\f$
     *
     * @param T: the input temperature
     * @return real_t \f$h_{eff}(T)\f$
     */
    real_t getheff(const real_t& T);

    /**
     * @brief Computes \f$\sqrt{g_{eff}^\ast(T)}\f$
     *
     * @param T: the input temperature
     * @return real_t \f$\sqrt{g_{eff}^\ast(T)}\f$
     */
    real_t getsqrtgstar(const real_t& T);

    /**
     * @brief Computes \f$g_{eff}(T)\f$
     *
     * @param T: the input temperature
     * @return real_t \f$g_{eff}(T)\f$
     */
    real_t getgeff(const real_t& T);

    /**
     * @brief Calls getheff(T)
     */
    inline real_t heff(const real_t& T) { return getheff(T); };

    /**
     * @brief Calls getsqrtgstar(T)
     */
    inline real_t sqrtgstar(const real_t& T) { return getsqrtgstar(T); };

    /**
     * @brief Calls getgeff(T)
     */
    inline real_t geff(const real_t& T) { return getgeff(T); };

    /**
        @brief This function computes the dark density in a modified cosmological scenario,
        which is not the decaying scalar field scenario.
        If the scenario is the standard cosmological model, or the decaying scalar
        field, this function returns 0.

        @param T The temperature at which to compute the dark energy density.
        @return The dark energy density at the given temperature. Returns 0 for the
                standard cosmological model or the decaying scalar field scenario.
    */
    real_t dark_density(const real_t& T);

    real_t dark_density_pressure(const real_t& T);

private:
    /**
     * @brief Function that appears in the denominator of the integral in dark_entropy.
     *
     * \f$ \sigma_s(T)=\exp\left( - \displaystyle \int_0^T \mathrm{d}\log(\tau) \frac{d\tilde\Sigma^\ast}{dT}(\tau)
     * \frac{\log(h_{eff}(\tau)(\tau)^3)}{(1 - \tilde \Sigma^\ast(\tau) )^2} \right) \f$
     *
     * @param T input temperature
     * @return real_t value of the function
     * @note It is equal to 1 in case of no radiation entropy injection.
     */
    real_t sigma_entropy(const real_t& T);

public:
    /**
     * @brief This function computes the total dark entropy \f$s_D(T)\f$ in the defined cosmological scenario.
     *        The scenario can be:
     *        - pure entropy production
     *        - combination of "standard" and "dark" entropy injection. In this second case, the function returns
     *          \f$s_D = \displaystyle M_{Pl}\sqrt{\frac{45}{4\pi^3}}(h_{eff}(T)T^3)^\frac{1}{1-\tilde\Sigma^\ast(T)} \sigma_s(T)
     *          \int_0^T d\ln(\tau) \frac{\sqrt{g_\ast}(\tau)\Sigma_D(\tau)}{\tilde H(\tau) (1- \tilde\Sigma^\ast(\tau))(h_{eff}(\tau)\tau^3)^\frac{2-\tilde\Sigma^\ast(\tau)}{1-\tilde\Sigma^\ast(\tau)}\sigma_s(\tau)}
     *          \f$,
     * 
     *          where \f$\sigma_s(\tau)\f$ is defined in sigma_entropy.
     *
     * @param T The temperature
     * @return real_t the total dark entropy \f$s_D(T)\f$
     * @note There might be a bug in the integral, and we need to check the formula for the combination.
     */
    real_t dark_entropy(const real_t& T);

    /**
     * @brief This function computes the derivative total dark entropy in the defined cosmological scenario.
     *
     * The scenario can be:
     * 1. pure entropy production
     * 2. combination of "standard" and "dark" entropy injection.
     * @param T The temperature
     * @return real_t the derivative total dark entropy
     */
    real_t dark_entropy_derivative(const real_t& T);

    /**
      @brief Computes the dark entropy injection \f$\Sigma_D\f$ at a given temperature.

      This function calculates the value of \f$\Sigma_D\f$ in the defined cosmological scenario.
      It returns a non-null value if the scenario is
       1. A combination of "dark entropy production" and "standard entropy injection", with a generalised version of
      equation (38) of the SuperIso Relic Manual
       2. Pure "dark entropy injection scenario", with equation (A8) of the SuperIso Relic 4 Manual

      @param T The temperature at which to compute \f$\Sigma_D\f$.
      @return The value of \f$\Sigma_D\f$ at the given temperature.
      @see SuperIso Relic Manual: https://arxiv.org/abs/1811.12736
      @see SuperIso Relic 4 Manual: https://arxiv.org/abs/1806.11489
      @note Need to check the result of point 1.
    */
    real_t dark_entropy_Sigmad(const real_t& T);

    // real_t dark_entropy_Sigmarad(const real_t &T);

    /**
      @brief Computes the standard entropy injection \f$\Sigma_{rad}\f$ at a given temperature.
      This function calculates the value of \f$\Sigma_{rad}\f$ as defined in Eq. (A11) of
      the SuperIso Relic 4 Manual. It corresponds to the scenario of "standard entropy injection".
      @param T The temperature at which to compute \f$\Sigma_{rad}\f$.
      @return The value of \f$\Sigma_{rad}\f$ at the given temperature.
      @see SuperIso Relic 4 Manual: https://arxiv.org/abs/1806.11489

    */
    real_t entropy_Sigmarad(const real_t& T);

    real_t nonthermal(const real_t& T);

    /**
      @brief Computes the neutrino density, including effects from neutrino degeneracy.
      This function calculates the neutrino density at a given neutrino temperature,
      taking into account any potential effects from neutrino degeneracy.
      @param Tnu The neutrino temperature.
      @return The neutrino density at the given temperature.
    */
    real_t neutdens(const real_t& Tnu);

    /**
      @brief Computes the temperature derivative of the neutrino energy density.
      This function calculates the rate of change of the neutrino energy density
      with respect to the neutrino temperature.
      @param Tnu The neutrino temperature.
      @return The derivative of the neutrino energy density with respect to temperature.
    */
    real_t neutdens_deriv(const real_t& Tnu);
  };

  /**
      @brief Computes the round N(z) function for neutrinos.

      This function calculates the round N(z) function for neutrinos as described in
      Pisanti et al., equation (A24).
      @param T The temperature.
      @return The value of the round N(z) function at the given temperature.
      @see Pisanti et al. (2007): https://arxiv.org/abs/0705.0290
  */
  real_t neutN(const real_t& T);

} //  end of namespace __SPEC_LIB_NAME__


/**
  @brief Stream insertion operator for DataWithTemperature_t objects.
  This operator allows DataWithTemperature_t objects to be easily output using C++ streams.
  @param out The output stream to write to.
  @param var The DataWithTemperature_t object to output.
  @return Reference to the output stream.
*/
inline std::ostream& operator<<(std::ostream& out, const __SPEC_LIB_NAME__::DataWithTemperature_t& var)
{
  var.print(out);
  return out;
}

/**
  @brief Stream insertion operator for Relicparam_t objects.
  This operator allows Relicparam_t objects to be easily output using C++ streams.
  @param out The output stream to write to.
  @param var The Relicparam_t object to output.
  @return Reference to the output stream.
*/
inline std::ostream& operator<<(std::ostream& out, const __SPEC_LIB_NAME__::Relicparam_t& var)
{
  var.print_relicparam(out);
  return out;
}
