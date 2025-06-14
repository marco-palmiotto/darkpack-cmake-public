inline const int SIZEPHYSICALSM = 20;
inline const int SIZEPHYSICALBSM = 28;
inline const int TOTAL_PARTICLES = 48;
inline const std::array<int, SIZEPHYSICALBSM> bsm_particles{
    sc_L,    su_L, sc_R, su_R, sd_R, ss_R,   sd_L,   ss_L, se_R, smu_R, se_L, smu_L, snu_e, snu_mu,
    snu_tau, st_1, st_2, sb_1, sb_2, stau_1, stau_2, sG,   N_1,  N_2,   N_3,  N_4,   C_1,   C_2};
inline const std::array<int, SIZEPHYSICALSM> sm_particles{G,      W, A,  Z,   A0, Hp, h, H0, nu_e, nu_mu,
                                                          nu_tau, e, mu, tau, u,  d,  s, c,  t,    b};
