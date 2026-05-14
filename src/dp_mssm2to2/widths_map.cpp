#include "correspondance.hpp"
namespace mssm2to2::corr
{
  const std::unordered_map<std::string, Entry_t> widths_map{{{
                                                                 EMPTYCHAR + corr::h,
                                                             },
                                                             {&c_width_h, 1, 1}},
                                                            {{
                                                                 EMPTYCHAR + corr::H0,
                                                             },
                                                             {&c_width_H0, 1, 1}},
                                                            {{
                                                                 EMPTYCHAR + corr::Hp,
                                                             },
                                                             {&c_width_Hp, 1, 1}},
                                                            {{
                                                                 EMPTYCHAR + corr::A0,
                                                             },
                                                             {&c_width_A0, 1, 1}}};
}
