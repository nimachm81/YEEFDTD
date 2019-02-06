
#include "WedgeGridArrayManipulator.h"

void test_run_fdtd_2d_pec_wedge_from_json() {
    FPNumber y0 = -3.0;
    FPNumber y1 = 3.0;
    FPNumber z0 = -3.0;
    FPNumber z1 = 3.0;
    std::size_t ny = 600;
    std::size_t nz = 600;
    FPNumber dy = (y1 - y0)/(FPNumber)(ny);
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber stabilityFactor = 0.99;
    FPNumber dt = (FPNumber)1.0/std::sqrt((FPNumber)1.0/(dy*dy) + (FPNumber)1.0/(dz*dz))*stabilityFactor;

    const FPNumber eps_r = 1.0;     // json file should be updated for eps_r != 1

    FPNumber z_j = -1.5;
    std::size_t indzJ = std::round(std::real((z_j - z0)/dz));
    FPNumber j_center_y = (y0 + y1)/(FPNumber)2.0;
    FPNumber j_decay_rate_y = 1.0;
    FPNumber j_center_t = 1.5;
    FPNumber jm_center_t = j_center_t + dt/(FPNumber)2.0*std::sqrt(eps_r);
    FPNumber jm_amplitude = -std::sqrt(eps_r);
    FPNumber j_mod_freq = 0.3;
    FPNumber j_mod_phase = M_PI/2.0;


    FPNumber wedgeAngle = 10.0/180.0*M_PI;
    FPNumber wedgeTipRadius = 0.01;
    FPNumber wedgeHeight = 5.0;
    std::array<FPNumber, 3> wedgeTipPosition = WedgeGridArrayManipulator::GetTipPositionGivenRoundedTipPosition(
                                                    wedgeAngle,
                                                    wedgeTipRadius,
                                                    {0.0, 0.0, 0.0});

    std::size_t numOfTimeSamples = 501;

    std::unordered_map<std::string, std::string> str_replacewith{
            {"\"_y0_\"", boost::lexical_cast<std::string>(std::real(y0))},
            {"\"_y1_\"", boost::lexical_cast<std::string>(std::real(y1))},
            {"\"_z0_\"", boost::lexical_cast<std::string>(std::real(z0))},
            {"\"_z1_\"", boost::lexical_cast<std::string>(std::real(z1))},
            {"\"_ny_\"", boost::lexical_cast<std::string>(ny)},
            {"\"_nz_\"", boost::lexical_cast<std::string>(nz)},
            {"\"_ny_p1_\"", boost::lexical_cast<std::string>(ny + 1)},
            {"\"_nz_p1_\"", boost::lexical_cast<std::string>(nz + 1)},
            {"\"_dy_\"", boost::lexical_cast<std::string>(std::real(dy))},
            {"\"_dz_\"", boost::lexical_cast<std::string>(std::real(dz))},
            {"\"_dt_\"", boost::lexical_cast<std::string>(std::real(dt))},
            {"\"_dt_dy_\"", boost::lexical_cast<std::string>(std::real(dt/dy))},
            {"\"_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-dt/dy))},
            {"\"_dt_dz_\"", boost::lexical_cast<std::string>(std::real(dt/dz))},
            {"\"_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-dt/dz))},
            {"\"_m_dt_dydz_\"", boost::lexical_cast<std::string>(std::real(-dt/(dy*dz)))},
            {"\"_indzJ_\"", boost::lexical_cast<std::string>(indzJ)},
            {"\"_indzJ_p1_\"", boost::lexical_cast<std::string>(indzJ + 1)},
            {"\"_indzJ_m1_\"", boost::lexical_cast<std::string>(indzJ - 1)},
            {"\"_J_center_y_\"", boost::lexical_cast<std::string>(std::real(j_center_y))},
            {"\"_J_decayrate_y_\"", boost::lexical_cast<std::string>(std::real(j_decay_rate_y))},
            {"\"_J_center_t_\"", boost::lexical_cast<std::string>(std::real(j_center_t))},
            {"\"_Jm_center_t_\"", boost::lexical_cast<std::string>(std::real(jm_center_t))},
            {"\"_Jm_amp_\"", boost::lexical_cast<std::string>(std::real(jm_amplitude))},
            {"\"_J_mod_freq_\"", boost::lexical_cast<std::string>(std::real(j_mod_freq))},
            {"\"_J_mod_phase_\"", boost::lexical_cast<std::string>(std::real(j_mod_phase))},
            {"\"_wedge_angle_\"", boost::lexical_cast<std::string>(std::real(wedgeAngle))},
            {"\"_wedge_tip_radius_\"", boost::lexical_cast<std::string>(std::real(wedgeTipRadius))},
            {"\"_wedge_height_\"", boost::lexical_cast<std::string>(std::real(wedgeHeight))},
            {"\"_wedge_tip_y_\"", boost::lexical_cast<std::string>(std::real(wedgeTipPosition[1]))},
            {"\"_wedge_tip_z\"", boost::lexical_cast<std::string>(std::real(wedgeTipPosition[2]))},
            {"\"_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee2D_PEC_Wedge.json",
                "instructions/processed/MaxwellYee2D_PEC_Wedge_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/processed/MaxwellYee2D_PEC_Wedge_processed.json");
    fileTranslator.Translate();
}


