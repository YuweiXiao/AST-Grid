#pragma once
#include <Clara/clara.hpp>

#include "general.h"
using namespace clara;

namespace Omni {

struct sim_args {
    bool showGUI{true};
    int num_frame{600};
    bool saveObj{false};      // for 3d fluid obj output
    bool dumpDensity{false};  // for 3d smoke density dump
    bool fixedTimeStep{false};
    bool saveStatus{false};
    bool openTilt{false};
    bool doViscosity{false};
    Real height{0.0};
    int substep{1};
    int resx{128};
    int resy{128};
    int resz{128};
    int dynamicEType{0};
    Real fixedAngle{-1};
};

sim_args getSimArgs(int argc, char *argv[]) {
    auto args = sim_args{};
    auto parser = Opt(args.showGUI, "gui")["-g"]["--gui"]("gui") |
                  Opt(args.num_frame, "num_frame")["-nf"]["--num_frame"]("num_frame") |
                  Opt(args.saveObj, "save_obj")["-so"]["--save_obj"]("save_obj") |
                  Opt(args.dumpDensity, "dump_density")["-dn"]["--dump_density"]("dump_density") |
                  Opt(args.fixedTimeStep, "fixed_timestep")["-t"]["--fixed_timestep"]("fixed_timestep") |
                  Opt(args.saveStatus, "save_status")["-st"]["--save_status"]("save_status") |
                  Opt(args.openTilt, "open_tilt")["-ot"]["--open_tilt"]("open_tilt") |
                  Opt(args.height, "height")["-h"]["--height"]("height") |
                  Opt(args.doViscosity, "viscosity")["-v"]["--viscosity"]("viscosity") |
                  Opt(args.substep, "substep")["-h"]["--substep"]("substep") |
                  Opt(args.resx, "resx")["-rx"]["--resx"]("resx") |
                  Opt(args.resy, "resy")["-ry"]["--resy"]("resy") |
                  Opt(args.resz, "resz")["-rz"]["--resz"]("resz") |
                  Opt(args.dynamicEType, "dynamicEType")["-de"]["--dynamicE"]("dynamicE") |
                  Opt(args.fixedAngle, "fixedAngle")["-fa"]["--angle"]("angle");
    auto const result = parser.parse(clara::Args(argc, argv));

    if (!result) {
        std::cerr << "Error in command line: " << result.errorMessage() << std::endl;
        exit(1);
    }
    return args;
}

}  // namespace Omni